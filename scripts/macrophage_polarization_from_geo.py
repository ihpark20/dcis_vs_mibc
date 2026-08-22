"""Macrophage polarization in DCIS and microinvasive cores.

Macrophages sit on a continuum between an inflammatory M1 state and a
tissue-remodelling, immunosuppressive M2 one, and which way a tumor's macrophages lean
says something about the environment it is building. The split comes from the cell states
themselves: a macrophage is M1 or M2 because the subclustering named its subcluster
`macrophage_M1` or `macrophage_M2`, scoring these programmes against the rest of the myeloid
pool's:

    M1   CD86, CD80, IL1B, TNF
    M2   CD163, MRC1, MARCO, ARG1, IL10, TREM2

Macrophages the subclustering placed elsewhere — proliferating, monocyte — carry no
polarization call: they count towards `n_mac` but towards neither side. Both programme
scores are still written per cell, for anyone who wants the continuum rather than the
split.

Each core then contributes the share of its macrophages on each side, their ratio, and how
many of each it holds — per thousand cells, and per square millimetre of tissue, tissue
area being what the core QC measured as actually covered by cells. DCIS is compared with
mDCIS by Mann-Whitney over cores, with Benjamini-Hochberg across the measures.

Usage:
    python scripts/use_published_subclusters.py
    python scripts/macrophage_polarization_from_geo.py

Output:
    03.data_processed/macrophage/polarization_cells.csv    per macrophage
    03.data_processed/macrophage/polarization_percore.csv  per core
    03.data_processed/macrophage/polarization_stats.csv    DCIS vs mDCIS
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_clusters.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
CORE_QC = ROOT / "02.tma_core_qc/core_qc_raw.csv"
OUT_DIR = ROOT / "03.data_processed/macrophage"

M1_MARKERS = ["CD86", "CD80", "IL1B", "TNF"]
M2_MARKERS = ["CD163", "MRC1", "MARCO", "ARG1", "IL10", "TREM2"]
MEASURES = [
    "M1_pct_of_mac",
    "M2_pct_of_mac",
    "M1_M2_ratio",
    "M1_per1k",
    "M2_per1k",
    "M1_per_mm2",
    "M2_per_mm2",
    "n_mac",
]
SEED = 42


def polarize(a, cell_ids):
    """Score each macrophage for both programmes and take the higher one."""
    b = a[a.obs_names.isin(cell_ids)].copy()
    b.X = b.layers["lognorm"].copy() if "lognorm" in b.layers else b.X
    for name, markers in (("M1", M1_MARKERS), ("M2", M2_MARKERS)):
        present = [g for g in markers if g in b.var_names]
        sc.tl.score_genes(b, present, score_name=f"score_{name}", random_state=SEED)
        print(f"  {name}: {len(present)}/{len(markers)} markers on the panel — {present}")
    b.obs["polarization"] = b.obs["score_M1"] - b.obs["score_M2"]

    return b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--core-qc", default=CORE_QC, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    args = ap.parse_args()

    major = pd.read_csv(args.major)  # cell_id, pool, leiden_sub, subtype_merged, major_celltype
    macrophage_rows = major[major["major_celltype"] == "Macrophage_Mono"]
    macrophages = set(macrophage_rows["cell_id"])
    state_of = dict(
        zip(
            macrophage_rows["cell_id"],
            macrophage_rows["subtype_merged"].map({"macrophage_M1": "M1", "macrophage_M2": "M2"}),
            strict=True,
        )
    )
    print(f"{len(macrophages):,} macrophage/monocyte cells")
    print(macrophage_rows["subtype_merged"].value_counts().to_string())

    a = sc.read_h5ad(args.clustered)
    b = polarize(a, macrophages)

    obs = b.obs.copy()
    obs["cell_id"] = b.obs_names
    obs["state"] = obs["cell_id"].map(state_of)
    obs["leiden_sub"] = obs["cell_id"].map(
        dict(zip(macrophage_rows["cell_id"], macrophage_rows["leiden_sub"], strict=True))
    )
    obs["core"] = obs["slide"].astype(str) + "_" + obs["core_id"].astype(str)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    obs[
        [
            "cell_id",
            "core",
            "sample_group",
            "leiden_sub",
            "score_M1",
            "score_M2",
            "polarization",
            "state",
        ]
    ].to_csv(args.out_dir / "polarization_cells.csv", index=False)
    print(
        f"\n{(obs.state == 'M1').sum():,} M1 | {(obs.state == 'M2').sum():,} M2 | "
        f"{obs.state.isna().sum():,} without a polarization call"
    )

    # cells per core, and the tissue area the QC measured
    everything = sc.read_h5ad(args.clustered, backed="r").obs
    core_of = everything["slide"].astype(str) + "_" + everything["core_id"].astype(str)
    n_cells = core_of.value_counts()
    qc = pd.read_csv(args.core_qc)
    qc["core"] = qc["slide"] + "_" + qc["core_id"]
    tissue_mm2 = (qc.set_index("core")["tissue_area_um2"] / 1e6).to_dict()

    rows = []
    for core, g in obs.groupby("core", observed=True):
        n_m1 = int((g["state"] == "M1").sum())
        n_m2 = int((g["state"] == "M2").sum())
        polarized = n_m1 + n_m2
        if polarized == 0:
            continue
        total = int(n_cells.get(core, len(g)))
        area = tissue_mm2.get(core, np.nan)
        rows.append(
            {
                "core": core,
                "sample_group": g["sample_group"].iloc[0],
                "n_mac": len(g),
                "n_M1": n_m1,
                "n_M2": n_m2,
                "n_polarized": polarized,
                "M1_pct_of_mac": round(n_m1 / polarized * 100, 2),
                "M2_pct_of_mac": round(n_m2 / polarized * 100, 2),
                "M1_M2_ratio": round(n_m1 / max(n_m2, 1), 4),
                "M1_per1k": round(n_m1 / total * 1000, 3),
                "M2_per1k": round(n_m2 / total * 1000, 3),
                "M1_per_mm2": round(n_m1 / area, 1) if area == area else np.nan,
                "M2_per_mm2": round(n_m2 / area, 1) if area == area else np.nan,
            }
        )
    per_core = pd.DataFrame(rows)
    per_core.to_csv(args.out_dir / "polarization_percore.csv", index=False)

    stats = []
    for measure in MEASURES:
        dcis = per_core.loc[per_core["sample_group"] == "dcis", measure].dropna()
        mibc = per_core.loc[per_core["sample_group"] == "mibc", measure].dropna()
        stats.append(
            {
                "measure": measure,
                "dcis_median": round(dcis.median(), 3),
                "mibc_median": round(mibc.median(), 3),
                "direction": "mDCIS higher" if mibc.median() > dcis.median() else "DCIS higher",
                "p": mannwhitneyu(dcis, mibc).pvalue,
                "n_dcis": len(dcis),
                "n_mibc": len(mibc),
            }
        )
    stats = pd.DataFrame(stats)
    stats["fdr"] = multipletests(stats["p"], method="fdr_bh")[1]
    stats.to_csv(args.out_dir / "polarization_stats.csv", index=False)

    print(f"\nDCIS {stats.n_dcis.iloc[0]} cores vs mDCIS {stats.n_mibc.iloc[0]} cores")
    print(stats.round(4).to_string(index=False))
    print(f"\nwrote {args.out_dir}")


if __name__ == "__main__":
    main()
