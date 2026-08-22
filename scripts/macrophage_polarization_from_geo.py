"""Macrophage polarization in DCIS and microinvasive cores.

Macrophages sit on a continuum between an inflammatory M1 state and a
tissue-remodelling, immunosuppressive M2 one, and which way a tumor's macrophages lean
says something about the environment it is building. Each macrophage is scored for both
programmes and takes the higher one:

    M1   CD86, CD80, IL1B, TNF
    M2   CD163, MRC1, MARCO, ARG1, IL10, TREM2

Polarization is settled per subcluster, as the annotation does: a subcluster is a group of
cells the clustering found coherent, and averaging over it is steadier than scoring a single
cell with a handful of sparse markers. The two programme scores are compared with each other
only — the pool's other programmes decide what a subcluster is, this decides which way its
macrophages lean. Per-cell scores stay in the output for anyone who wants the continuum
rather than the split.

Each core then contributes the share of its macrophages on each side, their ratio, and how
many of each it holds — per thousand cells, and per square millimetre of tissue, tissue
area being what the core QC measured as actually covered by cells. DCIS is compared with
mDCIS by Mann-Whitney over cores, with Benjamini-Hochberg across the measures.

Usage:
    python scripts/major_celltypes_from_geo.py
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
MYELOID_POOL = ROOT / "03.data_processed/subclustered/myeloid.h5ad"
CLUSTERED = ROOT / "03.data_processed/integrated_qc_passed_from_geo.h5ad"
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

    by_sub = b.obs.groupby("leiden_sub", observed=True)[["score_M1", "score_M2"]].mean()
    by_sub["state"] = np.where(by_sub["score_M1"] > by_sub["score_M2"], "M1", "M2")
    by_sub["n_cells"] = b.obs["leiden_sub"].value_counts()
    print("\nsubcluster polarization")
    print(by_sub.round(3).sort_values("n_cells", ascending=False).to_string())
    b.obs["state"] = b.obs["leiden_sub"].map(by_sub["state"]).astype(str)
    return b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--myeloid-pool", default=MYELOID_POOL, type=Path)
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--core-qc", default=CORE_QC, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    args = ap.parse_args()

    major = pd.read_csv(args.major)
    macrophages = set(major.loc[major["major_celltype"] == "Macrophage_Mono", "cell_id"])
    print(f"{len(macrophages):,} macrophage/monocyte cells")

    a = sc.read_h5ad(args.myeloid_pool)
    b = polarize(a, macrophages)

    obs = b.obs.copy()
    obs["cell_id"] = b.obs_names
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
    print(f"\n{(obs.state == 'M1').sum():,} M1 | {(obs.state == 'M2').sum():,} M2")

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
        total = int(n_cells.get(core, len(g)))
        area = tissue_mm2.get(core, np.nan)
        rows.append(
            {
                "core": core,
                "sample_group": g["sample_group"].iloc[0],
                "n_mac": len(g),
                "n_M1": n_m1,
                "n_M2": n_m2,
                "M1_pct_of_mac": round(n_m1 / len(g) * 100, 2),
                "M2_pct_of_mac": round(n_m2 / len(g) * 100, 2),
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
