"""Where PD-1+ T cells meet PD-L1, relative to the tumor.

A checkpoint interaction needs the two cells to be close enough to touch, and where that
happens matters: a PD-1+ T cell engaged inside a tumor mass or at its edge is being held
off at the point of invasion, one engaged out in the stroma is not.

Each T cell is placed in one of four classes:

    engaged-inside      PD-1+, a PD-L1+ cell within 30 um, and at least 60 % of its own
                        neighbours epithelial
    engaged-boundary    the same, with 20-60 % epithelial neighbours
    engaged-outside     the same, below 20 %
    PD1pos_unengaged    PD-1+ with no PD-L1+ cell within 30 um

Cores are compared on the share of their PD-1+ T cells engaged inside or at the boundary,
and only cores with at least 20 PD-1+ T cells take part: below that the share moves in
steps too coarse to compare — a core with five PD-1+ cells can only return 0, 20, 40 %.
Pooling the cells of all cores instead would make the comparison significant on the
strength of a few large cores, so the core is the unit.

Usage:
    python scripts/major_celltypes_from_geo.py
    python scripts/checkpoint_engagement_from_geo.py

Output:
    03.data_processed/ccc/checkpoint_tcells.csv     counts per core and class
    03.data_processed/ccc/checkpoint_percore.csv    per core, with the share engaged
    03.data_processed/ccc/checkpoint_stats.csv      DCIS vs mDCIS
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from analysis_exclusions import EXCLUDED_CORES
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_qc_passed_from_geo.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
OUT_DIR = ROOT / "03.data_processed/ccc"

LIGAND, RECEPTOR = "CD274", "PDCD1"  # PD-L1 -> PD-1
RADIUS = 30.0
EPITHELIAL = {"Tumor", "Myoepithelial"}
INSIDE, BOUNDARY = 0.6, 0.2
MIN_PD1 = 20


def region_of(fraction_epithelial):
    if fraction_epithelial >= INSIDE:
        return "inside"
    return "boundary" if fraction_epithelial >= BOUNDARY else "outside"


def classify_core(pos, is_epithelial, is_tcell, ligand, receptor):
    """Sort the T cells of one core into the four classes."""
    tree = cKDTree(pos)
    receptor_idx = np.where(receptor)[0]

    region = {}
    if len(receptor_idx):
        for i, neighbours in zip(
            receptor_idx, tree.query_ball_point(pos[receptor_idx], RADIUS), strict=True
        ):
            own = 1 if is_epithelial[i] else 0
            share = (is_epithelial[neighbours].sum() - own) / max(len(neighbours) - 1, 1)
            region[i] = region_of(share)

    engaged = np.zeros(len(pos), dtype=bool)
    if ligand.sum() and len(receptor_idx):
        ligand_pos = pos[ligand]
        ligand_tree = cKDTree(ligand_pos)
        for i in receptor_idx:
            near = ligand_tree.query_ball_point(pos[i], RADIUS)
            near = [k for k in near if not np.allclose(ligand_pos[k], pos[i], atol=1e-6)]
            engaged[i] = len(near) > 0

    classes = []
    for i in np.where(is_tcell)[0]:
        if receptor[i] and engaged[i]:
            classes.append(f"engaged-{region[i]}")
        elif receptor[i]:
            classes.append("PD1pos_unengaged")
        else:
            classes.append("PD1neg")
    return classes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    ap.add_argument("--min-pd1", default=MIN_PD1, type=int)
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered)
    obs = a.obs.copy()
    obs["cell_id"] = a.obs_names
    obs["core"] = obs["slide"].astype(str) + "_" + obs["core_id"].astype(str)
    major = pd.read_csv(args.major, usecols=["cell_id", "major_celltype"])
    obs = obs.merge(major, on="cell_id", how="left")
    obs["major_celltype"] = obs["major_celltype"].fillna("unassigned")
    keep = ~obs["core"].isin(EXCLUDED_CORES)

    counts = a.layers["counts"] if "counts" in a.layers else a.X
    counts = sp.csc_matrix(counts) if not sp.issparse(counts) else counts.tocsc()
    gene = {g: i for i, g in enumerate(a.var_names)}
    ligand_all = np.asarray(counts[:, gene[LIGAND]].todense()).ravel() > 0
    receptor_all = np.asarray(counts[:, gene[RECEPTOR]].todense()).ravel() > 0

    obs = obs[keep].reset_index(drop=True)
    ligand_all, receptor_all = ligand_all[keep.to_numpy()], receptor_all[keep.to_numpy()]
    print(f"{len(obs):,} cells in {obs['core'].nunique()} cores")

    rows = []
    for core, idx in obs.groupby("core", observed=True).indices.items():
        sub = obs.iloc[idx]
        classes = classify_core(
            sub[["x_centroid", "y_centroid"]].to_numpy(),
            sub["major_celltype"].isin(EPITHELIAL).to_numpy(),
            (sub["major_celltype"] == "T_NK").to_numpy(),
            ligand_all[idx],
            receptor_all[idx],
        )
        for cls in classes:
            rows.append({"core": core, "sample_group": sub["sample_group"].iloc[0], "class": cls})

    cells = pd.DataFrame(rows)
    per_class = cells.groupby(["core", "sample_group", "class"]).size().rename("n").reset_index()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    per_class.to_csv(args.out_dir / "checkpoint_tcells.csv", index=False)

    wide = per_class.pivot_table(
        index=["core", "sample_group"], columns="class", values="n", fill_value=0
    ).reset_index()
    for column in ("engaged-inside", "engaged-boundary", "engaged-outside", "PD1pos_unengaged"):
        if column not in wide:
            wide[column] = 0
    wide["n_pd1_pos"] = (
        wide["engaged-inside"]
        + wide["engaged-boundary"]
        + wide["engaged-outside"]
        + wide["PD1pos_unengaged"]
    )
    wide["pct_engaged_inside_boundary"] = (
        (wide["engaged-inside"] + wide["engaged-boundary"]) / wide["n_pd1_pos"] * 100
    ).round(2)
    wide.to_csv(args.out_dir / "checkpoint_percore.csv", index=False)

    tested = wide[wide["n_pd1_pos"] >= args.min_pd1]
    dcis = tested.loc[tested["sample_group"] == "dcis", "pct_engaged_inside_boundary"]
    mibc = tested.loc[tested["sample_group"] == "mibc", "pct_engaged_inside_boundary"]
    p = mannwhitneyu(dcis, mibc).pvalue
    stats = pd.DataFrame(
        [
            {
                "measure": "pct of PD-1+ T cells engaged inside or at the boundary",
                "min_pd1_pos": args.min_pd1,
                "n_dcis_cores": len(dcis),
                "n_mibc_cores": len(mibc),
                "dcis_median": round(dcis.median(), 2),
                "mibc_median": round(mibc.median(), 2),
                "direction": "mDCIS higher" if mibc.median() > dcis.median() else "DCIS higher",
                "p": p,
            }
        ]
    )
    stats.to_csv(args.out_dir / "checkpoint_stats.csv", index=False)

    pooled = cells[cells["class"] != "PD1neg"]
    share = pooled.groupby("sample_group")["class"].value_counts(normalize=True).unstack() * 100
    print("\npooled over all cores, % of PD-1+ T cells")
    print(share.round(1).to_string())
    print(
        f"\ncores with at least {args.min_pd1} PD-1+ T cells: "
        f"DCIS {len(dcis)}, mDCIS {len(mibc)} (of {len(wide)})"
    )
    print(stats.round(4).to_string(index=False))
    print(f"\nwrote {args.out_dir}")


if __name__ == "__main__":
    main()
