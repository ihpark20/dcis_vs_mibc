"""Myoepithelial coverage of the tumor boundary, and how it differs between DCIS and mDCIS.

A duct in situ is wrapped in a myoepithelial sheath; losing that wrapping is what lets
tumor cells meet the stroma directly, and microinvasion is the point at which they do.
Each boundary cell is therefore sorted by what covers it:

    myoep-sheath      the boundary cell is itself myoepithelial
    myoep-lined       a tumor cell with a myoepithelial cell within 15 um
    myoep-deficient   a tumor cell facing the microenvironment with no myoepithelium

Cores are the unit of comparison, not cells: each core contributes the percentage of its
boundary in each class, and DCIS is compared with mDCIS by Mann-Whitney with
Benjamini-Hochberg correction over the four measures.

Shares alone can mislead — a core with twice the tumor can hold twice the deficient
boundary at the same percentage — so the counts are reported too, both absolute and per 100
tumor cells, which is what asks whether a lesion carries more exposed boundary for the
amount of tumor it has.

Usage:
    python scripts/boundary_from_geo.py
    python scripts/boundary_myoep_from_geo.py

Output:
    03.data_processed/boundary/boundary_myoep_cells.csv     per boundary cell
    03.data_processed/boundary/boundary_myoep_percore.csv   per core
    03.data_processed/boundary/boundary_myoep_overall.csv
    03.data_processed/boundary/boundary_myoep_stats.csv     DCIS vs mDCIS, shares
    03.data_processed/boundary/boundary_counts_percore.csv  counts, absolute and per tumor
    03.data_processed/boundary/boundary_counts_stats.csv    DCIS vs mDCIS, counts
"""

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_clusters.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
BOUNDARY_DIR = ROOT / "03.data_processed/boundary"

MYOEP_RADIUS = 15.0
MIN_BOUNDARY_CELLS = 5
TYPES = ["myoep-sheath", "myoep-lined", "myoep-deficient"]


def classify_core(sub, boundary_ids):
    """Label each boundary cell of one core by what covers it."""
    pos = sub[["x_centroid", "y_centroid"]].to_numpy()
    major = sub["major_celltype"].to_numpy()
    is_boundary = sub["cell_id"].isin(boundary_ids).to_numpy()
    if is_boundary.sum() < MIN_BOUNDARY_CELLS:
        return None

    myoep = pos[major == "Myoepithelial"]
    tree = cKDTree(myoep) if len(myoep) else None

    labels = []
    for i in is_boundary.nonzero()[0]:
        if major[i] == "Myoepithelial":
            labels.append("myoep-sheath")
        elif tree is not None and tree.query_ball_point(pos[i], MYOEP_RADIUS):
            labels.append("myoep-lined")
        else:
            labels.append("myoep-deficient")

    idx = is_boundary.nonzero()[0]
    return pd.DataFrame(
        {
            "cell_id": sub["cell_id"].to_numpy()[idx],
            "core": sub["core"].iloc[0],
            "sample_group": sub["sample_group"].iloc[0],
            "x_centroid": pos[idx, 0],
            "y_centroid": pos[idx, 1],
            "boundary_type": labels,
        }
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--boundary-dir", default=BOUNDARY_DIR, type=Path)
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered, backed="r")
    cells = a.obs[["slide", "core_id", "sample_group", "x_centroid", "y_centroid"]].copy()
    cells["cell_id"] = a.obs_names
    cells["core"] = cells["slide"].astype(str) + "_" + cells["core_id"].astype(str)
    major = pd.read_csv(args.major, usecols=["cell_id", "major_celltype"])
    cells = cells.merge(major, on="cell_id", how="left")
    cells["major_celltype"] = cells["major_celltype"].fillna("unassigned")

    pool = pd.read_csv(args.boundary_dir / "pool_boundary_cells.csv")
    boundary_ids = set(pool.loc[pool["zone"] == "boundary", "cell_id"])
    cells = cells[cells["core"].isin(set(pool["core"]))]
    print(f"{len(boundary_ids):,} boundary cells over {cells['core'].nunique()} cores")

    frames = [
        f for _, sub in cells.groupby("core") if (f := classify_core(sub, boundary_ids)) is not None
    ]
    per_cell = pd.concat(frames, ignore_index=True)
    per_cell.to_csv(args.boundary_dir / "boundary_myoep_cells.csv", index=False)

    overall = per_cell["boundary_type"].value_counts().rename("n").reindex(TYPES).reset_index()
    overall["pct"] = (overall["n"] / overall["n"].sum() * 100).round(1)
    overall.to_csv(args.boundary_dir / "boundary_myoep_overall.csv", index=False)
    print("\n" + overall.to_string(index=False))

    rows = []
    for core, sub in per_cell.groupby("core"):
        share = sub["boundary_type"].value_counts(normalize=True) * 100
        rec = {
            "core": core,
            "sample_group": sub["sample_group"].iloc[0],
            "n_boundary": len(sub),
        }
        for t in TYPES:
            rec[t] = round(share.get(t, 0.0), 2)
        rec["lined_total"] = round(rec["myoep-sheath"] + rec["myoep-lined"], 2)
        rows.append(rec)
    per_core = pd.DataFrame(rows)
    per_core.to_csv(args.boundary_dir / "boundary_myoep_percore.csv", index=False)

    tumor_per_core = cells[cells["major_celltype"] == "Tumor"].groupby("core", observed=True).size()
    count_rows = []
    for core, sub in per_cell.groupby("core", observed=True):
        n_tumor = int(tumor_per_core.get(core, 0))
        for category in [*TYPES, "all"]:
            n = len(sub) if category == "all" else int((sub["boundary_type"] == category).sum())
            count_rows.append(
                {
                    "core": core,
                    "sample_group": sub["sample_group"].iloc[0],
                    "category": category,
                    "count": n,
                    "per100tumor": round(n / n_tumor * 100, 3) if n_tumor else float("nan"),
                    "n_tumor": n_tumor,
                }
            )
    counts = pd.DataFrame(count_rows)
    counts.to_csv(args.boundary_dir / "boundary_counts_percore.csv", index=False)

    count_stats = []
    for (category, measure), g in counts.melt(
        id_vars=["core", "sample_group", "category"],
        value_vars=["count", "per100tumor"],
        var_name="measure",
    ).groupby(["category", "measure"], observed=True):
        dcis = g.loc[g["sample_group"] == "dcis", "value"].dropna()
        mibc = g.loc[g["sample_group"] == "mibc", "value"].dropna()
        count_stats.append(
            {
                "category": category,
                "measure": measure,
                "dcis_median": round(dcis.median(), 2),
                "mibc_median": round(mibc.median(), 2),
                "direction": "mDCIS higher" if mibc.median() > dcis.median() else "DCIS higher",
                "p": mannwhitneyu(dcis, mibc).pvalue,
            }
        )
    count_stats = pd.DataFrame(count_stats)
    count_stats["fdr"] = multipletests(count_stats["p"], method="fdr_bh")[1]
    count_stats.to_csv(args.boundary_dir / "boundary_counts_stats.csv", index=False)
    print("\ncounts per core, absolute and per 100 tumor cells")
    print(count_stats.round(4).to_string(index=False))

    stats = []
    for measure in [*TYPES, "lined_total"]:
        dcis = per_core.loc[per_core["sample_group"] == "dcis", measure]
        mibc = per_core.loc[per_core["sample_group"] == "mibc", measure]
        stats.append(
            {
                "measure": measure,
                "dcis_median": round(dcis.median(), 2),
                "mibc_median": round(mibc.median(), 2),
                "direction": "mDCIS higher" if mibc.median() > dcis.median() else "DCIS higher",
                "p": mannwhitneyu(dcis, mibc).pvalue,
            }
        )
    stats = pd.DataFrame(stats)
    stats["fdr"] = multipletests(stats["p"], method="fdr_bh")[1]
    stats.to_csv(args.boundary_dir / "boundary_myoep_stats.csv", index=False)

    print(
        f"\nDCIS {(per_core.sample_group == 'dcis').sum()} cores vs "
        f"mDCIS {(per_core.sample_group == 'mibc').sum()} cores (median % of boundary)"
    )
    print(stats.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
