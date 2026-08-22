"""What the microenvironment is made of, and how densely it is packed.

Two questions about the same cells, answered separately because they can disagree. The
first is compositional: among microenvironment cells, does the mix of major types differ
between DCIS and microinvasive cores? The second is about crowding: is a given cell type
packed more densely in the tissue, regardless of how the rest of the mix moves?

Tumor and myoepithelial cells are left out and the composition renormalised within the ten
immune and stromal types, so the question is about the microenvironment rather than about
how much tumor a core happens to contain. Neutrophils are out as well: the panel cannot
resolve them.

Composition is compositional data — the proportions of a core sum to one, so a rise in any
one type forces the others down and testing them as independent percentages invents
differences. Counts are therefore centre-log-ratio transformed against the geometric mean
of the ten (with a pseudocount of 0.5 for the zeros), which frees them from the constraint.
The groups are compared as a whole by PERMANOVA on the Aitchison distance, then type by
type on the CLR coordinates with Benjamini-Hochberg correction.

Density is cells of a type per square millimetre of microenvironment tissue, the tissue
being the summed area of the microenvironment cells themselves rather than the area of the
core, so that a core with more open stroma is not counted as sparser.

Usage:
    python scripts/use_published_subclusters.py
    python scripts/tme_composition_from_geo.py

Output:
    03.data_processed/tme/composition_global.csv   PERMANOVA over the ten types
    03.data_processed/tme/composition_pertype.csv  CLR medians, per type
    03.data_processed/tme/composition_percore.csv
    03.data_processed/tme/density_pertype.csv      cells per mm2 of microenvironment
    03.data_processed/tme/density_cellstate.csv    the same, by cell state
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from analysis_exclusions import EXCLUDED_CORES
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_clusters.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
OUT_DIR = ROOT / "03.data_processed/tme"

TME_MAJORS = [
    "T_NK",
    "B_cell",
    "Plasma",
    "Macrophage_Mono",
    "Dendritic",
    "Mast",
    "CAF_Fibroblast",
    "Endothelial",
    "Pericyte",
    "Adipocyte",
]
PSEUDOCOUNT = 0.5
MIN_TME_CELLS = 50
N_PERMUTATIONS = 999
SEED = 0


def clr(counts):
    """Centre-log-ratio: each count against the geometric mean of its core."""
    p = counts + PSEUDOCOUNT
    p = p / p.sum(axis=1, keepdims=True)
    log_p = np.log(p)
    return log_p - log_p.mean(axis=1, keepdims=True)


def permanova(coords, groups, rng, n_permutations=N_PERMUTATIONS):
    """PERMANOVA on Euclidean (= Aitchison) distance between CLR coordinates."""
    n = len(coords)
    squared = ((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
    total = squared.sum() / (2 * n)

    def within(labels):
        return sum(
            squared[np.ix_(np.where(labels == g)[0], np.where(labels == g)[0])].sum()
            / (2 * (labels == g).sum())
            for g in np.unique(labels)
        )

    k = len(np.unique(groups))
    w = within(groups)
    f = ((total - w) / (k - 1)) / (w / (n - k))
    exceeded = 1
    for _ in range(n_permutations):
        wp = within(rng.permutation(groups))
        exceeded += (((total - wp) / (k - 1)) / (wp / (n - k))) >= f
    return float(f), exceeded / (n_permutations + 1)


def compare(frame, value_column, group_column="sample_group"):
    """Mann-Whitney per group of the frame, with BH correction across them."""
    rows = []
    for name, g in frame.groupby(level=0) if isinstance(frame, pd.Series) else frame:
        dcis = g.loc[g[group_column] == "dcis", value_column]
        mibc = g.loc[g[group_column] == "mibc", value_column]
        rows.append(
            {
                "name": name,
                "dcis_median": round(dcis.median(), 3),
                "mibc_median": round(mibc.median(), 3),
                "delta": round(mibc.median() - dcis.median(), 3),
                "p": mannwhitneyu(dcis, mibc).pvalue,
            }
        )
    out = pd.DataFrame(rows)
    out["fdr"] = multipletests(out["p"], method="fdr_bh")[1]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered, backed="r")
    cells = a.obs[["slide", "core_id", "sample_group", "cell_area"]].copy()
    cells["cell_id"] = a.obs_names
    cells["core"] = cells["slide"].astype(str) + "_" + cells["core_id"].astype(str)
    major = pd.read_csv(args.major, usecols=["cell_id", "major_celltype", "subtype_merged"])
    cells = cells.merge(major, on="cell_id", how="inner")
    cells = cells[~cells["core"].isin(EXCLUDED_CORES)]
    tme = cells[cells["major_celltype"].isin(TME_MAJORS)].copy()
    print(f"{len(tme):,} microenvironment cells in {tme['core'].nunique()} cores")

    counts = (
        tme.groupby(["core", "major_celltype"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=TME_MAJORS, fill_value=0)
    )
    counts = counts[counts.sum(axis=1) >= MIN_TME_CELLS]
    group_of = tme.groupby("core")["sample_group"].first().reindex(counts.index)
    coords = clr(counts.to_numpy(float))

    f, p = permanova(
        coords, (group_of.to_numpy() == "mibc").astype(int), np.random.default_rng(SEED)
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "n_types": len(TME_MAJORS),
                "n_dcis_cores": int((group_of == "dcis").sum()),
                "n_mibc_cores": int((group_of == "mibc").sum()),
                "permanova_F": round(f, 3),
                "permanova_p": p,
            }
        ]
    ).to_csv(args.out_dir / "composition_global.csv", index=False)
    print(f"\nPERMANOVA over the ten types: F = {f:.3f}, p = {p:.3f}")

    long = pd.DataFrame(coords, index=counts.index, columns=TME_MAJORS).stack().rename("clr")
    long = long.reset_index().rename(columns={"level_1": "major_celltype"})
    long["sample_group"] = long["core"].map(group_of)
    long["pct_of_tme"] = (counts.div(counts.sum(axis=1), axis=0).stack().to_numpy() * 100).round(3)
    long.to_csv(args.out_dir / "composition_percore.csv", index=False)

    per_type = compare(long.groupby("major_celltype"), "clr").rename(
        columns={"name": "major_celltype"}
    )
    per_type.to_csv(args.out_dir / "composition_pertype.csv", index=False)
    print(per_type.round(4).to_string(index=False))

    tissue_mm2 = tme.groupby("core")["cell_area"].sum() / 1e6
    for level, column, out_name in (
        ("major_celltype", "major_celltype", "density_pertype.csv"),
        ("cell state", "subtype_merged", "density_cellstate.csv"),
    ):
        counted = tme.groupby(["core", column]).size().unstack(fill_value=0)
        density = counted.div(tissue_mm2, axis=0)
        frame = density.stack().rename("per_mm2").reset_index()
        frame.columns = ["core", column, "per_mm2"]
        frame["sample_group"] = frame["core"].map(group_of)
        if column == "subtype_merged":
            frame["major_celltype"] = frame[column].map(
                tme.drop_duplicates(column).set_index(column)["major_celltype"]
            )
        stats = compare(frame.groupby(column), "per_mm2").rename(columns={"name": column})
        stats.to_csv(args.out_dir / out_name, index=False)
        if column == "major_celltype":
            print(f"\ndensity per mm2 of microenvironment tissue, by {level}")
            print(stats.round(3).to_string(index=False))
        else:
            dendritic = tme.loc[tme["major_celltype"] == "Dendritic", column].unique()
            print("\ndendritic cell states")
            print(stats[stats[column].isin(dendritic)].round(3).to_string(index=False))

    print(f"\nwrote {args.out_dir}")


if __name__ == "__main__":
    main()
