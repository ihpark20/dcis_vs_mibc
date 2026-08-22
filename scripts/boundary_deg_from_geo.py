"""Differential expression at the tumor boundary, within each coverage category.

Boundary cells of the three coverage classes are not the same population, so comparing DCIS
with mDCIS across all of them at once would mix a change in expression with a change in what
the boundary is made of. Each class is therefore tested on its own: the cells of a class are
summed per core, and cores are the replicates of a DESeq2 model.

The comparison is then run twice. Once as `~ pathology`, and once with the local tumor
density around each boundary cell added as a covariate, `~ tumor_z + pathology`. That
covariate is the point of the exercise: transcripts leak from neighbouring cells (see the
spillover section), so a gene that looks differential at a boundary may only be reporting how
much tumor sits next to it. A hit that survives the adjustment is not explained by tumor
proximity; one that does not, is.

Tumor density is the count of tumor cells within 30 um of a boundary cell, per 1000 um2 of
that disc, with the cell itself discounted where it is a tumor cell. A core contributes a
class when it holds at least 20 cells of it, and a class is tested when at least three cores
per group remain.

Usage:
    python scripts/boundary_myoep_from_geo.py
    python scripts/boundary_deg_from_geo.py

Output:
    03.data_processed/boundary/boundary_deg_by_category.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_clusters.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
BOUNDARY_DIR = ROOT / "03.data_processed/boundary"

RADIUS = 30.0
DISC_K = np.pi * RADIUS**2 / 1000.0  # 1000 um2 units of the 30 um disc
CATEGORIES = ["myoep-sheath", "myoep-lined", "myoep-deficient"]
MIN_CELLS_PER_CORE = 20
MIN_CORES_PER_GROUP = 3
MIN_GENE_COUNT = 10


def tumor_density(boundary, tumor_positions):
    """Tumor cells within 30 um of each boundary cell, per 1000 um2."""
    density = np.full(len(boundary), np.nan)
    for core, idx in boundary.groupby("core", observed=True).indices.items():
        sub = boundary.iloc[idx]
        positions = tumor_positions.get(core)
        if positions is None or not len(positions):
            density[idx] = 0.0
            continue
        tree = cKDTree(positions)
        near = np.array(
            [
                len(x)
                for x in tree.query_ball_point(sub[["x_centroid", "y_centroid"]].to_numpy(), RADIUS)
            ],
            dtype=float,
        )
        # a lined or deficient boundary cell is itself a tumor cell; do not count it
        itself = (sub["boundary_type"] != "myoep-sheath").to_numpy().astype(float)
        density[idx] = np.maximum(near - itself, 0) / DISC_K
    return density


def deseq(counts, metadata, design):
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        ref_level=["sample_group", "dcis"],
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["sample_group", "mibc", "dcis"], quiet=True)
    stats.summary()
    return stats.results_df[["log2FoldChange", "pvalue", "padj"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--boundary-dir", default=BOUNDARY_DIR, type=Path)
    args = ap.parse_args()

    boundary = pd.read_csv(args.boundary_dir / "boundary_myoep_cells.csv")
    a = sc.read_h5ad(args.clustered)
    obs = a.obs.copy()
    obs["cell_id"] = a.obs_names
    obs["core"] = obs["slide"].astype(str) + "_" + obs["core_id"].astype(str)
    major = pd.read_csv(args.major, usecols=["cell_id", "major_celltype"])
    obs = obs.merge(major, on="cell_id", how="left")
    tumor = obs[obs["major_celltype"] == "Tumor"]
    tumor_positions = {
        core: g[["x_centroid", "y_centroid"]].to_numpy()
        for core, g in tumor.groupby("core", observed=True)
    }
    boundary["tumor_density"] = tumor_density(boundary, tumor_positions)
    print(
        f"{len(boundary):,} boundary cells, tumor density "
        f"{boundary['tumor_density'].median():.2f} per 1000 um2 at the median"
    )

    a = a[boundary["cell_id"].to_numpy()].copy()
    counts_matrix = a.layers["counts"] if "counts" in a.layers else a.X
    counts_matrix = (
        sp.csr_matrix(counts_matrix) if not sp.issparse(counts_matrix) else counts_matrix.tocsr()
    )

    results = []
    for category in CATEGORIES:
        selected = (boundary["boundary_type"] == category).to_numpy()
        sub = boundary[selected]
        codes, cores = pd.factorize(sub["core"].to_numpy())
        indicator = sp.csr_matrix(
            (np.ones(len(sub)), (codes, np.arange(len(sub)))), shape=(len(cores), len(sub))
        )
        pseudobulk = np.asarray((indicator @ counts_matrix[selected]).todense())
        n_cells = np.asarray(indicator.sum(1)).ravel()

        per_core = pd.DataFrame(
            {
                "core": sub["core"].to_numpy(),
                "sample_group": sub["sample_group"].to_numpy(),
                "tumor_density": sub["tumor_density"].to_numpy(),
            }
        )
        group = per_core.groupby("core", sort=False)["sample_group"].first().reindex(cores)
        density = per_core.groupby("core", sort=False)["tumor_density"].mean().reindex(cores)

        keep = n_cells >= MIN_CELLS_PER_CORE
        counts = pd.DataFrame(
            np.rint(pseudobulk[keep]).astype(int), index=cores[keep], columns=a.var_names
        )
        metadata = pd.DataFrame(
            {"sample_group": group.to_numpy()[keep], "tumor_density": density.to_numpy()[keep]},
            index=cores[keep],
        )
        metadata["sample_group"] = pd.Categorical(
            metadata["sample_group"], categories=["dcis", "mibc"]
        )
        metadata["tumor_z"] = (
            metadata["tumor_density"] - metadata["tumor_density"].mean()
        ) / metadata["tumor_density"].std()
        n_dcis = int((metadata["sample_group"] == "dcis").sum())
        n_mibc = int((metadata["sample_group"] == "mibc").sum())
        if min(n_dcis, n_mibc) < MIN_CORES_PER_GROUP:
            print(f"  {category}: DCIS {n_dcis}, mDCIS {n_mibc} cores — too few, skipped")
            continue
        counts = counts.loc[:, counts.sum(0) >= MIN_GENE_COUNT]
        print(f"  {category}: DCIS {n_dcis} vs mDCIS {n_mibc} cores, {counts.shape[1]} genes")

        unadjusted = deseq(counts, metadata[["sample_group"]], "~sample_group")
        adjusted = deseq(
            metadata=metadata[["sample_group", "tumor_z"]],
            counts=counts,
            design="~tumor_z + sample_group",
        )
        merged = (
            unadjusted.rename(columns=lambda c: c + "_unadjusted")
            .join(adjusted.rename(columns=lambda c: c + "_adjusted"))
            .reset_index(names="gene")
        )
        merged["category"] = category
        merged["n_dcis_cores"] = n_dcis
        merged["n_mibc_cores"] = n_mibc
        merged["tumor_density_dcis"] = round(
            float(metadata.loc[metadata["sample_group"] == "dcis", "tumor_density"].median()), 3
        )
        merged["tumor_density_mibc"] = round(
            float(metadata.loc[metadata["sample_group"] == "mibc", "tumor_density"].median()), 3
        )
        results.append(merged)

    out = pd.concat(results, ignore_index=True)
    out.to_csv(args.boundary_dir / "boundary_deg_by_category.csv", index=False)

    hits = out[(out["padj_unadjusted"] < 0.05) | (out["padj_adjusted"] < 0.05)]
    print("\ngenes reaching FDR < 0.05 before or after the adjustment")
    print(
        hits[["category", "gene", "log2FoldChange_adjusted", "padj_unadjusted", "padj_adjusted"]]
        .sort_values("padj_adjusted")
        .round(4)
        .to_string(index=False)
    )
    print(f"\nwrote {args.boundary_dir / 'boundary_deg_by_category.csv'}")


if __name__ == "__main__":
    main()
