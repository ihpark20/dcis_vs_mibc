"""Differential expression within each major cell type, DCIS versus microinvasive.

The compartment comparison asks whether cancer, stroma or immune cells differ as a whole;
this asks it of each cell type on its own, which is where a change confined to one lineage
shows up. The cost is thirteen comparisons instead of three, on smaller populations, so
fewer genes clear correction — the two are complementary rather than one replacing the
other.

Counts are summed per core within a cell type and cores are the replicates of a DESeq2
model, `~ sample_group`, contrasting microinvasive against DCIS. A core contributes a cell
type when it holds at least 20 of its cells, and a cell type is tested when at least three
cores per group remain — which is what excludes the rarer types rather than any judgement
about them.

As in the compartment comparison, every result carries the spillover index of its gene.
A microinvasive core gives every non-epithelial cell more tumor to sit beside, so a gene
leaking out of tumor rises in those cells whether or not they express it; `spillover_from_
cancer` and `top_spillover_source` are there to be read before such a hit is believed.

Usage:
    python scripts/spillover_from_geo.py
    python scripts/celltype_deg_from_geo.py

Output:
    03.data_processed/compartment/celltype_deg.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from analysis_exclusions import EXCLUDED_CORES
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_clusters.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
SPILLOVER = ROOT / "03.data_processed/spillover/spillover_by_source.csv"
OUT_DIR = ROOT / "03.data_processed/compartment"

CANCER_TYPES = ["Tumor", "Myoepithelial"]
DROP_TYPES = {"Neutrophil"}
MIN_CELLS_PER_CORE = 20
MIN_CORES_PER_GROUP = 3
MIN_GENE_COUNT = 10


def deseq(counts, metadata):
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~sample_group",
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
    ap.add_argument("--spillover", default=SPILLOVER, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered)
    obs = a.obs.copy()
    obs["cell_id"] = a.obs_names
    obs["core"] = obs["slide"].astype(str) + "_" + obs["core_id"].astype(str)
    major = pd.read_csv(args.major, usecols=["cell_id", "major_celltype"])
    obs = obs.merge(major, on="cell_id", how="left")
    keep = (
        ~obs["core"].isin(EXCLUDED_CORES)
        & obs["major_celltype"].notna()
        & ~obs["major_celltype"].isin(DROP_TYPES)
    ).to_numpy()

    counts_matrix = a.layers["counts"] if "counts" in a.layers else a.X
    counts_matrix = (
        sp.csr_matrix(counts_matrix) if not sp.issparse(counts_matrix) else counts_matrix.tocsr()
    )
    counts_matrix, obs = counts_matrix[keep], obs[keep].reset_index(drop=True)

    spillover = pd.read_csv(args.spillover)
    cancer_leak = (
        spillover[spillover["source"].isin(CANCER_TYPES)].groupby("gene")["spillover_index"].max()
    )
    top_source = spillover.loc[spillover.groupby("gene")["spillover_index"].idxmax()].set_index(
        "gene"
    )

    results, skipped = [], []
    for celltype in sorted(obs["major_celltype"].unique()):
        selected = (obs["major_celltype"] == celltype).to_numpy()
        sub = obs[selected]
        codes, cores = pd.factorize(sub["core"].to_numpy())
        indicator = sp.csr_matrix(
            (np.ones(len(sub)), (codes, np.arange(len(sub)))), shape=(len(cores), len(sub))
        )
        pseudobulk = np.asarray((indicator @ counts_matrix[selected]).todense())
        n_cells = np.asarray(indicator.sum(1)).ravel()
        group = sub.groupby("core", sort=False)["sample_group"].first().reindex(cores)

        enough = n_cells >= MIN_CELLS_PER_CORE
        counts = pd.DataFrame(
            np.rint(pseudobulk[enough]).astype(int), index=cores[enough], columns=a.var_names
        )
        metadata = pd.DataFrame({"sample_group": group.to_numpy()[enough]}, index=cores[enough])
        metadata["sample_group"] = pd.Categorical(
            metadata["sample_group"], categories=["dcis", "mibc"]
        )
        n_dcis = int((metadata["sample_group"] == "dcis").sum())
        n_mibc = int((metadata["sample_group"] == "mibc").sum())
        if min(n_dcis, n_mibc) < MIN_CORES_PER_GROUP:
            skipped.append(f"{celltype} (DCIS {n_dcis}, mDCIS {n_mibc} cores)")
            continue
        counts = counts.loc[:, counts.sum(0) >= MIN_GENE_COUNT]
        print(
            f"  {celltype:16s} {selected.sum():>7,} cells | DCIS {n_dcis:2d} vs mDCIS {n_mibc:2d} cores"
            f" | {counts.shape[1]} genes"
        )

        result = deseq(counts, metadata).reset_index(names="gene")
        result["celltype"] = celltype
        result["n_dcis_cores"] = n_dcis
        result["n_mibc_cores"] = n_mibc
        result["spillover_from_cancer"] = result["gene"].map(cancer_leak).fillna(0).round(3)
        result["top_spillover_source"] = result["gene"].map(top_source["source"])
        result["top_spillover_index"] = result["gene"].map(top_source["spillover_index"]).fillna(0)
        results.append(result)

    out = pd.concat(results, ignore_index=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_dir / "celltype_deg.csv", index=False)

    if skipped:
        print("\ntoo few cores to test: " + ", ".join(skipped))
    hits = out[out["padj"] < 0.05].sort_values(["celltype", "padj"])
    print(f"\n{len(hits)} genes at FDR < 0.05 across {out['celltype'].nunique()} cell types")
    print(
        hits[
            [
                "celltype",
                "gene",
                "log2FoldChange",
                "padj",
                "spillover_from_cancer",
                "top_spillover_source",
            ]
        ]
        .round(3)
        .to_string(index=False)
    )
    own = hits[hits["top_spillover_source"] == hits["celltype"]]
    print(
        f"\nof these, {len(own)} leak most strongly out of the very cell type they were "
        "found in, and so are not explained by spillover from elsewhere"
    )
    print(f"\nwrote {args.out_dir / 'celltype_deg.csv'}")


if __name__ == "__main__":
    main()
