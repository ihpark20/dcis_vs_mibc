"""Differential expression within each compartment, DCIS versus microinvasive.

Cells are grouped into the three compartments the tissue is built from and each is compared
on its own, so that a difference in expression is not confused with a difference in what a
core is made of — a microinvasive core holds more tumor, and testing everything together
would report that as a change in expression.

    cancer   Tumor, Myoepithelial
    stroma   CAF_Fibroblast, Endothelial, Pericyte, Adipocyte
    immune   T_NK, B_cell, Plasma, Macrophage_Mono, Dendritic, Mast

Counts are summed per core within a compartment and cores are the replicates of a DESeq2
model, `~ sample_group`, contrasting microinvasive against DCIS. A core contributes a
compartment when it holds at least 20 of its cells, and a compartment is tested when at
least three cores per group remain.

Every result carries the spillover index of its gene alongside it, from
`spillover/spillover_by_source.csv`. Segmentation assigns a transcript to whichever cell it
falls in, so an epithelial gene can be counted in the fibroblast next to a tumor cell — and
a microinvasive core has more tumor for stromal and immune cells to sit beside. A gene that
comes out higher in mDCIS stroma while leaking strongly from tumor is not evidence of
stromal expression; the columns `spillover_from_cancer` and `top_spillover_source` are there
to be read before any such hit is believed.

Usage:
    python scripts/spillover_from_geo.py
    python scripts/compartment_deg_from_geo.py

Output:
    03.data_processed/compartment/compartment_deg.csv
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

COMPARTMENTS = {
    "cancer": ["Tumor", "Myoepithelial"],
    "stroma": ["CAF_Fibroblast", "Endothelial", "Pericyte", "Adipocyte"],
    "immune": ["T_NK", "B_cell", "Plasma", "Macrophage_Mono", "Dendritic", "Mast"],
}
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
    keep = (~obs["core"].isin(EXCLUDED_CORES)).to_numpy()

    counts_matrix = a.layers["counts"] if "counts" in a.layers else a.X
    counts_matrix = (
        sp.csr_matrix(counts_matrix) if not sp.issparse(counts_matrix) else counts_matrix.tocsr()
    )
    counts_matrix, obs = counts_matrix[keep], obs[keep].reset_index(drop=True)

    spillover = pd.read_csv(args.spillover)
    cancer_leak = (
        spillover[spillover["source"].isin(COMPARTMENTS["cancer"])]
        .groupby("gene")["spillover_index"]
        .max()
    )
    top_source = spillover.loc[spillover.groupby("gene")["spillover_index"].idxmax()].set_index(
        "gene"
    )

    results = []
    for compartment, types in COMPARTMENTS.items():
        selected = obs["major_celltype"].isin(types).to_numpy()
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
            print(f"  {compartment}: DCIS {n_dcis}, mDCIS {n_mibc} cores — too few, skipped")
            continue
        counts = counts.loc[:, counts.sum(0) >= MIN_GENE_COUNT]
        print(
            f"  {compartment}: {selected.sum():,} cells, DCIS {n_dcis} vs mDCIS {n_mibc} cores, "
            f"{counts.shape[1]} genes"
        )

        result = deseq(counts, metadata).reset_index(names="gene")
        result["compartment"] = compartment
        result["n_dcis_cores"] = n_dcis
        result["n_mibc_cores"] = n_mibc
        result["spillover_from_cancer"] = result["gene"].map(cancer_leak).fillna(0).round(3)
        result["top_spillover_source"] = result["gene"].map(top_source["source"])
        result["top_spillover_index"] = result["gene"].map(top_source["spillover_index"]).fillna(0)
        results.append(result)

    out = pd.concat(results, ignore_index=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_dir / "compartment_deg.csv", index=False)

    hits = out[out["padj"] < 0.05].sort_values(["compartment", "padj"])
    print(f"\n{len(hits)} genes at FDR < 0.05")
    print(
        hits[
            [
                "compartment",
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
    leaky = hits[(hits["compartment"] != "cancer") & (hits["spillover_from_cancer"] > 1)]
    if len(leaky):
        print(
            f"\n{len(leaky)} of them leak from tumor or myoepithelium (index > 1) — "
            "read those against the spillover table before believing them"
        )
    print(f"\nwrote {args.out_dir / 'compartment_deg.csv'}")


if __name__ == "__main__":
    main()
