"""Stage-1 clustering, starting from the GEO submission (GSE343808).

Takes the object written by `build_anndata_from_geo.py` and runs the initial clustering
of the paper on it: the two slides are integrated with Harmony and partitioned by Leiden,
with the parameters used for the published clustering.

Two of the 41 cores are dropped before integration, as in the paper — TMA1 D8 for a median
of 0 transcripts per cell, TMA1 D4 for carrying 0.9 % tumor cells. Both verdicts are read
from `02.tma_core_qc/tma_core_qc.csv`, so the cell set matches the published one; pass
--all-cores to skip that filter and cluster every core.

Clustering is stochastic in practice: Harmony and Leiden depend on floating-point
accumulation order, so a rerun on different hardware reproduces the cluster structure
rather than each cell's label exactly.

Usage:
    python scripts/load_from_geo.py --token <reviewer token>
    python scripts/build_anndata_from_geo.py
    python scripts/cluster_from_geo.py

Reclustering takes about an hour and does not land on exactly the paper's partition; pass
--labels to skip it and take the published assignment from
`03.data_processed/integrated_cluster_labels.csv.gz`, which covers all 520,506 cells of the
analysis. Everything downstream works the same either way.

Clusters are renamed to the paper's CL0-CL14 by matching each one's mean expression to the
reference profiles in `03.data_processed/integrated_cluster_profiles.csv`, so the numbering
means the same thing here as in the paper; `03.data_processed/integrated_cluster_annotation_manuscript.csv`
carries what each CL was annotated as and the evidence behind it. Leiden's own numbering is kept in `obs["leiden"]`.

Output:
    03.data_processed/integrated_qc_passed_from_geo.h5ad
    03.data_processed/cluster_naming.csv
"""

import argparse
from pathlib import Path

import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.optimize import linear_sum_assignment

# resolved from this file, so a clone works wherever it sits
ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "03.data_processed/geo_slides.h5ad"
PROFILES = ROOT / "03.data_processed/integrated_cluster_profiles.csv"
LABELS = ROOT / "03.data_processed/integrated_cluster_labels.csv.gz"
QC = ROOT / "02.tma_core_qc/tma_core_qc.csv"
OUT = ROOT / "03.data_processed/integrated_qc_passed_from_geo.h5ad"

N_HVG = 300
N_PCS = 30
N_NEIGHBORS = 15
RESOLUTION = 0.5
SEED = 42


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=IN, type=Path)
    ap.add_argument("--out", default=OUT, type=Path)
    ap.add_argument(
        "--qc", default=QC, type=Path, help="core QC verdicts from qc_cores_composition.py"
    )
    ap.add_argument("--profiles", default=PROFILES, type=Path, help="reference CL profiles")
    ap.add_argument(
        "--labels",
        type=Path,
        help="skip clustering and take the published labels from this file",
    )
    ap.add_argument("--max-cells", type=int, default=0, help="subsample, for a smoke test")
    ap.add_argument(
        "--all-cores",
        action="store_true",
        help="cluster every core, without applying the recorded core QC",
    )
    args = ap.parse_args()

    import pandas as pd

    a = sc.read_h5ad(args.input)
    if args.all_cores:
        a = a[a.obs["core_id"].astype(str) != "unassigned"].copy()
        print(f"all cores: {a.obs['core_id'].nunique()}  ->  {a.n_obs:,} cells")
        return cluster(a, args)

    qc = pd.read_csv(args.qc)
    keep = qc[qc["analysis_include"]]
    passing = set(keep["slide"] + "_" + keep["core_id"])
    dropped = (
        qc.loc[~qc["analysis_include"], "slide"] + "_" + qc.loc[~qc["analysis_include"], "core_id"]
    ).tolist()
    key = a.obs["slide"].astype(str) + "_" + a.obs["core_id"].astype(str)
    a = a[key.isin(passing)].copy()
    print(f"cores kept: {len(passing)} (dropped {', '.join(dropped)})  ->  {a.n_obs:,} cells")

    return cluster(a, args)


def name_clusters(a, profiles_path, cluster_key="leiden"):
    """Rename Leiden clusters to the paper's CL numbers.

    Leiden numbers its clusters by size, so a rerun renumbers everything. The published
    clusters are matched by what they express instead: each new cluster is summarised by
    its mean expression, correlated against the reference profile of every CL, and the
    one-to-one assignment maximising total correlation wins. That keeps CL0-CL14 meaning
    the same population it means in the paper.
    """
    reference = pd.read_csv(profiles_path, index_col=0)
    genes = [g for g in reference.columns if g in a.var_names]

    lognorm = a.copy()
    lognorm.X = lognorm.layers["counts"].copy()
    sc.pp.normalize_total(lognorm, target_sum=1e4)
    sc.pp.log1p(lognorm)
    x = lognorm[:, genes].X
    frame = pd.DataFrame(
        x.todense() if hasattr(x, "todense") else x,
        index=a.obs[cluster_key].astype(str).values,
        columns=genes,
    )
    observed = frame.groupby(level=0).mean()

    corr = pd.DataFrame(
        np.corrcoef(observed.to_numpy(), reference[genes].to_numpy())[
            : len(observed), len(observed) :
        ],
        index=observed.index,
        columns=reference.index,
    )
    rows, cols = linear_sum_assignment(-corr.to_numpy())
    mapping = {corr.index[r]: corr.columns[c] for r, c in zip(rows, cols, strict=True)}
    table = pd.DataFrame(
        {
            "leiden": list(mapping),
            "cluster": [mapping[k] for k in mapping],
            "correlation": [round(corr.loc[k, mapping[k]], 3) for k in mapping],
            "n_cells": [int((a.obs[cluster_key].astype(str) == k).sum()) for k in mapping],
        }
    ).sort_values("cluster", key=lambda c: c.str.removeprefix("CL").astype(int))

    a.obs["cluster"] = a.obs[cluster_key].astype(str).map(mapping).astype("category")
    return table


def attach_published(a, labels_path):
    """Take the manuscript's cluster labels instead of clustering again."""
    labels = pd.read_csv(labels_path).set_index("cell_id")["cluster"]
    known = a.obs_names.intersection(labels.index)
    print(f"{len(known):,} of {a.n_obs:,} cells carry a published label")
    a = a[known].copy()
    a.obs["cluster"] = labels.reindex(a.obs_names).astype("category")
    a.obs["leiden"] = a.obs["cluster"].astype(str).str.removeprefix("CL").astype("category")
    print(a.obs["cluster"].value_counts().sort_index().to_string())
    return a


def cluster(a, args):
    if args.labels:
        a = attach_published(a, args.labels)
        a.write_h5ad(args.out, compression="gzip")
        print(f"wrote {args.out} (published clustering, nothing recomputed)")
        return None

    if args.max_cells and a.n_obs > args.max_cells:
        sc.pp.subsample(a, n_obs=args.max_cells, random_state=0)
        print(f"subsampled to {a.n_obs:,} cells")

    sc.pp.filter_cells(a, min_counts=10)
    sc.pp.filter_cells(a, min_genes=5)
    print(f"after cell filters: {a.n_obs:,} cells")

    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, batch_key="slide", flavor="seurat")
    sc.pp.scale(a, max_value=10)
    # random_state=0 is scanpy's default and what the published run used
    sc.tl.pca(a, n_comps=N_PCS, use_highly_variable=True, random_state=0)

    ho = hm.run_harmony(
        a.obsm["X_pca"], a.obs, vars_use=["slide"], random_state=SEED, max_iter_harmony=30
    )
    a.obsm["X_pca_harmony"] = np.array(ho.Z_corr)

    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep="X_pca_harmony")
    sc.tl.umap(a, random_state=SEED)
    sc.tl.leiden(a, resolution=RESOLUTION, random_state=SEED, key_added="leiden")

    if args.profiles.exists():
        table = name_clusters(a, args.profiles)
        table.to_csv(args.out.parent / "cluster_naming.csv", index=False)
        print("\nLeiden clusters matched to the published numbering:")
        print(table.to_string(index=False))
    else:
        print(f"{args.profiles} not found; clusters keep their Leiden numbers")
        a.obs["cluster"] = "CL" + a.obs["leiden"].astype(str)
    print("\n" + a.obs["cluster"].value_counts().sort_index().to_string())

    a.X = a.layers["counts"].copy()
    a.write_h5ad(args.out, compression="gzip")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
