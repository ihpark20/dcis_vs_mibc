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

Output:
    03.data_processed/integrated_qc_passed_from_geo.h5ad
"""

import argparse
from pathlib import Path

import harmonypy as hm
import numpy as np
import scanpy as sc

# resolved from this file, so a clone works wherever it sits
ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "03.data_processed/geo_slides.h5ad"
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


def cluster(a, args):
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
    sc.tl.pca(a, n_comps=N_PCS, use_highly_variable=True, random_state=SEED)

    ho = hm.run_harmony(
        a.obsm["X_pca"], a.obs, vars_use=["slide"], random_state=SEED, max_iter_harmony=30
    )
    a.obsm["X_pca_harmony"] = np.array(ho.Z_corr)

    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep="X_pca_harmony")
    sc.tl.umap(a, random_state=SEED)
    sc.tl.leiden(a, resolution=RESOLUTION, random_state=SEED, key_added="leiden")
    print(a.obs["leiden"].value_counts().sort_index().to_string())

    a.X = a.layers["counts"].copy()
    a.write_h5ad(args.out, compression="gzip")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
