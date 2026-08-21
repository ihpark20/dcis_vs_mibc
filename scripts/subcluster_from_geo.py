"""Pool-wise subclustering of the integrated clustering.

A single clustering of half a million cells over a 374-gene panel separates lineages but
not the states inside them. The initial clusters are therefore grouped into seven lineage
pools and each pool is re-clustered on its own markers, which is what gives the cell
states used in the paper. Per-pool settings live in `pool_config.py`, taken from the
original subcluster scripts.

For every pool:

  1. cells whose cluster belongs to the pool enter, provided they carry a non-zero count
     for at least one of the pool's lineage genes;
  2. those lineage markers become the feature space for PCA, in place of highly variable
     genes — with 374 panel genes, HVG selection inside one lineage is dominated by
     transcripts that leak in from neighbouring cells;
  3. scale, PCA, Harmony over slide, a 15-neighbour graph, UMAP, Leiden at the pool's
     resolution;
  4. each cell is scored for the pool's programmes and every subcluster takes the name of
     its highest-scoring one.

Pools are matched to clusters by lineage score, not by cluster number: a rerun numbers its
clusters differently, so each cluster is assigned to the pool whose lineage genes it
expresses most strongly.

Usage:
    python scripts/cluster_from_geo.py
    python scripts/subcluster_from_geo.py

Output:
    03.data_processed/subclustered/<pool>.h5ad
    03.data_processed/subclustered/subcluster_assignments.csv
    03.data_processed/subclustered/pool_of_cluster.csv
"""

import argparse
from pathlib import Path

import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc
from pool_config import POOLS

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "03.data_processed/integrated_qc_passed_from_geo.h5ad"
OUT_DIR = ROOT / "03.data_processed/subclustered"

N_NEIGHBORS = 15
SEED = 42


def assign_pools(a, cluster_key):
    """Give every cluster to the pool whose lineage genes it expresses most strongly."""
    b = a.copy()
    b.X = b.layers["counts"].copy() if "counts" in b.layers else b.X
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    for name, pool in POOLS.items():
        genes = [g for g in pool["filter_genes"] if g in b.var_names]
        sc.tl.score_genes(b, genes, score_name=f"pool_{name}", random_state=SEED)

    cols = [f"pool_{n}" for n in POOLS]
    means = b.obs.groupby(cluster_key, observed=True)[cols].mean()
    winner = means.idxmax(axis=1).str.replace("pool_", "", regex=False)
    table = pd.DataFrame(
        {
            "cluster": means.index,
            "pool": winner.values,
            "n_cells": a.obs[cluster_key].value_counts().reindex(means.index).values,
            "best_score": means.max(axis=1).round(3).values,
        }
    )
    return winner.to_dict(), table


def subcluster(a, name, pool):
    """The shared recipe, run on one pool."""
    sub = a.copy()
    sub.X = sub.layers["counts"].copy() if "counts" in sub.layers else sub.X

    present = [g for g in pool["filter_genes"] if g in sub.var_names]
    if present:
        counts = sub[:, present].X
        keep = np.asarray((counts > 0).sum(axis=1)).ravel() > 0
        print(f"  lineage filter on {len(present)} genes: {keep.sum():,}/{sub.n_obs:,} cells kept")
        sub = sub[keep].copy()

    sub.layers["counts"] = sub.X.copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sub.layers["lognorm"] = sub.X.copy()

    features = [g for g in pool["feature_genes"] if g in sub.var_names]
    sub.var["highly_variable"] = sub.var_names.isin(features)
    print(f"  features: {len(features)}/{len(pool['feature_genes'])}")

    sc.pp.scale(sub, max_value=10)
    n_comps = min(pool["n_pcs"], len(features) - 1)
    sc.tl.pca(sub, n_comps=n_comps, use_highly_variable=True, random_state=SEED)

    harmony = hm.run_harmony(
        sub.obsm["X_pca"], sub.obs, vars_use=["slide"], random_state=SEED, max_iter_harmony=30
    )
    sub.obsm["X_pca_harmony"] = np.array(harmony.Z_corr)

    sc.pp.neighbors(sub, n_neighbors=N_NEIGHBORS, n_pcs=n_comps, use_rep="X_pca_harmony")
    sc.tl.umap(sub, random_state=SEED)
    sc.tl.leiden(sub, resolution=pool["resolution"], random_state=SEED, key_added="leiden_sub")

    sub.X = sub.layers["lognorm"].copy()
    for programme, genes in pool["programs"].items():
        present = [g for g in genes if g in sub.var_names]
        if present:
            sc.tl.score_genes(sub, present, score_name=f"score_{programme}", random_state=SEED)

    score_cols = [c for c in sub.obs.columns if c.startswith("score_")]
    means = sub.obs.groupby("leiden_sub", observed=True)[score_cols].mean()
    dominant = means.idxmax(axis=1).str.replace("score_", "", regex=False)
    sub.obs["program"] = sub.obs["leiden_sub"].map(dominant).astype(str)
    sub.X = sub.layers["counts"].copy()

    print(f"  {sub.n_obs:,} cells -> {sub.obs['leiden_sub'].nunique()} subclusters")
    print(sub.obs["program"].value_counts().to_string())
    return sub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=IN, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    ap.add_argument("--cluster-key", default="leiden")
    ap.add_argument("--pools", nargs="*", default=list(POOLS), help="subset of pools to run")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    a = sc.read_h5ad(args.input)
    print(f"{a.n_obs:,} cells x {a.n_vars} genes, {a.obs[args.cluster_key].nunique()} clusters")

    pool_of, table = assign_pools(a, args.cluster_key)
    table.to_csv(args.out_dir / "pool_of_cluster.csv", index=False)
    print("\ncluster -> pool")
    print(table.to_string(index=False))

    assignments = []
    for name in args.pools:
        clusters = [c for c, p in pool_of.items() if p == name]
        print(f"\n=== {name} ({POOLS[name]['label']}) from clusters {sorted(clusters)} ===")
        if not clusters:
            print("  no cluster assigned, skipped")
            continue
        pool_cells = a[a.obs[args.cluster_key].isin(clusters)].copy()
        sub = subcluster(pool_cells, name, POOLS[name])
        sub.write_h5ad(args.out_dir / f"{name}.h5ad", compression="gzip")
        assignments.append(
            pd.DataFrame(
                {
                    "cell_id": sub.obs_names,
                    "pool": name,
                    "leiden_sub": sub.obs["leiden_sub"].astype(str).values,
                    "program": sub.obs["program"].values,
                }
            )
        )

    out = pd.concat(assignments, ignore_index=True)
    out.to_csv(args.out_dir / "subcluster_assignments.csv", index=False)
    print(f"\n{len(out):,} cells assigned across {out['pool'].nunique()} pools")
    print(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
