"""Assign clusters to the seven lineage pools, the way the paper did.

The paper groups its initial clusters by hand into seven pools — CL3, CL4 and CL10 make the
tumor pool, CL2, CL11 and CL12 the myeloid one, and so on — and those groupings are in
`pool_config.py`. Because `cluster_from_geo.py` has already matched every cluster to its
published counterpart, the same grouping can be applied directly here, which is what makes
the pools the paper's pools rather than something rederived.

Rederiving them is still useful as a check, so each cluster is also scored against the
lineage genes of all seven pools and the best-scoring pool recorded beside the assigned
one. The two agreeing is evidence that the cluster matching held; they disagree when a
cluster failed to match cleanly, and the script says so instead of letting it pass.

Usage:
    python scripts/cluster_from_geo.py
    python scripts/assign_pools_from_geo.py

Output:
    03.data_processed/pool_assignment.csv
"""

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
from pool_config import POOLS

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_qc_passed_from_geo.h5ad"
NAMING = ROOT / "03.data_processed/cluster_naming.csv"
OUT = ROOT / "03.data_processed/pool_assignment.csv"

SEED = 42
CLUSTER_OF_POOL = {f"CL{c}": pool for pool, cfg in POOLS.items() for c in cfg["clusters"]}


def lineage_scores(a, cluster_key):
    """Mean score of every pool's lineage genes, per cluster."""
    b = a.copy()
    b.X = b.layers["counts"].copy() if "counts" in b.layers else b.X
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    for name, pool in POOLS.items():
        genes = [g for g in pool["filter_genes"] if g in b.var_names]
        sc.tl.score_genes(b, genes, score_name=f"pool_{name}", random_state=SEED)
    cols = [f"pool_{n}" for n in POOLS]
    return b.obs.groupby(cluster_key, observed=True)[cols].mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--naming", default=NAMING, type=Path)
    ap.add_argument("--out", default=OUT, type=Path)
    ap.add_argument("--cluster-key", default="cluster", help="obs column holding the CL names")
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered)
    if args.cluster_key not in a.obs:
        raise SystemExit(
            f"{args.clustered} has no obs['{args.cluster_key}'] — rerun cluster_from_geo.py, "
            "which matches the clusters to the published numbering"
        )

    scores = lineage_scores(a, args.cluster_key)
    by_score = scores.idxmax(axis=1).str.replace("pool_", "", regex=False)

    table = pd.DataFrame(
        {
            "cluster": scores.index,
            "pool": [CLUSTER_OF_POOL.get(c, "unassigned") for c in scores.index],
            "pool_by_lineage_score": by_score.values,
            "n_cells": a.obs[args.cluster_key].value_counts().reindex(scores.index).values,
        }
    )
    table["agrees"] = table["pool"] == table["pool_by_lineage_score"]

    if args.naming.exists():
        naming = pd.read_csv(args.naming)[["cluster", "leiden", "correlation"]]
        table = table.merge(naming, on="cluster", how="left")

    table = table.sort_values("cluster", key=lambda c: c.str.removeprefix("CL").astype(int))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, index=False)
    print(table.to_string(index=False))

    unassigned = table[table["pool"] == "unassigned"]
    if len(unassigned):
        print(f"\nnot covered by any pool: {', '.join(unassigned['cluster'])}")
    disagreeing = table[~table["agrees"]]
    if len(disagreeing):
        print("\nassigned pool differs from the lineage score — check the cluster matching:")
        print(disagreeing.to_string(index=False))
    else:
        print(f"\nall {len(table)} clusters: the paper's grouping and the lineage score agree")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
