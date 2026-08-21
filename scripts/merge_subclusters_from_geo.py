"""Clean up and merge the pool subclusters.

Leiden at a fixed resolution splits a pool into more pieces than there are cell states.
This step merges them by programme and marks the ones too small to stand on their own:

    fragment   fewer cells than max(200, 0.5 % of the pool)
    clean      everything else

Fragments keep their programme label — they are small, not wrong.

Subclusters are not screened for cells of a foreign lineage. The original analysis did that
against the per-slide annotation of the initial clustering, which the GEO deposit does not
carry, and scoring a handful of lineage markers is not specific enough to stand in for it:
pericytes and myofibroblasts score like endothelium often enough that real cells would be
thrown away. What the scores say is reported per subcluster (`dominant_lineage`,
`dominant_pct`) and left as a judgement for the reader.

Usage:
    python scripts/subcluster_from_geo.py
    python scripts/merge_subclusters_from_geo.py

Output:
    03.data_processed/subclustered/subcluster_merged_mapping.csv   per subcluster
    03.data_processed/subclustered/cell_states.csv                 per cell
"""

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
from pool_config import POOLS

ROOT = Path(__file__).resolve().parent.parent
SUB_DIR = ROOT / "03.data_processed/subclustered"

FRAG_ABS = 200
FRAG_REL = 0.005
SEED = 42


def lineage_of_cells(a):
    """For each cell, the pool whose lineage genes it expresses most strongly."""
    b = a.copy()
    b.X = b.layers["lognorm"].copy() if "lognorm" in b.layers else b.X
    for name, pool in POOLS.items():
        genes = [g for g in pool["filter_genes"] if g in b.var_names]
        sc.tl.score_genes(b, genes, score_name=f"pool_{name}", random_state=SEED)
    cols = [f"pool_{n}" for n in POOLS]
    return b.obs[cols].idxmax(axis=1).str.replace("pool_", "", regex=False)


def classify(obs, pool_name, pool_total):
    """qc_status and merged label for every subcluster of one pool."""
    frag_thr = max(FRAG_ABS, FRAG_REL * pool_total)
    rows = []
    for sub, cells in obs.groupby("leiden_sub", observed=True):
        n = len(cells)
        programme = cells["program"].iloc[0]
        share = cells["lineage"].value_counts(normalize=True) * 100
        dominant, dominant_pct = share.index[0], share.iloc[0]

        status = "fragment" if n < frag_thr else "clean"
        merged = programme
        rows.append(
            {
                "pool": pool_name,
                "leiden_sub": str(sub),
                "n_cells": n,
                "program": programme,
                "dominant_lineage": dominant,
                "dominant_pct": round(dominant_pct, 1),
                "qc_status": status,
                "subtype_merged": merged,
            }
        )
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub-dir", default=SUB_DIR, type=Path)
    ap.add_argument("--pools", nargs="*", default=list(POOLS))
    args = ap.parse_args()

    mappings, cells_out = [], []
    for name in args.pools:
        path = args.sub_dir / f"{name}.h5ad"
        if not path.exists():
            print(f"{name}: {path} 없음, 건너뜀")
            continue
        a = sc.read_h5ad(path)
        a.obs["lineage"] = lineage_of_cells(a).values
        frag_thr = max(FRAG_ABS, FRAG_REL * a.n_obs)
        table = classify(a.obs, name, a.n_obs)
        mappings.append(table)

        by_status = table.groupby("qc_status")["n_cells"].agg(["size", "sum"])
        print(
            f"\n=== {name} ({POOLS[name]['label']}) — {a.n_obs:,} cells, "
            f"fragment 기준 {frag_thr:.0f} cells ==="
        )
        print(by_status.to_string())
        foreign = table[table["dominant_lineage"] != name]
        if len(foreign):
            print(
                "다른 계통 점수가 높은 서브클러스터 (판정 아님, 참고용):",
                ", ".join(
                    f"sub{r.leiden_sub}->{r.dominant_lineage} ({r.n_cells:,})"
                    for r in foreign.nlargest(5, "n_cells").itertuples()
                ),
            )

        merged = dict(zip(table["leiden_sub"], table["subtype_merged"], strict=True))
        status = dict(zip(table["leiden_sub"], table["qc_status"], strict=True))
        key = a.obs["leiden_sub"].astype(str)
        cells_out.append(
            pd.DataFrame(
                {
                    "cell_id": a.obs_names,
                    "pool": name,
                    "leiden_sub": key.values,
                    "program": a.obs["program"].values,
                    "subtype_merged": key.map(merged).values,
                    "qc_status": key.map(status).values,
                }
            )
        )

    mapping = pd.concat(mappings, ignore_index=True)
    cells = pd.concat(cells_out, ignore_index=True)
    mapping.to_csv(args.sub_dir / "subcluster_merged_mapping.csv", index=False)
    cells.to_csv(args.sub_dir / "cell_states.csv", index=False)

    print(f"\n{len(mapping)} subclusters -> {mapping['subtype_merged'].nunique()} labels")
    print(mapping.groupby("qc_status")["n_cells"].agg(["size", "sum"]).to_string())
    print(f"\n{len(cells):,} cells written to {args.sub_dir / 'cell_states.csv'}")


if __name__ == "__main__":
    main()
