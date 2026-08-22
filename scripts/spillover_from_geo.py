"""How far each cell type's transcripts leak into its neighbours.

Xenium assigns a transcript to whichever segmented cell it falls in, and segmentation is not
perfect: a transcript from a tumor cell can be counted in the fibroblast beside it. The
effect is systematic, not random — it follows whatever the neighbour is — so a gene that a
cell type expresses strongly turns up in the cells around it and can be mistaken for a real
signal there. This measures that leakage per source cell type and gene, from two quantities
that have to agree for spillover to be the explanation:

    source_spec(S, g) = log2( mean of g in S / mean of g outside S )
        how much the gene belongs to the source in the first place

    prox_fc(S, g)     = log2( mean of g in non-S cells surrounded by S
                              / mean of g in non-S cells away from S )
        whether cells of other types carry more of it when S is nearby

    spillover(S, g)   = sqrt(source_spec x prox_fc)  when both are positive, else 0

The geometric mean asks for both at once: a gene the source barely expresses cannot leak,
and a gene that does not rise near the source is not leaking. A cell counts as near S when
at least half of its neighbours within 30 um are S, and as far when fewer than a tenth are.

Neutrophils are excluded, as they are everywhere else here, so the index rests on the same
cell population as the rest of the analysis.

Usage:
    python scripts/use_published_subclusters.py
    python scripts/spillover_from_geo.py

Output:
    03.data_processed/spillover/spillover_by_source.csv   source x gene
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from analysis_exclusions import EXCLUDED_CORES
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_clusters.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
OUT_DIR = ROOT / "03.data_processed/spillover"

RADIUS = 30.0
NEAR, FAR = 0.5, 0.1
MIN_NEAR, MIN_FAR = 80, 200
MIN_SOURCE = 3000
EPS = 1e-3
DROP_TYPES = {"Neutrophil"}


def neighbour_fractions(positions, codes, cores, n_types):
    """For every cell, the share of its 30 um neighbours of each cell type."""
    fractions = np.zeros((len(positions), n_types))
    for _, idx in pd.Series(cores).groupby(cores).indices.items():
        idx = np.asarray(idx)
        tree = cKDTree(positions[idx])
        pairs = tree.query_pairs(RADIUS, output_type="ndarray")
        onehot = sp.csr_matrix(
            (np.ones(len(idx)), (np.arange(len(idx)), codes[idx])), shape=(len(idx), n_types)
        )
        rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
        cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
        adjacency = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(idx), len(idx)))
        counts = adjacency @ onehot
        total = np.asarray(counts.sum(1)).ravel()
        total[total == 0] = 1
        fractions[idx] = np.asarray(counts.todense()) / total[:, None]
    return fractions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered)
    obs = a.obs.copy()
    obs["cell_id"] = a.obs_names
    obs["core"] = obs["slide"].astype(str) + "_" + obs["core_id"].astype(str)
    major = pd.read_csv(args.major, usecols=["cell_id", "major_celltype"])
    obs = obs.merge(major, on="cell_id", how="left")

    counts = a.layers["counts"] if "counts" in a.layers else a.X
    counts = sp.csr_matrix(counts) if not sp.issparse(counts) else counts.tocsr()

    keep = (
        ~obs["core"].isin(EXCLUDED_CORES)
        & obs["major_celltype"].notna()
        & ~obs["major_celltype"].isin(DROP_TYPES)
    ).to_numpy()
    counts, obs = counts[keep], obs[keep].reset_index(drop=True)
    print(f"{len(obs):,} cells, {obs['major_celltype'].nunique()} cell types")

    types = sorted(obs["major_celltype"].unique())
    index_of = {t: i for i, t in enumerate(types)}
    codes = obs["major_celltype"].map(index_of).to_numpy()
    fractions = neighbour_fractions(
        obs[["x_centroid", "y_centroid"]].to_numpy(), codes, obs["core"].to_numpy(), len(types)
    )

    def mean_of(mask):
        return np.asarray(counts[mask].mean(axis=0)).ravel() if mask.sum() else np.zeros(a.n_vars)

    rows = []
    for source in types:
        i = index_of[source]
        is_source = codes == i
        if is_source.sum() < MIN_SOURCE:
            continue
        other = ~is_source
        near = other & (fractions[:, i] >= NEAR)
        far = other & (fractions[:, i] < FAR)
        if near.sum() < MIN_NEAR or far.sum() < MIN_FAR:
            print(f"  {source}: near {near.sum()}, far {far.sum()} — too few, skipped")
            continue

        source_spec = np.log2((mean_of(is_source) + EPS) / (mean_of(other) + EPS))
        prox_fc = np.log2((mean_of(near) + EPS) / (mean_of(far) + EPS))
        spillover = np.where((source_spec > 0) & (prox_fc > 0), np.sqrt(source_spec * prox_fc), 0.0)
        for gene, spec, prox, index in zip(
            a.var_names, source_spec, prox_fc, spillover, strict=True
        ):
            rows.append(
                {
                    "source": source,
                    "gene": gene,
                    "source_spec": round(float(spec), 3),
                    "prox_fc": round(float(prox), 3),
                    "spillover_index": round(float(index), 3),
                    "n_near": int(near.sum()),
                    "n_far": int(far.sum()),
                }
            )

    out = pd.DataFrame(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_dir / "spillover_by_source.csv", index=False)

    print(f"\n{out['source'].nunique()} sources x {a.n_vars} genes")
    print("\nstrongest leak per source")
    top = out.sort_values("spillover_index", ascending=False).groupby("source").head(3)
    print(
        top.sort_values(["source", "spillover_index"], ascending=[True, False]).to_string(
            index=False
        )
    )
    print(f"\nwrote {args.out_dir / 'spillover_by_source.csv'}")


if __name__ == "__main__":
    main()
