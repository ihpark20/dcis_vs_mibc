"""Tumor boundary, from how much microenvironment each epithelial cell touches.

The epithelial pool is every cell called Tumor or Myoepithelial. For each of its cells,
f_env is the fraction of its neighbours within 30 um that lie outside the pool — stroma,
immune cells, vessels. A cell deep inside a duct or a tumor mass touches only its own kind
and has f_env near zero; a cell at the edge touches the microenvironment and f_env rises.
Boundary is therefore a continuum, not a category, and the threshold is a choice: f_env has
no natural gap. TAU = 0.20 is the value used in the paper, and the sensitivity table
records what other thresholds would have given.

Two conditions keep the measure meaningful:

    n_nbr >= 3       a cell with almost no neighbours has an unstable fraction
    comp_size >= 30  the cell must sit in a connected epithelial mass of some size,
                     not in a handful of scattered cells

Cores are dropped here, not earlier: TMA1 D6 supplies most of two tumor subclusters on its
own (see the README), and boundaries need the whole epithelial compartment of a core, so
no subcluster-level exclusion is applied.

Usage:
    python scripts/major_celltypes_from_geo.py
    python scripts/boundary_from_geo.py

Output:
    03.data_processed/boundary/pool_boundary_cells.csv   per epithelial cell
    03.data_processed/boundary/tau_sensitivity.csv       boundary fraction by threshold
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from sklearn.neighbors import radius_neighbors_graph

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_qc_passed_from_geo.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
OUT_DIR = ROOT / "03.data_processed/boundary"

EPITHELIAL = {"Tumor", "Myoepithelial"}
EXCLUDE_CORES = {"TMA1_D6"}
RADIUS = 30.0
MIN_NBR = 3
MIN_MASS = 30
TAU = 0.20
TAU_GRID = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]


def core_boundary(sub):
    """f_env, neighbour count and epithelial mass for one core."""
    pos = sub[["x_centroid", "y_centroid"]].to_numpy()
    is_epi = sub["major_celltype"].isin(EPITHELIAL).to_numpy()
    if is_epi.sum() < MIN_MASS:
        return None

    epi_pos = pos[is_epi]
    epi_idx = np.where(is_epi)[0]
    tree = cKDTree(pos)
    neighbours = tree.query_ball_point(epi_pos, RADIUS)

    f_env = np.full(is_epi.sum(), np.nan)
    n_nbr = np.zeros(is_epi.sum(), dtype=int)
    for i, (self_idx, neigh) in enumerate(zip(epi_idx, neighbours, strict=True)):
        neigh = [k for k in neigh if k != self_idx]
        n_nbr[i] = len(neigh)
        if neigh:
            f_env[i] = (~is_epi[neigh]).mean()

    graph = radius_neighbors_graph(epi_pos, RADIUS, include_self=False)
    _, component = connected_components(graph, directed=False)
    comp_size = np.bincount(component)[component]

    return pd.DataFrame(
        {
            "cell_id": sub["cell_id"].to_numpy()[is_epi],
            "core": sub["core"].iloc[0],
            "sample_group": sub["sample_group"].iloc[0],
            "x_centroid": epi_pos[:, 0],
            "y_centroid": epi_pos[:, 1],
            "major_celltype": sub["major_celltype"].to_numpy()[is_epi],
            "f_env": np.round(f_env, 4),
            "n_nbr": n_nbr,
            "comp_size": comp_size,
        }
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    ap.add_argument("--tau", default=TAU, type=float)
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered, backed="r")
    cells = a.obs[["slide", "core_id", "sample_group", "x_centroid", "y_centroid"]].copy()
    cells["cell_id"] = a.obs_names
    cells["core"] = cells["slide"].astype(str) + "_" + cells["core_id"].astype(str)

    major = pd.read_csv(args.major, usecols=["cell_id", "major_celltype"])
    cells = cells.merge(major, on="cell_id", how="left")
    cells["major_celltype"] = cells["major_celltype"].fillna("unassigned")
    cells = cells[~cells["core"].isin(EXCLUDE_CORES)]
    print(
        f"{len(cells):,} cells in {cells['core'].nunique()} cores, "
        f"{cells['major_celltype'].isin(EPITHELIAL).sum():,} epithelial"
    )

    frames = [f for _, sub in cells.groupby("core") if (f := core_boundary(sub)) is not None]
    res = pd.concat(frames, ignore_index=True)

    valid = (res["n_nbr"] >= MIN_NBR) & (res["comp_size"] >= MIN_MASS) & res["f_env"].notna()
    res["valid"] = valid
    res["zone"] = np.where(
        valid, np.where(res["f_env"] >= args.tau, "boundary", "interior"), "excluded"
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out_dir / "pool_boundary_cells.csv", index=False)

    fe = res.loc[valid, "f_env"]
    print(
        f"\nepithelial pool {len(res):,} | valid (nbr>={MIN_NBR}, mass>={MIN_MASS}) {valid.sum():,}"
    )
    print(
        "f_env percentiles: "
        + ", ".join(f"p{q}={np.percentile(fe, q):.3f}" for q in (10, 25, 50, 75, 90))
    )
    print(f"fully interior (f_env=0): {(fe == 0).mean() * 100:.1f}%")

    sens = pd.DataFrame(
        [
            {
                "tau": t,
                "n_boundary": int((fe >= t).sum()),
                "boundary_pct": round((fe >= t).mean() * 100, 1),
            }
            for t in TAU_GRID
        ]
    )
    sens.to_csv(args.out_dir / "tau_sensitivity.csv", index=False)
    print("\n" + sens.to_string(index=False))
    print(
        f"\nat tau={args.tau}: {(res['zone'] == 'boundary').sum():,} boundary, "
        f"{(res['zone'] == 'interior').sum():,} interior"
    )
    print(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
