"""Attach the paper's cluster assignment to the cells built from GEO.

The initial clustering of the study is provided rather than recomputed. Reclustering half a
million cells takes about an hour and lands close to but not on the published partition —
the same cells, fifteen clusters matching one to one, most cells in the same one — and every
step downstream would then rest on an approximation of the paper rather than on the paper.
The assignment itself is small enough to travel with the code, so it does.

This drops the cores the QC excludes, attaches the published cluster of every remaining
cell, and writes the object the rest of the pipeline reads. Nothing is computed.

Usage:
    python scripts/build_anndata_from_geo.py
    python scripts/qc_cores_composition.py
    python scripts/attach_published_clusters.py

Output:
    03.data_processed/integrated_clusters.h5ad
"""

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "03.data_processed/geo_slides.h5ad"
LABELS = ROOT / "03.data_processed/integrated_cluster_labels.csv.gz"
QC = ROOT / "02.tma_core_qc/tma_core_qc.csv"
OUT = ROOT / "03.data_processed/integrated_clusters.h5ad"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=IN, type=Path)
    ap.add_argument("--labels", default=LABELS, type=Path)
    ap.add_argument("--qc", default=QC, type=Path)
    ap.add_argument("--out", default=OUT, type=Path)
    args = ap.parse_args()

    a = sc.read_h5ad(args.input)
    qc = pd.read_csv(args.qc)
    keep = qc[qc["analysis_include"]]
    passing = set(keep["slide"] + "_" + keep["core_id"])
    dropped = (
        qc.loc[~qc["analysis_include"], "slide"] + "_" + qc.loc[~qc["analysis_include"], "core_id"]
    ).tolist()
    core = a.obs["slide"].astype(str) + "_" + a.obs["core_id"].astype(str)
    a = a[core.isin(passing)].copy()
    print(f"cores kept: {len(passing)} (dropped {', '.join(dropped)})  ->  {a.n_obs:,} cells")

    labels = pd.read_csv(args.labels).set_index("cell_id")["cluster"]
    known = a.obs_names.intersection(labels.index)
    print(f"{len(known):,} of {a.n_obs:,} cells carry a published cluster")
    a = a[known].copy()
    a.obs["cluster"] = labels.reindex(a.obs_names).astype("category")

    if "counts" not in a.layers:
        a.layers["counts"] = a.X.copy()
    print(a.obs["cluster"].value_counts().sort_index().to_string())

    a.write_h5ad(args.out, compression="gzip")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
