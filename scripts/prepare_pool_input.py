"""Build the object the original subclustering scripts expect.

Those scripts were written against the analysis object of the study and read a few columns
this pipeline names differently or has not carried this far: the initial cluster as a plain
Leiden number, the DCIS/mDCIS group under the name `pathology`, and the per-slide cell type
each cell had before the integrated clustering. This assembles them so the scripts can run
unchanged, and drops the cores the analysis excludes.

Usage:
    python scripts/cluster_from_geo.py --labels 03.data_processed/integrated_cluster_labels.csv.gz
    python scripts/prepare_pool_input.py

Output:
    03.data_processed/pool_input.h5ad
"""

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
from analysis_exclusions import drop_excluded

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_qc_passed_from_geo.h5ad"
QC_DIR = ROOT / "02.tma_core_qc"
OUT = ROOT / "03.data_processed/pool_input.h5ad"

SLIDES = ["TMA1", "TMA2"]
CANCER_TYPES = {"cancer"}
IMMUNE_TYPES = {
    "B",
    "CD8T",
    "NKT",
    "T",
    "Treg",
    "perivascular NK",
    "plasma cell",
    "mast cell",
    "M1",
    "M2",
    "moDC",
    "cDC",
    "pDC",
}
STROMA_TYPES = {"CAF", "endothelial cell", "myoepithelial cell"}


def celltype_group(celltype):
    if celltype in CANCER_TYPES:
        return "cancer"
    if celltype in IMMUNE_TYPES:
        return "immune"
    if celltype in STROMA_TYPES:
        return "stroma"
    return "other"


def prior_celltypes(qc_dir):
    """The cell type each cell carried in the per-slide annotation."""
    frames = []
    for slide in SLIDES:
        clusters = pd.read_csv(qc_dir / f"{slide}.clusters.csv")
        clusters.columns = ["cell_id", "cluster_id"]
        annotation = pd.read_csv(qc_dir / f"{slide}.cluster_celltype_annotation.csv")
        celltype = annotation.set_index("cluster_id")[annotation.columns[-2]].to_dict()
        clusters["cell_id"] = slide + "_" + clusters["cell_id"]
        clusters["celltype_final"] = clusters["cluster_id"].map(celltype).fillna("unknown")
        frames.append(clusters[["cell_id", "celltype_final"]])
    return pd.concat(frames, ignore_index=True).set_index("cell_id")["celltype_final"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--qc-dir", default=QC_DIR, type=Path)
    ap.add_argument("--out", default=OUT, type=Path)
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered)
    core = a.obs["slide"].astype(str) + "_" + a.obs["core_id"].astype(str)
    a = a[drop_excluded(core.to_frame("core"), "core").index].copy()

    a.obs["leiden"] = a.obs["cluster"].astype(str).str.removeprefix("CL").astype("category")
    a.obs["pathology"] = a.obs["sample_group"].astype(str)
    prior = prior_celltypes(args.qc_dir)
    a.obs["celltype_final"] = prior.reindex(a.obs_names).fillna("unknown").to_numpy()
    a.obs["celltype_group"] = [celltype_group(c) for c in a.obs["celltype_final"]]
    if "counts" not in a.layers:
        a.layers["counts"] = a.X.copy()

    print(f"{a.n_obs:,} cells x {a.n_vars} genes in {a.obs['leiden'].nunique()} clusters")
    print(a.obs["celltype_group"].value_counts().to_string())
    a.write_h5ad(args.out, compression="gzip")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
