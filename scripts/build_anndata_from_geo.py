"""Assemble the downloaded GEO files into the object the clustering pipeline starts from.

Reads the three files per slide fetched by `load_from_geo.py` and writes one AnnData with
raw counts, spatial coordinates, slide and TMA core assignment. Cells outside every core
bounding box are kept and labelled `unassigned`; the clustering step drops them along with
the cores that fail QC.

Usage:
    python scripts/load_from_geo.py --token <reviewer token>
    python scripts/build_anndata_from_geo.py

Output:
    03.data_processed/geo_slides.h5ad
"""

import argparse
import gzip
import io
from pathlib import Path

import anndata
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent.parent
GEO_DIR = ROOT / "01.data_raw/GSE343808"
OUT = ROOT / "03.data_processed/geo_slides.h5ad"

SLIDES = ["TMA1", "TMA2"]
CELL_COLS = ["x_centroid", "y_centroid", "transcript_counts", "cell_area", "nucleus_area"]


def find(geo_dir, slide, suffix):
    hits = sorted(geo_dir.glob(f"*_{slide}_{suffix}")) + sorted(
        geo_dir.glob(f"*_{slide}_{suffix}.gz")
    )
    if not hits:
        raise FileNotFoundError(
            f"{slide}: no file matching *_{slide}_{suffix} in {geo_dir} "
            "— run scripts/load_from_geo.py first"
        )
    return hits[0]


def read_table(path):
    """Read a parquet or csv file that may or may not be gzipped."""
    raw = gzip.open(path, "rb").read() if path.suffix == ".gz" else path.read_bytes()
    name = path.name[:-3] if path.suffix == ".gz" else path.name
    buf = io.BytesIO(raw)
    return pd.read_parquet(buf) if name.endswith(".parquet") else pd.read_csv(buf)


def assign_cores(cells, boxes):
    """Label each cell with the TMA core whose bounding box contains it."""
    cells["core_id"] = "unassigned"
    cells["sample_group"] = "unassigned"
    for _, b in boxes.iterrows():
        inside = cells["x_centroid"].between(b["x_min_um"], b["x_max_um"]) & cells[
            "y_centroid"
        ].between(b["y_min_um"], b["y_max_um"])
        cells.loc[inside, "core_id"] = b["core_id"]
        cells.loc[inside, "sample_group"] = b.get("sample_group", "unassigned")
    return cells


def load_slide(geo_dir, slide):
    adata = sc.read_10x_h5(find(geo_dir, slide, "cell_feature_matrix.h5"))
    adata.var_names_make_unique()

    cells = read_table(find(geo_dir, slide, "cells.parquet")).set_index("cell_id")
    boxes = read_table(find(geo_dir, slide, "core_bounding_box.csv"))
    cells = assign_cores(cells.reset_index(), boxes).set_index("cell_id")

    adata.obs_names = [f"{slide}_{c}" for c in adata.obs_names]
    cells.index = [f"{slide}_{c}" for c in cells.index]
    shared = adata.obs_names[adata.obs_names.isin(cells.index)]
    adata = adata[shared].copy()
    cells = cells.loc[shared]

    for col in CELL_COLS + ["core_id", "sample_group"]:
        adata.obs[col] = cells[col].values
    adata.obsm["X_spatial"] = cells[["x_centroid", "y_centroid"]].values
    adata.obs["slide"] = slide
    print(
        f"[{slide}] {adata.n_obs:,} cells x {adata.n_vars:,} genes, "
        f"{(adata.obs['core_id'] != 'unassigned').sum():,} inside a core"
    )
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geo-dir", default=GEO_DIR, type=Path, help="the downloaded GEO files")
    ap.add_argument("--out", default=OUT, type=Path)
    args = ap.parse_args()

    adata = anndata.concat([load_slide(args.geo_dir, s) for s in SLIDES], join="inner")
    adata.obs_names_make_unique()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")
    print(f"wrote {args.out}  ({adata.n_obs:,} cells x {adata.n_vars:,} genes)")


if __name__ == "__main__":
    main()
