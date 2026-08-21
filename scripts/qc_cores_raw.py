"""Stage 1 core QC: technical quality, from the GEO files alone.

Judges each TMA core on what the raw output can tell us — how many cells were segmented,
how much signal each carries, how densely the tissue is covered — without any notion of
which cells are tumor. Cell-type composition is a property of the clustering and is
assessed afterwards by `qc_cores_composition.py`.

A core passes when all three hold:

    n_cells >= 50            enough cells to analyse
    median_transcripts >= 20 enough signal per cell
    tissue_density >= 50     cells per mm2 of occupied tissue, not of the whole core

Tissue area is the area actually covered by cells: cell centroids are binned on a 50 um
grid and the occupied bins are counted, so an empty rim or a torn core does not count as
tissue. The other columns (median genes, cell area, low-quality fraction, control-probe
fraction) are reported for inspection and do not gate.

Usage:
    python scripts/qc_cores_raw.py                    # uses 01.data_raw/GSE343808
    python scripts/qc_cores_raw.py --geo-dir <dir>

Output:
    02.tma_core_qc/core_qc_raw.csv
"""

import argparse
import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent.parent
GEO_DIR = ROOT / "01.data_raw/GSE343808"
OUT = ROOT / "02.tma_core_qc/core_qc_raw.csv"

SLIDES = ["TMA1", "TMA2"]
GRID_UM = 50.0
MIN_CELLS = 50
MIN_MEDIAN_TX = 20
MIN_DENSITY = 50
LOW_QUALITY_TX = 10
CONTROL_COLS = [
    "control_probe_counts",
    "genomic_control_counts",
    "control_codeword_counts",
    "unassigned_codeword_counts",
    "deprecated_codeword_counts",
]


def read_table(path):
    raw = gzip.open(path, "rb").read() if path.suffix == ".gz" else path.read_bytes()
    name = path.name[:-3] if path.suffix == ".gz" else path.name
    buf = io.BytesIO(raw)
    return pd.read_parquet(buf) if name.endswith(".parquet") else pd.read_csv(buf)


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


def tissue_area_um2(x, y, grid_um=GRID_UM):
    """Area actually covered by cells: occupied bins of a grid over the centroids."""
    if len(x) == 0:
        return np.nan
    xs = np.floor((x - x.min()) / grid_um).astype(int)
    ys = np.floor((y - y.min()) / grid_um).astype(int)
    return len(set(zip(xs, ys, strict=True))) * grid_um**2


def genes_per_cell(geo_dir, slide, cell_ids):
    """Detected genes per cell, from the count matrix."""
    a = sc.read_10x_h5(find(geo_dir, slide, "cell_feature_matrix.h5"))
    a.var_names_make_unique()
    n_genes = pd.Series(np.asarray((a.X > 0).sum(axis=1)).ravel(), index=a.obs_names)
    return n_genes.reindex(cell_ids)


def core_metrics(cells, boxes, slide):
    records = []
    for _, box in boxes.iterrows():
        sub = cells[cells["core_id"] == box["core_id"]]
        n = len(sub)
        tissue = tissue_area_um2(sub["x_centroid"].to_numpy(), sub["y_centroid"].to_numpy())
        control = sub[[c for c in CONTROL_COLS if c in sub]].to_numpy().sum()
        records.append(
            {
                "slide": slide,
                "core_id": box["core_id"],
                "sample_group": box.get("sample_group", ""),
                "n_cells": n,
                "core_area_um2": box["area_um2"],
                "cell_density_per_mm2": round(n / box["area_um2"] * 1e6, 2),
                "tissue_area_um2": round(tissue, 1),
                "tissue_density_per_mm2": round(n / tissue * 1e6, 2),
                "median_transcripts": sub["transcript_counts"].median(),
                "median_genes": sub["n_genes"].median(),
                "median_cell_area_um2": round(sub["cell_area"].median(), 2),
                "pct_low_quality": round(
                    (sub["transcript_counts"] < LOW_QUALITY_TX).mean() * 100, 2
                ),
                "pct_control_counts": round(control / max(sub["total_counts"].sum(), 1) * 100, 3),
            }
        )
    return pd.DataFrame(records)


def assign_cores(cells, boxes):
    cells["core_id"] = "unassigned"
    for _, b in boxes.iterrows():
        inside = cells["x_centroid"].between(b["x_min_um"], b["x_max_um"]) & cells[
            "y_centroid"
        ].between(b["y_min_um"], b["y_max_um"])
        cells.loc[inside, "core_id"] = b["core_id"]
    return cells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geo-dir", default=GEO_DIR, type=Path)
    ap.add_argument("--out", default=OUT, type=Path)
    args = ap.parse_args()

    frames = []
    for slide in SLIDES:
        cells = read_table(find(args.geo_dir, slide, "cells.parquet"))
        boxes = read_table(find(args.geo_dir, slide, "core_bounding_box.csv"))
        cells = assign_cores(cells, boxes)
        cells["n_genes"] = genes_per_cell(args.geo_dir, slide, cells["cell_id"]).to_numpy()
        frames.append(core_metrics(cells, boxes, slide))

    qc = pd.concat(frames, ignore_index=True)
    qc["pass_ncells"] = qc["n_cells"] >= MIN_CELLS
    qc["pass_transcripts"] = qc["median_transcripts"] >= MIN_MEDIAN_TX
    qc["pass_density"] = qc["tissue_density_per_mm2"] >= MIN_DENSITY
    qc["qc_pass_raw"] = qc[["pass_ncells", "pass_transcripts", "pass_density"]].all(axis=1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(args.out, index=False)
    failed = qc.loc[~qc["qc_pass_raw"], ["slide", "core_id"]].agg("_".join, axis=1).tolist()
    print(f"{len(qc)} cores | pass {qc['qc_pass_raw'].sum()} | fail {(~qc['qc_pass_raw']).sum()}")
    if failed:
        print("failed:", ", ".join(failed))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
