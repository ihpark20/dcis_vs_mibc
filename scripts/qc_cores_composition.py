"""Stage 2 core QC: what each core is made of.

Stage 1 (`qc_cores_raw.py`) asks whether a core was measured well. This asks what is in it:
a core can be technically flawless and still carry almost no tumor, which makes it useless
for a DCIS versus microinvasion comparison but is not a quality failure. Keeping the two
apart means the technical filter never depends on knowing which cells are tumor.

Composition is read from the graph-based clustering Xenium produces for each slide, mapped
to cell types through the annotation of those clusters. Both travel with this repository —
`02.tma_core_qc/TMA{1,2}.clusters.csv` and `TMA{1,2}.cluster_celltype_annotation.csv` —
because the GEO deposit carries neither.

A core is flagged when

    pct_cancer < 3       too little tumor to contribute to the comparison
    n_celltypes < 3      too little cellular diversity to be intact tissue

Cores failing either stage are excluded from the analysis, with the reason recorded.

Usage:
    python scripts/load_from_geo.py --token <reviewer token>
    python scripts/qc_cores_raw.py
    python scripts/qc_cores_composition.py

Output:
    02.tma_core_qc/tma_core_qc.csv
"""

import argparse
import gzip
import io
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
GEO_DIR = ROOT / "01.data_raw/GSE343808"
QC_DIR = ROOT / "02.tma_core_qc"
RAW_QC = QC_DIR / "core_qc_raw.csv"
OUT = QC_DIR / "tma_core_qc.csv"

SLIDES = ["TMA1", "TMA2"]
MIN_CANCER_PCT = 3.0
MIN_CELLTYPES = 3

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


def assign_cores(cells, boxes):
    cells["core_id"] = "unassigned"
    for _, b in boxes.iterrows():
        inside = cells["x_centroid"].between(b["x_min_um"], b["x_max_um"]) & cells[
            "y_centroid"
        ].between(b["y_min_um"], b["y_max_um"])
        cells.loc[inside, "core_id"] = b["core_id"]
    return cells


def composition(geo_dir, qc_dir, slide):
    """Per-core cell-type composition, from Xenium's own clustering."""
    cells = read_table(find(geo_dir, slide, "cells.parquet"))
    boxes = read_table(find(geo_dir, slide, "core_bounding_box.csv"))
    cells = assign_cores(cells, boxes)

    clusters = pd.read_csv(qc_dir / f"{slide}.clusters.csv")
    clusters.columns = ["cell_id", "cluster_id"]
    annotation = pd.read_csv(qc_dir / f"{slide}.cluster_celltype_annotation.csv")
    celltype_col = annotation.columns[-2]
    celltype = annotation.set_index("cluster_id")[celltype_col].to_dict()

    cells = cells.merge(clusters, on="cell_id", how="left")
    cells["celltype"] = cells["cluster_id"].map(celltype).fillna("unknown")
    cells["group"] = cells["celltype"].map(celltype_group)

    rows = []
    for core_id, sub in cells[cells["core_id"] != "unassigned"].groupby("core_id"):
        share = sub["group"].value_counts(normalize=True) * 100
        rows.append(
            {
                "slide": slide,
                "core_id": core_id,
                "n_celltypes": sub["celltype"].nunique(),
                "pct_cancer": round(share.get("cancer", 0.0), 1),
                "pct_immune": round(share.get("immune", 0.0), 1),
                "pct_stroma": round(share.get("stroma", 0.0), 1),
                "pct_other": round(share.get("other", 0.0), 1),
            }
        )
    return pd.DataFrame(rows)


def reason(row):
    """One line saying why a core is out, empty for the cores that stay in."""
    if not row["qc_pass_raw"]:
        failed = []
        if not row["pass_ncells"]:
            failed.append(f"{int(row['n_cells'])} cells")
        if not row["pass_transcripts"]:
            failed.append(f"median {row['median_transcripts']:.0f} transcripts/cell")
        if not row["pass_density"]:
            failed.append(f"{row['tissue_density_per_mm2']:.0f} cells/mm2 of tissue")
        return "technical: " + ", ".join(failed)
    if row["content_flag"]:
        failed = []
        if not row["pass_cancer_pct"]:
            failed.append(f"{row['pct_cancer']:.1f}% tumor cells")
        if not row["pass_celltypes"]:
            failed.append(f"{int(row['n_celltypes'])} cell types")
        return "content: " + ", ".join(failed)
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geo-dir", default=GEO_DIR, type=Path)
    ap.add_argument("--qc-dir", default=QC_DIR, type=Path)
    ap.add_argument("--raw-qc", default=RAW_QC, type=Path)
    ap.add_argument("--out", default=OUT, type=Path)
    args = ap.parse_args()

    comp = pd.concat([composition(args.geo_dir, args.qc_dir, s) for s in SLIDES], ignore_index=True)
    qc = pd.read_csv(args.raw_qc).merge(comp, on=["slide", "core_id"], how="left", validate="1:1")

    qc["pass_cancer_pct"] = qc["pct_cancer"] >= MIN_CANCER_PCT
    qc["pass_celltypes"] = qc["n_celltypes"] >= MIN_CELLTYPES
    qc["content_flag"] = ~(qc["pass_cancer_pct"] & qc["pass_celltypes"])
    qc["analysis_include"] = qc["qc_pass_raw"] & ~qc["content_flag"]
    qc["exclude_reason"] = [reason(row) for _, row in qc.iterrows()]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(args.out, index=False)

    print(f"{len(qc)} cores | included {qc['analysis_include'].sum()}")
    for _, row in qc[~qc["analysis_include"]].iterrows():
        print(f"  {row['slide']}_{row['core_id']}: {row['exclude_reason']}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
