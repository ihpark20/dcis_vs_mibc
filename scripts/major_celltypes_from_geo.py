"""Fold the pool cell states back into major cell types.

The pools give fine-grained states; most analyses need the coarse lineages behind them.
This maps every (pool, programme) pair to one major cell type. Two pools split: the tumor
pool holds both tumor and myoepithelial cells, and pericytes appear in both the
endothelial and the CAF pool, so the programme decides in those cases.

The panel carries almost no neutrophil markers — of the usual set only S100A8, CEACAM8 and
ITGAM are on it, none of them specific — so no pool resolves neutrophils and they are not
among the major types here.

Usage:
    python scripts/merge_subclusters_from_geo.py
    python scripts/major_celltypes_from_geo.py

Output:
    03.data_processed/subclustered/major_celltypes.csv   per cell
    03.data_processed/subclustered/major_celltype_counts.csv
"""

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SUB_DIR = ROOT / "03.data_processed/subclustered"

# the major type a pool gives, unless the programme says otherwise
POOL_MAJOR = {
    "tnk": "T_NK",
    "endothelial": "Endothelial",
    "myeloid": "Macrophage_Mono",
    "tumor": "Tumor",
    "bplasma": "B_cell",
    "caf": "CAF_Fibroblast",
    "mast": "Mast",
}

PROGRAM_MAJOR = {
    ("tumor", "myoepithelial"): "Myoepithelial",
    ("endothelial", "pericyte_mural"): "Pericyte",
    ("caf", "pericyte_mural"): "Pericyte",
    ("caf", "adipocyte"): "Adipocyte",
    ("caf", "preadipo"): "Adipocyte",
    ("bplasma", "plasma"): "Plasma",
    ("bplasma", "igG"): "Plasma",
    ("bplasma", "igM"): "Plasma",
    ("myeloid", "cDC1"): "Dendritic",
    ("myeloid", "cDC2"): "Dendritic",
    ("myeloid", "pDC"): "Dendritic",
    ("myeloid", "mregDC"): "Dendritic",
    ("myeloid", "langerhans"): "Dendritic",
}


def major_of(pool, programme):
    return PROGRAM_MAJOR.get((pool, programme), POOL_MAJOR[pool])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub-dir", default=SUB_DIR, type=Path)
    args = ap.parse_args()

    cells = pd.read_csv(args.sub_dir / "cell_states.csv")
    cells["major_celltype"] = [
        major_of(p, s) for p, s in zip(cells["pool"], cells["subtype_merged"], strict=True)
    ]
    cells.to_csv(args.sub_dir / "major_celltypes.csv", index=False)

    counts = (
        cells.groupby(["major_celltype", "pool"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values("n_cells", ascending=False)
    )
    counts.to_csv(args.sub_dir / "major_celltype_counts.csv", index=False)

    print(cells["major_celltype"].value_counts().to_string())
    print(f"\n{cells['major_celltype'].nunique()} major cell types over {len(cells):,} cells")
    print(f"wrote {args.sub_dir / 'major_celltypes.csv'}")


if __name__ == "__main__":
    main()
