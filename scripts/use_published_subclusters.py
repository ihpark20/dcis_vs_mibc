"""Take the cell states the paper assigned, instead of subclustering again.

Subclustering the seven pools takes the better part of an hour and lands close to, but not
exactly on, the paper's partition. The assignment it produced travels with this repository,
so an analysis that wants the published cell states can start from them directly. The
major cell types come with them, so this also writes the table the boundary, polarization
and communication steps read; everything after this point runs the same either way.

Usage:
    python scripts/use_published_subclusters.py

Output:
    03.data_processed/subclustered/cell_states.csv
    03.data_processed/subclustered/major_celltypes.csv
"""

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LABELS = ROOT / "03.data_processed/subcluster_labels.csv.gz"
OUT_DIR = ROOT / "03.data_processed/subclustered"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default=LABELS, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    args = ap.parse_args()

    cells = pd.read_csv(args.labels)
    cells["program"] = cells["subtype_merged"]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cells.to_csv(args.out_dir / "cell_states.csv", index=False)

    print(f"{len(cells):,} cells across {cells['pool'].nunique()} pools")
    print(
        cells.groupby("pool", observed=True)["subtype_merged"]
        .nunique()
        .rename("cell states")
        .to_string()
    )
    print(f"\nqc_status: {cells['qc_status'].value_counts().to_dict()}")
    major = cells[["cell_id", "pool", "leiden_sub", "subtype_merged"]].copy()
    major["major_celltype"] = cells["celltype_level1"]
    major.to_csv(args.out_dir / "major_celltypes.csv", index=False)
    print("\n" + major["major_celltype"].value_counts().to_string())
    print(f"\nwrote {args.out_dir / 'cell_states.csv'} and major_celltypes.csv")


if __name__ == "__main__":
    main()
