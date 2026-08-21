"""Download the GSE343808 files that the analysis needs, from GEO.

The cell-level expression matrix is not part of this repository; it is deposited at GEO.
Three files per slide are enough for everything downstream — about 60 MB per slide, not
the 18 GB morphology image or the 1 GB transcript table.

While the series is private, pass the reviewer access token with --token: GEO then opens a
session in which the supplementary files download normally. The token is used for the
requests only — it is never written to disk.

Usage:
    python scripts/load_from_geo.py --token <reviewer token>   # while private
    python scripts/load_from_geo.py                            # once GSE343808 is public

Downloads next to the scripts, into 01.data_raw/GSE343808/ (override with --geo-dir):
    GSM9963125_TMA1_cell_feature_matrix.h5      cell x gene counts
    GSM9963125_TMA1_cells.parquet.gz            centroids and per-cell QC metrics
    GSM9963125_TMA1_core_bounding_box.csv.gz    TMA core extents and sample group
    GSM9963126_TMA2_...                         the same three files for TMA2

Then build the AnnData object with `build_anndata_from_geo.py`.
"""

import argparse
import urllib.parse
import urllib.request
from http.cookiejar import CookieJar
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GEO_DIR = ROOT / "01.data_raw/GSE343808"

SERIES = "GSE343808"
SAMPLES = {"TMA1": "GSM9963125", "TMA2": "GSM9963126"}
NEEDED = ["cell_feature_matrix.h5", "cells.parquet.gz", "core_bounding_box.csv.gz"]
ACC_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
DL_URL = "https://www.ncbi.nlm.nih.gov/geo/download/"


def geo_session(token=None):
    """Open a GEO session; with a token it also unlocks a private series."""
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(CookieJar()))
    query = {"acc": SERIES}
    if token:
        query["token"] = token
    opener.open(f"{ACC_URL}?{urllib.parse.urlencode(query)}", timeout=120).read()
    return opener


def download(geo_dir, token=None, force=False):
    """Fetch the files the analysis needs, skipping ones already present."""
    opener = geo_session(token)
    geo_dir.mkdir(parents=True, exist_ok=True)
    print(f"downloading to {geo_dir}")
    for slide, acc in SAMPLES.items():
        for base in NEEDED:
            name = f"{acc}_{slide}_{base}"
            dst = geo_dir / name
            if dst.exists() and dst.stat().st_size > 0 and not force:
                print(f"  {name:46s} {dst.stat().st_size / 1e6:8.1f} MB  (already present)")
                continue
            query = urllib.parse.urlencode({"acc": acc, "format": "file", "file": name})
            with opener.open(f"{DL_URL}?{query}", timeout=1800) as response:
                body = response.read()
            if body[:1] == b"<":  # GEO answers with an HTML error page
                raise RuntimeError(f"{name}: GEO refused the download — is the token right?")
            dst.write_bytes(body)
            print(f"  {name:46s} {len(body) / 1e6:8.1f} MB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token", help="GEO reviewer access token, while the series is private")
    ap.add_argument("--geo-dir", default=GEO_DIR, type=Path, help="where the GEO files go")
    ap.add_argument("--force", action="store_true", help="download again even if present")
    args = ap.parse_args()

    download(args.geo_dir, token=args.token, force=args.force)
    print(f"files are in {args.geo_dir.resolve()}")


if __name__ == "__main__":
    main()
