# Spatial transcriptomics of DCIS and microinvasive breast carcinoma

Code and data to reproduce the findings in '[Spatial transcriptomics reveals tumor-stroma interface remodeling in HER2-positive ductal carcinoma in situ with microinvasion]' (under review). 
Two breast tissue microarrays (41 cores) were profiled with 10x Genomics Xenium using a
374-gene custom panel (280-gene Human Breast panel + 94-gene immuno-oncology add-on) and
compared between pure DCIS and microinvasive carcinoma (mDCIS).

## Environment

All computation ran in a single conda environment (`py312_r44`): Python 3.12.8,
scanpy 1.11.1, anndata 0.11.4, harmonypy 0.2.0, PyDESeq2 0.5.4, scikit-learn 1.5.2,
SciPy 1.15.1, NumPy 2.2.6, pandas 2.2.3, statsmodels 0.14.4.
Figures were rendered in R 4.4.3 with ggplot2, patchwork, ggalluvial and ComplexHeatmap.

## Getting the data (GEO GSE343808)

The cell-level expression matrix is not in this repository; it is deposited at GEO under
accession **GSE343808**, one sample per tissue microarray
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343808).

**The series is still private while the manuscript is under review.** Reaching the files
needs a reviewer access token: editors and reviewers receive one with the manuscript, and
anyone else can request it from the corresponding author. Opening the accession URL with
the token starts a session in which the per-sample supplementary file links download
normally. This note will be removed once the series is released.

Download these files per slide (`GSM9963125_TMA1_*` and `GSM9963126_TMA2_*`) into a
single directory:

| File | What it is | Needed for |
|---|---|---|
| `*_TMA{1,2}_cell_feature_matrix.h5` | cell x gene counts (374 genes) | everything |
| `*_TMA{1,2}_cells.parquet.gz` | centroids, cell and nucleus area, transcript counts | everything |
| `*_TMA{1,2}_core_bounding_box.csv.gz` | TMA core extents and DCIS/mDCIS group | assigning cells to cores |
| `*_TMA{1,2}_transcripts.parquet.gz` | individual transcript coordinates | not used by this pipeline |
| `*_TMA{1,2}_cell_boundaries.parquet.gz`, `*_nucleus_boundaries.parquet.gz` | segmentation polygons | not used by this pipeline |
| `*_TMA{1,2}_morphology.ome.tif.gz` | morphology image (18 GB) | not used by this pipeline |

The analysis needs only the first three files of each slide — about 60 MB per slide, not
the 18 GB image or the 1 GB transcript table.

Two scripts cover this:

```bash
# 1. download the six files into 01.data_raw/GSE343808/
python scripts/load_from_geo.py --token <reviewer access token>

# 2. assemble them: raw counts, spatial coordinates, slide, TMA core
python scripts/build_anndata_from_geo.py
#    -> 03.data_processed/geo_slides.h5ad
```

`load_from_geo.py` opens a GEO session with the token and fetches only the three files per
slide listed above, skipping any already present (`--force` re-fetches, `--geo-dir` puts
them elsewhere). The token is used for the requests only and is never written to disk. Once
the series is public it can be omitted. Both scripts resolve their paths from their own
location, so the data lands inside the repository wherever it is cloned.

This was checked against the deposit: the six files assemble to 546,021 cells x 374 genes,
of which 333,964 (TMA1) and 210,460 (TMA2) fall inside a TMA core. Restricted to the
QC-passing cores and the cell filters of the paper (`min_counts=10`, `min_genes=5`), this
is exactly the 520,506 cells of the published object — same cell identifiers, same raw
counts.

## QC

Two of the 41 cores are left out of every analysis: **TMA1 D8**, whose cells carry a median
of 0 transcripts, and **TMA1 D4**, in which only 0.9 % of the cells are tumor — one a
technical failure, the other a core with nothing to contribute to a DCIS versus
microinvasion comparison.

Only the outcome is provided here, as `02.tma_core_qc/tma_core_qc.csv`: one row
per core with the metrics the decision rests on, an `analysis_include` flag and a one-line
`exclude_reason`. Downstream steps read that flag rather than recomputing anything, so the
cell set matches the published one.

The procedure that produced the table is not included. Part of it cannot be rerun from the
deposit: the cell-type composition of each core came from the per-slide graph-based
clustering that Xenium writes beside the raw bundle, and that clustering is not among the
deposited files.

## Integrated Clustering Across TMA Batches (CL0-CL14)

The two microarrays were imaged as separate slides, so slide is the batch to correct for.
`cluster_from_geo.py` integrates them with Harmony and partitions the result with Leiden,
giving the fifteen clusters (CL0-CL14) the annotation is built on.

```bash
python scripts/cluster_from_geo.py
#    -> 03.data_processed/integrated_qc_passed_from_geo.h5ad
```

It keeps the cores flagged `analysis_include` in the QC table, then follows the published
settings:

| step | setting |
|---|---|
| cell filters | `min_counts=10`, `min_genes=5` |
| normalization | counts per 10,000, `log1p` |
| feature selection | 300 highly variable genes, `flavor="seurat"`, `batch_key="slide"` |
| scaling, PCA | `max_value=10`, 30 components |
| batch correction | Harmony on `slide`, `max_iter_harmony=30` |
| graph, embedding | 15 neighbors on the Harmony components, UMAP |
| clustering | Leiden, `resolution=0.5` |
| seed | 42 throughout |

The written object keeps raw counts in `X`, the Harmony embedding in
`obsm["X_pca_harmony"]`, the UMAP coordinates and the cluster labels in `obs["leiden"]`.
Pass `--all-cores` to skip the QC filter and cluster every core instead.

Expect around an hour on half a million cells; the Leiden step alone takes a good part of
it and runs on a single core. Harmony and Leiden both accumulate floating point in an
order that depends on the machine, so a rerun reproduces the cluster structure rather than
every individual cell label.

## Subclustering of Pooled Clusters (7 Pools)

## Data not included here

* The cell-level count matrix and the raw Xenium output are at GEO (GSE343808), as above.
* The immunohistochemistry images used for Figure 3B (p63, SMMHC) are not redistributed here.

