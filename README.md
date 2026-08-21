# Spatial transcriptomics of DCIS and microinvasive breast carcinoma

Code and source data behind the figures and supplementary tables of the manuscript.
Two breast tissue microarrays (41 cores) were profiled with 10x Genomics Xenium using a
374-gene custom panel (280-gene Human Breast panel + 94-gene immuno-oncology add-on) and
compared between pure DCIS and microinvasive carcinoma (mDCIS).

# Environment
All computation ran in a single conda environment (py312_r44): Python 3.12.8, scanpy 1.11.1, anndata 0.11.4, harmonypy 0.2.0, PyDESeq2 0.5.4, scikit-learn 1.5.2, SciPy 1.15.1, NumPy 2.2.6, pandas 2.2.3, statsmodels 0.14.4, PyMuPDF (figure assembly). Figures were rendered in R 4.4.3 with ggplot2, patchwork, ggalluvial and ComplexHeatmap.


## Getting the data (GEO GSE343808)

The cell-level expression matrix is not in this repository; it is deposited at GEO under
accession **GSE343808**, one sample per tissue microarray
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343808).

**The series is still private while the manuscript is under review.** Reaching the files
needs a reviewer access token: editors and reviewers receive one with the manuscript, and
anyone else can request it from the corresponding author. Append the token to the accession
URL as GEO instructs, and the supplementary files become downloadable. This note will be
removed once the series is released.

Download these files per slide (`GSM9963125_TMA1_*` and the corresponding `GSM*_TMA2_*`)
into a single directory:

| File | What it is | Needed for |
|---|---|---|
| `*_TMA1_cell_feature_matrix.h5` | cell x gene counts (374 genes) | everything |
| `*_TMA1_cells.parquet.gz` | centroids, cell and nucleus area, transcript counts | everything |
| `*_TMA1_core_bounding_box.csv.gz` | TMA core extents and DCIS/mDCIS group | assigning cells to cores |
| `*_TMA1_transcripts.parquet.gz` | individual transcript coordinates | not used by this pipeline |
| `*_TMA1_cell_boundaries.parquet.gz`, `*_nucleus_boundaries.parquet.gz` | segmentation polygons | not used by this pipeline |
| `*_TMA1_morphology.ome.tif.gz` | morphology image (18 GB) | not used by this pipeline |

The analysis needs only the first three files of each slide — about 60 MB per slide, not
the 18 GB image or the 1 GB transcript table.

