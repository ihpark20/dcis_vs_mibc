# Spatial transcriptomics of DCIS and microinvasive breast carcinoma

Code and source data behind the figures and supplementary tables of the manuscript.
Two breast tissue microarrays (41 cores) were profiled with 10x Genomics Xenium using a
374-gene custom panel (280-gene Human Breast panel + 94-gene immuno-oncology add-on) and
compared between pure DCIS and microinvasive carcinoma (mDCIS).

Everything needed to regenerate a published panel or supplementary table is here: the
preprocessing and annotation pipeline, the per-figure preparation scripts, the rendering
scripts, and the intermediate tables the figures are drawn from. Exploratory analyses that
did not reach the manuscript are not included.

## Environment

All computation ran in a single conda environment (`py312_r44`): Python 3.12.8,
scanpy 1.11.1, anndata 0.11.4, harmonypy 0.2.0, PyDESeq2 0.5.4, scikit-learn 1.5.2,
SciPy 1.15.1, NumPy 2.2.6, pandas 2.2.3, statsmodels 0.14.4, PyMuPDF (figure assembly).
Figures were rendered in R 4.4.3 with ggplot2, patchwork, ggalluvial and ComplexHeatmap.
Python code is linted and formatted with Ruff (`ruff.toml` included).

Random seeds are fixed at every stochastic step (Harmony and UMAP `random_state=42`;
Leiden `random_state=42`; k-means `random_state=0`; NumPy `default_rng(0)`; PERMANOVA
permutation generator `default_rng(0)`).

## Layout

```
scripts/                  Python: preprocessing, annotation, per-figure preparation, figure assembly
21.figures_hires_scripts/ R: one script per published panel, plus theme_hires.R
scripts_r/                R: shared publication theme
00.panel_info/            374-gene panel definition
01.data_raw/              per-slide cluster annotation used to seed cell typing
02.slide_tma_core/        TMA core coordinates and the DCIS/mDCIS sample map
04.tables/                master per-cell annotation, subcluster summaries, programme scores
11.tables_publication/    the tables each figure and supplementary table is drawn from
docs/                     figure legends and Supplementary Texts S1-S3 (LaTeX source)
MANIFEST.json             every file in this repository, by category
```

Scripts address files by absolute path through a `PROJ` constant
(`/BiO2/codes/breast_xenium_v2_pjt1`). To run them elsewhere, point `PROJ` at the
repository root — the tree above mirrors the original layout, so nothing else changes.

Tables larger than 5 MB are gzipped (`*.csv.gz`); `pandas.read_csv` and `readr::read_csv`
open them without any change to the call.

## Pipeline order

1. `scripts/load_from_geo.py` — build the analysis object from the GEO download
   (see "Getting the data" below). `scripts/qc_tma_cores.py` is the per-core QC that ran on
   the original bundles; two cores (TMA1 D4, D8) fail and are dropped.
2. `scripts/cluster_from_geo.py` — normalization, HVG selection, PCA, Harmony integration
   over slide, Leiden clustering into 15 clusters. `integrate_slides_harmony.py` and
   `integrate_qc_passed.py` are the equivalent scripts as originally run.
3. `scripts/subcluster_*.py` — seven lineage pools re-clustered on curated lineage markers
   (see `docs/SupplementaryTextS3.tex` for the rationale and the per-pool parameters).
4. `scripts/merge_subclusters.py` -> `scripts/build_master_annotation.py` ->
   `scripts/add_celltype_level1.py` — subcluster QC flags, the per-cell master annotation,
   and the 13 major cell types.
5. `scripts/prep_*.py` — one script per analysis, each writing tables into
   `11.tables_publication/`.
6. `21.figures_hires_scripts/*.R` — one script per panel.
7. `scripts/compose_figure*.py` — assemble the panels into the manuscript figures
   (PyMuPDF; panels stay vector).

`scripts/core_exclusion.py` is the single source of truth for which cores and subclusters
are excluded from analysis, and is imported wherever that matters.

## Figure index

Each rendering script is named for the panel it draws. The panel PDF it writes keeps the
internal working name (`Fig01m_cohort_maps` and so on), which is how `compose_figure*.py`
picks the panels up.

| Figure | Panel | Rendering script (`21.figures_hires_scripts/`) | Source table (`11.tables_publication/`, `04.tables/`) | Preparation script (`scripts/`) |
|---|---|---|---|---|
| Figure 2 | A | `Figure02A_cohort_maps.R` | `master_cell_annotation.csv` | `add_celltype_level1.py`, `build_master_annotation.py` |
| Figure 2 | B | `Figure02B_composition_coda.R` | `fig01gh_coda_global.csv`, `fig01gh_coda_persub.csv` | `prep_fig01gh_coda.py` |
| Figure 2 | C | `Figure02C_cellstate_density.R` | `fig01l_tme_cellstate_density.csv`, `fig01l_tme_cellstate_stats.csv` | `prep_fig01kl_tme_density.py` |
| Figure 3 | A | `Figure03A_myoep_oxtr_limitation.R` | `fig03h_auc.csv` | `prep_fig03h_basal_vs_myoep.py` |
| Figure 3 | B | immunohistochemistry image, no script | — | — |
| Figure 3 | C | `Figure03C_fenv_distribution.R` | `fig04b_pool_boundary_cells.csv`, `fig04b_pool_tau_sensitivity.csv` | `prep_fig04b_pool_boundary.py` |
| Figure 3 | D | `Figure03D_boundary_myoep_dcis_vs_mdcis.R` | `fig04d_boundary_myoep_percore.csv`, `fig04d_boundary_myoep_stats.csv` | `prep_fig04d_boundary_myoep.py` |
| Figure 3 | E | `Figure03E_boundary_myoep_maps.R` | `fig04b_boundary_class.csv`, `fig04b_pool_boundary_cells.csv`, `fig04d_boundary_myoep_cells.csv` | `prep_fig04b_connected.py`, `prep_fig04b_pool_boundary.py`, `prep_fig04d_boundary_myoep.py` |
| Figure 3 | F | `Figure03F_boundary_deficient_counts.R` | `fig04p_boundary_counts_percore.csv`, `fig04p_boundary_counts_stats.csv` | `prep_fig04p_boundary_counts.py` |
| Figure 5 | A | `Figure05A_macrophage_m1m2.R` | `fig01j_cellstate_percore.csv` | `prep_fig01j_cellstate_percore.py` |
| Figure 5 | B | `Figure05B_lr_dcis_vs_mdcis.R` | `fig14_lr_stats.csv` | `prep_fig14_ccc.py` |
| Figure 5 | C | `Figure05C1_tcell_engage_position.R`, `Figure05C2_tcell_engage_pooled.R` | `fig14i_tcell_allcores.csv` | `prep_fig14i_tcell_allcores.py` |
| Figure 5 | D | `Figure05D_tcell_checkpoint_allcores.R` | `fig14i_tcell_allcores.csv` | `prep_fig14i_tcell_allcores.py` |
| Figure S1 | A–D | `FigureS1A_tumor_pool_umap.R`, `FigureS1B_tumor_pool_program_scores.R`, `FigureS1C_tumor_pool_core_dominance.R`, `FigureS1D_tumor_pool_excluded.R` | `cl3_4_10_combined_program_scores.csv`, `figS1_tumor_pool_core_dominance.csv`, `figS1a_tumor_pool_umap.csv`, `figS2_tumor_pool_subclusters.csv` | `prep_figS1S2_tumor_pool.py`, `prep_figS1a_tumor_pool_umap.py` |
| Figure S2 | A–D | `FigureS2_annotation_panels.R` | `figS2a_umap.csv`, `figS2c_major_markers.csv`, `figS2d_major_quality.csv` | `prep_figS2_annotation.py` |
| Figure S3 | — | `FigureS3_annotation_flow.R` | `figS3_annotation_flow.csv` | `prep_figS3_annotation_flow.py` |
| Figure S4 | A–G | `FigureS4_myoep_identification.R` | `fig03h_auc.csv`, `fig03h_basal_myoep_dotplot.csv`, `fig03m_basal_subcells.csv`, `figS4_boundary_distance_cells.csv`, `figS4_boundary_distance_percore.csv`, `figS4_group_expression.csv`, `figS5_basal_oxtr_accounting.csv`, `figS5_basal_subcluster_summary.csv` | `prep_fig03h_basal_vs_myoep.py`, `prep_fig03m_basal_subcluster.py`, `prep_figS4_boundary_distance.py`, `prep_figS5_basal_myoep_rule.py` |
| Figure S5 | A–D | `FigureS5_S8_niche_panels.R` | `fig27_cells_example.csv`, `fig27_coda.csv`, `fig27_coda_global.csv`, `fig27_niche_composition.csv`, `fig27_percore.csv`, `rev_cellularity_percore.csv`, `rev_cellularity_stats.csv` | `prep_cellularity_boundary.py`, `prep_fig02e_niche_percore_clr.py`, `prep_fig27_niche.py` |
| Figure S6 | A–B | `FigureS6_spillover_bysource.R` | `fig24_spillover_bysource.csv` | `prep_fig24_spillover_bysource.py` |
| Figure S7 | A–B | `FigureS7A_composition_coda.R`, `FigureS7B_composition_percore.R` | `fig01ef_coda_global.csv`, `fig01ef_coda_percore.csv`, `fig01ef_coda_persub.csv` | `prep_fig01ef_coda.py` |
| Figure S8 | — | `FigureS5_S8_niche_panels.R` | `fig27_cells_example.csv`, `fig27_coda.csv`, `fig27_coda_global.csv`, `fig27_niche_composition.csv`, `fig27_percore.csv`, `rev_cellularity_percore.csv`, `rev_cellularity_stats.csv` | `prep_cellularity_boundary.py`, `prep_fig02e_niche_percore_clr.py`, `prep_fig27_niche.py` |
| Figure S9 | A–C | `FigureS9A_boundary_schematic.R`, `FigureS9B_boundary_composition.R`, `FigureS9C_coverage_composition.R` | `fig04d_boundary_myoep_overall.csv`, `fig04d_boundary_myoep_percore.csv`, `fig04d_boundary_myoep_stats.csv` | `prep_fig04d_boundary_myoep.py` |
| Figure S10 | — | `FigureS10_deficient_neighborhood.R` | `fig04n_deficient_nbr_percore.csv`, `fig04n_deficient_nbr_stats.csv` | `prep_fig04n_deficient_neighborhood.py` |
| Figure S11 | A–B | `FigureS11_deg_across_cellsets.R` | `figS11_deg_genes.csv`, `figS11_deg_matrix.csv` | `prep_figS11_deg_across_cellsets.py` |

Preparation scripts keep their working names because several of them serve more than one
figure (`prep_fig04d_boundary_myoep.py` feeds Figure 3D, 3E and S9, for instance) and many
produce supplementary tables that belong to no numbered figure.

Supplementary tables (differential expression, compositional analysis, spillover,
cell-cell communication, QC and the annotation dictionaries) are the remaining CSVs in
`11.tables_publication/`; each is written by the `prep_*.py` script sharing its figure prefix.

## Patient identifiers

`pathology` in `02.slide_tma_core/mibc_tma_map.csv` and `03.data_processed/tma_core_qc.csv`
is a pseudonymous patient code (`P01`-`P41`), not the original pathology accession number.
Each core comes from a distinct patient, so the codes carry the same grouping information
as the accessions they replace.

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

```bash
# 1. from GEO downloads to the object the clustering starts from
python scripts/load_from_geo.py --geo-dir /path/to/GSE343808
#    -> 03.data_processed/geo_slides.h5ad

# 2. initial clustering: QC-passing cores, normalization, HVG, PCA, Harmony, Leiden
python scripts/cluster_from_geo.py
#    -> 03.data_processed/integrated_qc_passed_from_geo.h5ad
```

This path was checked end to end against the deposit: the six files assemble to
546,021 cells x 374 genes, and after the core QC table and the cell filters
(`min_counts=10`, `min_genes=5`) the result is exactly the 520,506 cells of the published
object — same cell identifiers, same raw counts.

Core-level QC is not recomputed by step 2: it uses the per-slide graph-based clustering
that Xenium writes beside the raw bundle, which is not part of the GEO deposit. Its outcome
is recorded in `03.data_processed/tma_core_qc.csv` (39 of 41 cores pass) and is applied
directly, so the cell set matches the published one.

`integrate_slides_harmony.py` and `integrate_qc_passed.py` are the scripts that produced
the published clustering from the original Xenium bundles; `load_from_geo.py` +
`cluster_from_geo.py` are the same steps with the same parameters, entered from the GEO
files. From there, `subcluster_*.py`, `merge_subclusters.py`, `build_master_annotation.py`
and the `prep_*.py` scripts run unchanged.

### Without downloading anything

Two summed count matrices travel with the code, written by
`scripts/build_pseudobulk_counts.py`:

* `04.tables/pseudobulk_core_gene_counts.csv` — 39 QC-passing cores x 374 genes.
* `04.tables/pseudobulk_core_celltype_gene_counts.csv` — core x major cell type x 374 genes.

With `04.tables/master_cell_annotation.csv.gz` (per-cell type, state, core and spatial
coordinates) these reproduce the differential expression, compositional and per-core
analyses — every pseudobulk result in the paper — from this repository alone. Analyses that
need single-cell expression (clustering, programme scoring, neighbourhood and spillover
measures) need the GEO download.

## Data not included here

* The cell-level count matrix and the raw Xenium output are at GEO (GSE343808), as above.
* The immunohistochemistry images used for Figure 3B (p63, SMMHC) are not redistributed here.

## Licence

Code (`scripts/`, `21.figures_hires_scripts/`, `scripts_r/`) is released under the MIT
License (`LICENSE`). Data files and documents are released under CC BY 4.0
(`LICENSE-DATA`). Please cite the accompanying publication when you use either.
