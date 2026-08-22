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

Cores are judged in two steps, kept apart so that the technical filter never depends on
knowing which cells are tumor.

```bash
python scripts/qc_cores_raw.py          # technical quality, from the GEO files alone
python scripts/qc_cores_composition.py  # what the core is made of
#    -> 02.tma_core_qc/core_qc_raw.csv, tma_core_qc.csv
```

**Stage 1 — was the core measured well?** Cell count, transcripts and genes per cell, and
how densely cells cover the tissue, all read from `cells.parquet` and the core bounding
boxes. A core passes with at least 50 cells, a median of 20 transcripts per cell, and 50
cells per mm2 of occupied tissue — occupied, not total, so an empty rim does not flatter a
torn core. Low-quality fraction, cell area and control-probe fraction are reported without
gating.

**Stage 2 — what is in it?** Composition comes from the graph-based clustering Xenium
produces per slide, mapped to cell types through the annotation of those clusters. The GEO
deposit carries neither, so both travel with this repository as
`02.tma_core_qc/TMA{1,2}.clusters.csv` and `TMA{1,2}.cluster_celltype_annotation.csv`. A
core is flagged when it holds under 3 % tumor cells or fewer than 3 cell types.

Two of the 41 cores are out, one at each stage: **TMA1 D8**, whose cells carry a median of
0 transcripts, and **TMA1 D4**, technically sound but only 0.9 % tumor — one a measurement
failure, the other a core with nothing to contribute to a DCIS versus microinvasion
comparison. `tma_core_qc.csv` records every metric, an `analysis_include` flag and a
one-line `exclude_reason`; the clustering step reads that flag rather than recomputing
anything.

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
`obsm["X_pca_harmony"]`, the UMAP coordinates, Leiden's own labels in `obs["leiden"]` and
the published numbering in `obs["cluster"]`. Pass `--all-cores` to skip the QC filter and
cluster every core instead.

**CL0-CL14 keep their meaning.** Leiden numbers clusters by size, so a rerun renumbers
everything and CL3 would no longer be the same population it is in the paper. Each new
cluster is therefore summarised by its mean expression across the panel, correlated against
the reference profile of every published cluster in
`03.data_processed/integrated_cluster_profiles.csv`, and the one-to-one assignment with the
highest total correlation is taken. The pairing and its correlations are written to
`03.data_processed/cluster_naming.csv`; matching this way is far more stable than the
clustering itself, because a cluster's average profile barely moves when a few per cent of
its boundary cells change hands. What each CL was annotated as is in
`03.data_processed/integrated_cluster_annotation_manuscript.csv`, one row per cluster: the
annotation, the pool it was sent to, and the three lines of evidence behind the call — the
cell types those cells carried in the per-slide annotation, the compartment composition,
and the top markers. Three clusters have a note explaining why the marker evidence
overruled the prior label.

Expect around an hour on half a million cells; the Leiden step alone takes a good part of
it and runs on a single core.

## Subclustering of Pooled Clusters

One clustering of half a million cells over 374 genes separates lineages but not the states
inside them: a small lineage, or a state that differs in a handful of genes, is absorbed by
the dominant axis of variation. The fifteen clusters are therefore grouped into seven
lineage pools and each pool is re-clustered on its own markers.

```bash
python scripts/subcluster_from_geo.py
#    -> 03.data_processed/subclustered/<pool>.h5ad, subcluster_assignments.csv

python scripts/merge_subclusters_from_geo.py
#    -> 03.data_processed/subclustered/cell_states.csv
```

Clusters are assigned to pools by lineage score rather than by number — a rerun numbers its
clusters differently — and each pool then follows the same recipe:

1. **lineage filter.** A cell enters the pool only if it has a non-zero count for at least
   one of the pool's lineage genes, which removes cells placed there on overall similarity
   with no lineage evidence of their own.
2. **curated feature space.** The pool's marker set becomes the feature space for PCA,
   in place of highly variable genes. With 374 panel genes, HVG selection inside a single
   lineage is dominated by the abundant transcripts that leak in from neighbouring cells;
   restricting the features makes the sub-structure drive the components.
3. **embedding and clustering.** Scale, PCA, Harmony over slide, a 15-neighbour graph,
   UMAP, Leiden at the pool's resolution.
4. **programme scoring.** Every cell is scored for the pool's programmes and each
   subcluster takes the name of its highest-scoring one.

Per-pool settings — which clusters, which lineage genes, which features, which programmes,
resolution and number of components — are in `scripts/pool_config.py`, one entry per pool.

`merge_subclusters_from_geo.py` then merges the subclusters by programme and marks those
holding fewer than max(200, 0.5 %) of the pool's cells as `fragment`; fragments keep their
label, being small rather than wrong. Subclusters are not screened for cells of a foreign
lineage — see below.

## Major cell types

The pools give fine-grained states; most analyses need the coarse lineages behind them.

```bash
python scripts/major_celltypes_from_geo.py
#    -> 03.data_processed/subclustered/major_celltypes.csv
```

Each (pool, programme) pair maps to one major type. A pool is usually one lineage, but two
need the programme to decide: the tumor pool holds myoepithelial cells alongside tumor, and
pericytes turn up in both the endothelial and the CAF pool. Plasma cells are split from B
cells, and the dendritic programmes from the macrophages. That gives twelve major types.

## Tumor boundary

Where a duct or a tumor mass meets its surroundings is not marked in the data; it has to be
inferred from who each cell sits next to. The epithelial pool is every cell called Tumor or
Myoepithelial, and for each of its cells `f_env` is the fraction of its neighbours within
30 um that fall outside the pool — stroma, immune cells, vessels.

```bash
python scripts/boundary_from_geo.py
#    -> 03.data_processed/boundary/pool_boundary_cells.csv, tau_sensitivity.csv
```

A cell deep inside an epithelial mass touches only its own kind and has `f_env` near zero —
38 % of the pool sits at exactly zero — while a cell at the edge touches the
microenvironment and the value rises. Two conditions keep the fraction meaningful: at least
3 neighbours, and membership of an epithelial component of at least 30 cells at the same
30 um radius, so that scattered cells are not called a boundary.

`f_env` is continuous and has no natural gap, so the cut is a choice rather than a
discovery. The analysis uses **f_env >= 0.20**, which puts about a third of the pool at the
boundary; `tau_sensitivity.csv` records what every other threshold from 0.05 to 0.50 would
have given, and the decision should be read with that table in hand.

Cores, not subclusters, are excluded at this step: a boundary needs the whole epithelial
compartment of a core to be present.

## Myoepithelial coverage of the boundary

A duct in situ is wrapped in a myoepithelial sheath, and microinvasion is the point at
which tumor cells get past it. Each boundary cell is therefore sorted by what covers it:

| class | meaning |
|---|---|
| `myoep-sheath` | the boundary cell is itself myoepithelial |
| `myoep-lined` | a tumor cell with a myoepithelial cell within 15 um |
| `myoep-deficient` | a tumor cell facing the microenvironment with no myoepithelium |

```bash
python scripts/boundary_myoep_from_geo.py
#    -> 03.data_processed/boundary/boundary_myoep_{cells,percore,overall,stats}.csv
```

The comparison is made across cores, not cells: each core contributes the percentage of its
boundary falling in each class, and DCIS is compared with mDCIS by Mann-Whitney with
Benjamini-Hochberg correction over the four measures. Treating cells as independent would
inflate the sample from tens of cores to tens of thousands of cells and turn any difference
significant.

The result is the one the paper reports — the boundary of a microinvasive lesion is
markedly less covered. In this rerun the median core has 53 % of its boundary myoepithelium
deficient in DCIS against 77 % in mDCIS, with sheath and lined shifted the other way by a
similar margin. The four raw p-values sit between 0.03 and 0.06 and do not survive
correction here, where the published run's did; with 38 cores the effect is stable but the
significance verdict is not, which is worth keeping in mind when rerunning any of the
core-level comparisons.

## Ligand-receptor communication

Two cells can only signal through a short-range ligand if they are close enough to share
the same space, which spatial data can see directly. For each ligand-receptor pair with
both partners on the panel — ten of them — the question is whether receptor-positive cells
sit next to ligand-positive cells more often than the core's own composition would give:

```
obs = ligand-positive neighbours of receptor-positive cells / all their neighbours
exp = ligand-positive cells / all cells in the core
enrichment = log2((obs + 1e-4) / (exp + 1e-4))
```

```bash
python scripts/ccc_from_geo.py
#    -> 03.data_processed/ccc/lr_{enrich,stats,celltypes}.csv
```

Neighbours are cells within 30 um; expression is a raw count above zero. Sums are divided
rather than per-cell fractions averaged, so a sparsely surrounded cell cannot swing the
result, and `exp` is the core's own density, so a core full of ligand-positive cells is not
credited for proximity that comes for free. A core contributes a pair only when it holds at
least 15 cells positive for each partner. Cores are the replicates: each pair is tested
against zero enrichment with a one-sample Wilcoxon corrected across pairs, and DCIS is
compared with mDCIS by Mann-Whitney.

The lymphoid chemokine axes come out strongest — CXCL13-CXCR5 and CCL19-CCR7 both above one
doubling of enrichment — and the checkpoint pairs CD80-CTLA4 and CD274-PDCD1 are the ones
that separate the groups, both higher in microinvasive cores.

What this measures is opportunity, not activity: cells positioned to signal, not signalling
observed.

## What the analysis leaves out

Two exclusions are made downstream of the clustering, and the evidence for both is visible
in the output of this pipeline.

**TMA1 D6.** Two subclusters of the tumor pool, the ones scoring highest for the
luminal-mature programme, draw 79 % and 93 % of their cells from this single core; no other
subcluster comes close. A cell state carried by one core cannot support a comparison
between DCIS and microinvasive disease, so D6 is dropped from the analyses. It stays in the
clustering itself, which is why the cell set still matches the published one.

**Neutrophils.** The panel was not designed to resolve them: of the markers normally used —
FCGR3B, CSF3R, S100A9, ELANE, MPO, CXCR2 — none is on it, and the three that are (S100A8,
CEACAM8, ITGAM) are not specific. No pool recovers a neutrophil state, and the neutrophils
annotated in the original run carry the lowest transcript count of any major type, about a
quarter of the median cell. They are excluded rather than analysed as a population.

## Reproducing this analysis

A rerun reproduces the shape of the published analysis — the same cells, the same lineages,
the same cell states in the same proportions — but not every individual cell label. Four
reasons, none of them avoidable here:

* **Clustering is not portable.** Harmony and Leiden accumulate floating point in an order
  that depends on the number of threads and on the BLAS build. Seeds are fixed, so a rerun
  on the same machine repeats itself, but cells sitting near a cluster boundary can land on
  either side elsewhere — and that difference carries into the pools built on top.
* **Cluster numbering is arbitrary.** The clusters come out in a different order each time,
  which is why pools are matched by lineage score instead of by cluster number.
* **One step is left out.** The original screened each subcluster for cells of a foreign
  lineage, using the per-slide annotation of the initial clustering. A marker-score proxy
  is not specific enough to stand in for it — pericytes and myofibroblasts score like
  endothelium often enough that real cells would be discarded — so the screen is omitted
  and the scores are reported per subcluster instead. (The per-slide clustering itself,
  which the GEO deposit does not carry, is included here for the QC step.)
* **Versions drift.** Defaults in scanpy, leidenalg and umap change between releases; the
  versions this was run with are listed under Environment.

## Data not included here

* The cell-level count matrix and the raw Xenium output are at GEO (GSE343808), as above.
* The immunohistochemistry images used for Figure 3B (p63, SMMHC) are not redistributed here.

