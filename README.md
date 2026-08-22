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
Harmony integrates them and Leiden partitions the result into the fifteen clusters
(CL0-CL14) the annotation is built on. There are two ways to get there.

### Either: take the clustering the paper used

```bash
python scripts/cluster_from_geo.py --labels 03.data_processed/integrated_cluster_labels.csv.gz
#    -> 03.data_processed/integrated_qc_passed_from_geo.h5ad
```

The published assignment for all 520,506 cells of the analysis travels with this repository
as `03.data_processed/integrated_cluster_labels.csv.gz`. This attaches it and computes
nothing, so it takes seconds instead of an hour and every downstream step starts from
exactly the partition the paper reports. **This is the path the rest of this repository
uses.** Clustering is worth rerunning to check that it reproduces, but not as the routine
entry point to everything else.

### Or: cluster it yourself

```bash
python scripts/cluster_from_geo.py
#    -> 03.data_processed/integrated_qc_passed_from_geo.h5ad
#    -> 03.data_processed/cluster_naming.csv
```

Around an hour on half a million cells, the Leiden step alone taking a good part of it on a
single core. It keeps the cores flagged `analysis_include` in the QC table, then follows the
published settings:

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

Expect a partition close to the published one but not identical to it. In our rerun the
same 520,506 cells came out in fifteen clusters that matched one-to-one, with 92 % of cells
landing in the same one; the epithelial clusters differed by well under a percent, while
cells moved between the three CAF clusters and between the three myeloid ones, which look
alike and share a pool anyway. That is what makes the numbering step below necessary.

### Your cluster numbers will not be the manuscript's

Leiden numbers clusters by size, so a rerun renumbers everything: in our own rerun Leiden's
cluster 1 was the manuscript's CL6, and its cluster 2 was CL1. Read a figure by the raw
Leiden number and you are looking at the wrong population. The numbers therefore have to be
matched to the manuscript rather than trusted. **If you cluster yourself, that matching is
a step you have to run** — `cluster_from_geo.py` does it before writing its output and
records the result in `03.data_processed/cluster_naming.csv`, which everything downstream
depends on. Taking the published labels skips the problem, since they arrive already named.

Each cluster is summarised by its mean expression across the panel and correlated against
the reference profile of every published cluster in
`03.data_processed/integrated_cluster_profiles.csv`; the one-to-one assignment with the
highest total correlation wins. The result goes to `03.data_processed/cluster_naming.csv`,
one row per cluster:

| leiden | cluster | correlation | n_cells |
|---|---|---|---|
| 1 | CL6 | 1.000 | 57867 |
| 2 | CL1 | 1.000 | 57012 |

`obs["cluster"]` then carries the manuscript's names and `obs["leiden"]` keeps the raw ones,
so anything downstream can speak in CL numbers. If you cluster some other way — a different
resolution, another tool — match your clusters to the same reference profiles and you can
still speak in CL numbers. Check the correlation column: ours ran above 0.98 for fourteen of
the fifteen clusters, and a low value means the cluster you found does not correspond
cleanly to any published one and should not borrow its name.

Matching this way is far steadier than the clustering itself, because a cluster's average
profile barely moves when a few per cent of its boundary cells change hands.

What each CL was annotated as is in
`03.data_processed/integrated_cluster_annotation_manuscript.csv`, one row per cluster: the
annotation, the pool it was sent to, and the three lines of evidence behind the call — the
cell types those cells carried in the per-slide annotation, the compartment composition,
and the top markers. Three clusters have a note explaining why the marker evidence
overruled the prior label.

## Subclustering of Pooled Clusters

One clustering of half a million cells over 374 genes separates lineages but not the states
inside them: a small lineage, or a state that differs in a handful of genes, is absorbed by
the dominant axis of variation. The fifteen clusters are therefore grouped into seven
lineage pools and each pool is re-clustered on its own markers. As with the clustering,
there are two ways to get there.

### Either: take the cell states the paper assigned

```bash
python scripts/use_published_subclusters.py
#    -> 03.data_processed/subclustered/cell_states.csv
```

The published assignment for all 477,681 cells travels with this repository as
`03.data_processed/subcluster_labels.csv.gz` — pool, subcluster, cell state and the
fragment flag for each cell. Seconds instead of an hour, and the cell states are exactly
the paper's.

### Or: run the paper's subclustering

```bash
python scripts/prepare_pool_input.py            # the object those scripts expect
python scripts/paper_subcluster_cl0_tcell_only.py
python scripts/paper_subcluster_cl1_endothelial.py
python scripts/paper_subcluster_cl2_11_12_myeloid_dc.py
python scripts/paper_subcluster_cl3_4_10_combined.py
python scripts/paper_subcluster_cl5_8_bplasma.py
python scripts/paper_subcluster_cl6_7_13_14_caf_adipo.py
python scripts/paper_subcluster_cl9_mast.py
python scripts/paper_merge_subclusters.py       # fragment and contamination flags, merge
python scripts/paper_build_master_annotation.py # one row per cell
python scripts/paper_add_celltype_level1.py     # major cell types
```

The `paper_*` scripts are the study's own code, one per pool, changed only in their path
constants — the filters, features, parameters, seeds, programme scoring, subtype rule,
fragment threshold and merge are as published. Around an hour for all seven pools.

There is also `subcluster_from_geo.py`, a single consolidated implementation of the same
recipe with the per-pool settings gathered into `pool_config.py`. It is easier to read and
to modify, and it is what the rest of this README describes; the `paper_*` scripts are the
reference it was written against.

### Which cluster goes into which pool

The paper groups its clusters by hand: CL3, CL4 and CL10 make the tumor pool, CL2, CL11 and
CL12 the myeloid one, and so on. Those groupings are in `scripts/pool_config.py`, and
because the clusters have already been matched to the published numbering they can be
applied as they stand — which is what makes these the paper's pools rather than something
rederived.

```bash
python scripts/assign_pools_from_geo.py
#    -> 03.data_processed/pool_assignment.csv
```

The assignment is also rederived from the data as a check: each cluster is scored against
the lineage genes of all seven pools, and the best-scoring pool is recorded next to the
assigned one. Agreement means the cluster matching held. Disagreement means a cluster did
not match its published counterpart cleanly and is being sent somewhere its expression does
not support — the script prints those rows rather than letting them pass. In our rerun all
fifteen agreed.

`subcluster_from_geo.py` reads that file; each pool then follows the same recipe.

TMA1 D6 is dropped before the pools are formed. The clustering itself ran on every
QC-passing core, D6 among them, and the core was set aside afterwards; the subclustering
was then redone without it, which is the version everything downstream uses. The exclusion
lives in `scripts/analysis_exclusions.py`, and the boundary and communication steps read it
from there too.


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
   subcluster takes the name of its highest-scoring one — chosen among the pool's lineage
   programmes only. A pan-lineage or interferon signature describes something every
   subcluster of the pool carries to some degree, so letting it compete lets it claim
   subclusters a specific programme should name; `pool_config.py` keeps those out of the
   running, as the original analysis does. It matters: with `core_macrophage` in the
   running the macrophage subclusters stop being M1 or M2, and with `stress_g1arrest` in it
   a quarter of the tumor pool is named after a single cell-cycle gene.

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
markedly less covered. The median core has 58 % of its boundary myoepithelium deficient in
DCIS against 82 % in mDCIS, where the paper has 66 % and 86 %, and sheath and lined shift
the other way by a similar margin. Three of the four measures survive correction here and
the fourth sits just outside it, against all four in the paper.

Coverage comes out a few points higher throughout because this rerun calls more cells
myoepithelial than the paper does — the one major cell type still well off, and the reason
the deficient share is 67 % overall against 74 %. The difference between DCIS and mDCIS is
what carries; the level does not, and with 38 cores the effect size travels between runs
more reliably than the significance verdict.

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

## Macrophage polarization

Macrophages run between an inflammatory M1 state and a tissue-remodelling, immunosuppressive
M2 one, and which way a lesion's macrophages lean says something about the environment it is
building. Each macrophage subcluster is scored for both programmes and takes the higher one;
the assignment is made per subcluster rather than per cell, a subcluster being a group the
clustering already found coherent and an average over it steadier than a handful of sparse
markers read off one cell.

| programme | markers |
|---|---|
| M1 | CD86, CD80, IL1B, TNF |
| M2 | CD163, MRC1, MARCO, ARG1, IL10, TREM2 |

```bash
python scripts/macrophage_polarization_from_geo.py
#    -> 03.data_processed/macrophage/polarization_{cells,percore,stats}.csv
```

Each core contributes the share of its macrophages on each side, their ratio, and how many
of each it holds — per thousand cells and per square millimetre of tissue, area being what
the core QC measured as actually covered by cells rather than the whole core. Cores are
compared by Mann-Whitney with Benjamini-Hochberg across the measures.

The directions reproduce: DCIS cores lean more M1, microinvasive cores more M2, and
microinvasive cores carry more M2 macrophages per thousand cells. Densities per square
millimetre separate the groups least of all — a microinvasive core holds more macrophages
but is also more cellular, so the difference in composition does not survive being expressed
per unit of tissue.

None of it is significant, in this rerun or in the paper, where the same comparisons sat at
p = 0.07 to 0.08 before correction. The closest is M2 per thousand cells, 44.5 in DCIS
against 56.0 in mDCIS at p = 0.09, where the paper had 38.3 against 48.0 at p = 0.08. The
balance point still sits higher than the paper's — a third of macrophages called M1 against
a fifth — which is what remains of the subcluster labelling once the pan-macrophage
programme is out of the running.

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

How far the divergence travels depends on where you start. Beginning from the published
clustering, all seven pools hold exactly the cells they hold in the paper — 108,594 in the
tumor pool, 96,841 in T/NK, and so on — 98 % of cells land in the same major cell type, and
the ligand-receptor enrichments come out identical to three decimals, that analysis reading
only coordinates and counts and so being fully determined once the cell set is. What still
moves is inside the pools: Leiden splits each one a little differently, which shifts where
the myoepithelial call lands and carries a few percent into the boundary composition.
Beginning from a fresh clustering instead, the major cell types agree for 87 % of cells,
the CAF and myeloid clusters trading members. Both runs give the same answer to the
biological question; they differ in how many cells sit on either side of it.

## Data not included here

* The cell-level count matrix and the raw Xenium output are at GEO (GSE343808), as above.
* The immunohistochemistry images used for Figure 3B (p63, SMMHC) are not redistributed here.

