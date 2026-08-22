# Spatial transcriptomics of DCIS and microinvasive breast carcinoma

Code and data to reproduce the findings in '[Spatial transcriptomics reveals tumor-stroma interface remodeling in HER2-positive ductal carcinoma in situ with microinvasion]' (under review). 
Two breast tissue microarrays (41 cores) were profiled with 10x Genomics Xenium using a
374-gene custom panel (280-gene Human Breast panel + 94-gene immuno-oncology add-on) and
compared between pure DCIS and microinvasive carcinoma (mDCIS).

## Which figure each analysis is

| panel | what it shows | script | output |
|---|---|---|---|
| Fig 2B | microenvironment composition, DCIS vs mDCIS, by CoDA | `tme_composition_from_geo.py` | `tme/composition_*.csv` |
| Fig 2C | cell-state density per unit of microenvironment tissue | `tme_composition_from_geo.py` | `tme/density_cellstate.csv` |
| Fig 3C | how much microenvironment each epithelial cell touches (`f_env`) | `boundary_from_geo.py` | `boundary/pool_boundary_cells.csv`, `tau_sensitivity.csv` |
| Fig 3D | myoepithelial coverage of the boundary, DCIS vs mDCIS | `boundary_myoep_from_geo.py` | `boundary/boundary_myoep_{percore,stats}.csv` |
| Fig 3E | the same, mapped over two example cores | `boundary_myoep_from_geo.py` | `boundary/boundary_myoep_cells.csv` |
| Fig 3F | deficient boundary per core, absolute and per 100 tumor cells | `boundary_myoep_from_geo.py` | `boundary/boundary_counts_{percore,stats}.csv` |
| Fig 5A | macrophage M1/M2 by core | `macrophage_polarization_from_geo.py` | `macrophage/polarization_*.csv` |
| Fig 5B | ligand-receptor proximity, DCIS vs mDCIS | `ccc_from_geo.py` | `ccc/lr_{enrich,stats}.csv` |
| Fig 5C | where PD-1+ T cells engage PD-L1 | `checkpoint_engagement_from_geo.py` | `ccc/checkpoint_*.csv` |
| Fig S6 | transcript spillover, by source cell type | `spillover_from_geo.py` | `spillover/spillover_by_source.csv` |

Figure 2A is a map of the annotated cells, and Figure 3B an immunohistochemistry image;
neither is an analysis. Figure 3A rests on a myoepithelial marker comparison that is not
reproduced here.

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

```bash
python scripts/attach_published_clusters.py
#    -> 03.data_processed/integrated_clusters.h5ad
```

The assignment for all 520,506 cells of the analysis travels with this repository as
`03.data_processed/integrated_cluster_labels.csv.gz`, and this attaches it to the cells
built from GEO. Nothing is recomputed.

That is deliberate. Reclustering takes about an hour and does not land on the published
partition: Harmony and Leiden accumulate floating point in an order that depends on the
number of threads and the BLAS build, so cells near a cluster boundary fall either way, and
Leiden renumbers its clusters by size, so CL3 would no longer name the population it names
in the paper. Every step downstream would then rest on an approximation of the paper rather
than on the paper. The settings that produced the clustering were: cells filtered at
`min_counts=10` and `min_genes=5`, counts normalised to 10,000 and log1p-transformed, 300
highly variable genes selected per slide, 30 principal components, Harmony over slide, a
15-neighbour graph and Leiden at resolution 0.5, all seeded at 42.

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
lineage pools — CL3, CL4 and CL10 make the tumor pool, CL2, CL11 and CL12 the myeloid one,
and so on — and each pool was re-clustered on its own lineage markers, in place of the
highly variable genes that a single lineage's variation is dominated by.

```bash
python scripts/use_published_subclusters.py
#    -> 03.data_processed/subclustered/cell_states.csv, major_celltypes.csv
```

The published assignment for all 477,681 cells travels with this repository as
`03.data_processed/subcluster_labels.csv.gz` — pool, subcluster, cell state, fragment flag
and major cell type for each cell. As with the clustering, it is attached rather than
recomputed, for the same reason.

## What the microenvironment is made of (Fig 2B, 2C)

Two questions about the same cells, kept apart because they can disagree: does the mix of
cell types differ between DCIS and microinvasive cores, and is any type packed more densely
in the tissue?

```bash
python scripts/tme_composition_from_geo.py
#    -> 03.data_processed/tme/composition_{global,pertype,percore}.csv
#    -> 03.data_processed/tme/density_{pertype,cellstate}.csv
```

Tumor and myoepithelial cells are set aside and the composition renormalised within the ten
immune and stromal types, so the question is about the microenvironment and not about how
much tumor a core happens to contain. Neutrophils are out as well, the panel being unable to
resolve them.

Proportions are compositional: they sum to one, so a rise in any one type forces the others
down, and testing them as independent percentages invents differences that are only the
constraint. Counts are therefore centre-log-ratio transformed against the geometric mean of
the ten, with a pseudocount of 0.5 for the zeros, which frees them; the groups are then
compared as a whole by PERMANOVA on the Aitchison distance over 999 permutations, and type
by type on the CLR coordinates with Benjamini-Hochberg correction.

The composition as a whole does not separate the groups — PERMANOVA gives F = 2.02 at
p = 0.11 — and one type survives correction on its own: dendritic cells, lower in
microinvasive cores at an FDR of 0.005. Reading the ten as raw percentages would have turned
several of the others significant.

Density asks the other question: cells of a type per square millimetre of microenvironment
tissue, where the tissue is the summed area of the microenvironment cells themselves rather
than the area of the core, so a core with open stroma is not counted as sparser. Dendritic
cells fall there too, and by cell state it is cDC2 (261 to 137 per mm2), mregDC (207 to 168)
and pDC (148 to 92) that thin out in microinvasive cores.

## Tumor boundary (Fig 3C)

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

## Myoepithelial coverage of the boundary (Fig 3D-F)

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
markedly less covered. Running on the published cell states, the median core has 66 % of its
boundary myoepithelium deficient in DCIS against 87 % in mDCIS, where the paper has 66 % and
86 %; the overall split is 7.6 % sheath, 17.7 % lined and 74.7 % deficient against 7.9 %,
18.5 % and 73.6 %. All four measures survive correction, as in the paper.

Shares alone can mislead: a core with twice the tumor can hold twice the deficient boundary
at the same percentage. The counts are therefore reported as well, absolute and per 100
tumor cells, and it is the normalised one that separates the groups — 17.1 deficient
boundary cells per 100 tumor cells in the median DCIS core against 32.0 in the median
microinvasive one (p = 0.02), where the raw count does not (499 against 707, p = 0.13). A
microinvasive lesion carries more exposed boundary for the amount of tumor it has, not
merely more tumor.

## Ligand-receptor communication (Fig 5B)

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

### Where PD-1 meets PD-L1 (Fig 5C)

That the checkpoint pairs come out enriched says the cells are close enough to engage; it
does not say where. A PD-1+ T cell held off inside a tumor mass or at its edge is being
stopped at the point of invasion, one engaged out in the stroma is not, so each T cell is
placed as well as classified:

```bash
python scripts/checkpoint_engagement_from_geo.py
#    -> 03.data_processed/ccc/checkpoint_{tcells,percore,stats}.csv
```

A T cell counts as engaged when it carries PDCD1 and a CD274-positive cell lies within
30 um; the position comes from its own surroundings — at least 60 % epithelial neighbours is
inside, 20 to 60 % the boundary, less than that outside.

The engagement moves inward in microinvasive cores. Pooled over all cores, 3.5 % of PD-1+
T cells are engaged inside the tumor and 6.4 % at its boundary in mDCIS, against 1.6 % and
3.5 % in DCIS, with the stromal share falling from 72 % to 65 %.

Cores are still the unit of comparison, and only those with at least 20 PD-1+ T cells take
part: below that the share moves in steps too coarse to compare, a core with five such cells
returning 0, 20 or 40 % and nothing between. That leaves 11 DCIS and 15 microinvasive cores,
where the median share engaged inside or at the boundary is 0 % against 5 % (p = 0.03).
Pooling the cells instead of the cores would make this look far more decisive than it is —
a handful of large cores would carry the result.

## Macrophage polarization (Fig 5A)

Macrophages run between an inflammatory M1 state and a tissue-remodelling, immunosuppressive
M2 one, and which way a lesion's macrophages lean says something about the environment it is
building. The split is the cell state itself: a macrophage is M1 or M2 because the
subclustering named its subcluster `macrophage_M1` or `macrophage_M2`, scoring these
programmes against the rest of the myeloid pool's.

| programme | markers |
|---|---|
| M1 | CD86, CD80, IL1B, TNF |
| M2 | CD163, MRC1, MARCO, ARG1, IL10, TREM2 |

Macrophages the subclustering placed elsewhere — proliferating, monocyte, about a quarter of
them — carry no polarization call. They count towards the macrophage total but towards
neither side, so the shares are of the polarized macrophages.

```bash
python scripts/macrophage_polarization_from_geo.py
#    -> 03.data_processed/macrophage/polarization_{cells,percore,stats}.csv
```

Each core contributes the share of its macrophages on each side, their ratio, and how many
of each it holds — per thousand cells and per square millimetre of tissue, area being what
the core QC measured as actually covered by cells rather than the whole core. Cores are
compared by Mann-Whitney with Benjamini-Hochberg across the measures.

Read this way the comparison reproduces exactly: 18.7 % of polarized macrophages are M1 in
the median DCIS core against 15.0 % in the median microinvasive one, at p = 0.074, and M2
macrophages per thousand cells go 38.3 against 48.0 at p = 0.084 — the paper's numbers to
three decimals. Neither is significant, there or here.

Density per square millimetre, which the paper did not report, separates the groups least of
all: 131 M2 macrophages per mm2 in DCIS against 156 in mDCIS at p = 0.76. A microinvasive
core holds more macrophages but is also more cellular, and the difference in composition
does not survive being expressed per unit of tissue.

## Transcript spillover (Fig S6)

Xenium assigns a transcript to whichever segmented cell it falls in, and segmentation is not
perfect: a transcript from a tumor cell can end up counted in the fibroblast beside it. The
error is not random — it follows whatever the neighbour happens to be — so a gene one cell
type expresses strongly appears in the cells around it and can be read as a real signal
there. This measures how far each cell type's transcripts leak.

```bash
python scripts/spillover_from_geo.py
#    -> 03.data_processed/spillover/spillover_by_source.csv
```

Two things have to hold at once for spillover to be the explanation, so the index is their
geometric mean:

```
source_spec(S, g) = log2( mean of g in S / mean of g outside S )
prox_fc(S, g)     = log2( mean of g in non-S cells surrounded by S
                          / mean of g in non-S cells away from S )
spillover(S, g)   = sqrt(source_spec x prox_fc), when both are positive
```

A gene the source barely expresses cannot leak, however it behaves near the source; a gene
that does not rise near the source is not leaking, however specific it is. Requiring both
keeps genuine expression in neighbouring cells from being mistaken for leakage. A cell
counts as near S when at least half its neighbours within 30 um are S, and as far when
fewer than a tenth are.

What comes out is each lineage's own signature turning up in its neighbours: TPSAB1, CTSG
and CPA3 around mast cells, IGHG1 and MZB1 around plasma cells, OXTR around myoepithelium,
CD3E and IL7R around T cells. The strongest of all are the epithelial genes around tumor —
FOXA1, ABCC11, CEACAM6, ERBB2 — which is why a differential expression result in
non-epithelial cells that reads as an epithelial gene deserves to be checked against this
table before it is believed.

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

Everything from the tumor boundary onward is computed here. What is not recomputed is the
clustering and the subclustering: those are taken from the paper, so that every result rests
on the partition it reports rather than on an approximation of it.

The reason is that clustering does not travel. Harmony and Leiden accumulate floating point
in an order that depends on the number of threads and on the BLAS build, so cells sitting
near a cluster boundary fall either way on another machine, and Leiden numbers its clusters
by size, so a rerun renames them. We did rerun both to see how far that goes: the same
520,506 cells came out in fifteen clusters matching the published ones one to one, with 92 %
of cells in the same one, and the pools built on them agreed on 87 % of major cell type
calls — close, but close is not the paper, and the difference compounds through every step
downstream.

Taking the published assignment instead, what remains is arithmetic on the deposited counts
and coordinates, and it reproduces: the boundary composition lands within a point of the
published one, the macrophage M1/M2 comparison and the ligand-receptor enrichments match to
three decimals.

One step of the original is not reproduced at all. The subclustering screened each subcluster
for cells of a foreign lineage against the per-slide annotation; that judgement is carried in
the published cell states rather than remade here.

Package versions matter for anything that is recomputed — defaults in scanpy and its
dependencies change between releases — and the versions used are listed under Environment.

## Data not included here

* The cell-level count matrix and the raw Xenium output are at GEO (GSE343808), as above.
* The immunohistochemistry images used for Figure 3B (p63, SMMHC) are not redistributed here.

