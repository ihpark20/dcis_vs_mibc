# Gene panel

`panel_gene_table.csv` describes the 374 genes measured in this study, one row per gene:

| column | meaning |
|---|---|
| `gene`, `ensembl` | gene symbol and Ensembl ID |
| `source` | `Breast` (200), `Breast+IO` (80 shared by both panels), `IO` (94 added) |
| `annotation_breast`, `annotation_io` | the 10x panel annotation for the gene |
| `programs`, `cell_lineages`, `program_states`, `marker_states` | programme signatures the gene belongs to in the pool subclustering (170 of 374 genes) |
| `pct_expressing`, `mean_lognorm`, `mean_counts` | expression across the analysed cells |
| `top_major`, `top3_majors`, `top_compartment`, `mean_*` | the cell types and compartments where the gene is highest |

The panel is the pre-designed 10x Xenium Human Breast panel (`hBreast_v1.1`, 280 targets)
plus a 94-target custom immuno-oncology add-on (`hBreast_94g`, design ID `TEKTBT`), giving
374 measured genes. The vendor panel definitions are distributed by 10x Genomics and are
not duplicated here; each gene's panel of origin is in the `source` column.

This table is Supplementary Table 1 of the paper.
