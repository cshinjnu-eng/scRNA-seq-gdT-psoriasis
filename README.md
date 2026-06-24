# CSMD1_gdT

Single-cell RNA-seq analysis of **γδ T cells and CSMD1 in psoriasis**.

This repository holds the complete analysis code accompanying the study (manuscript under
revision). It traces the full path from raw 10x matrices to the cell–cell communication,
trajectory, and gene-regulatory analyses that characterise a **GZMK⁺CSMD1⁺ γδ T** population.

## Data

- **Input:** 6 samples, 10x Genomics gene-expression matrices —
  `NPBMC_gdt1–3` and `PBMC_gdt1–3` (control vs. psoriasis PBMC γδ T sorts; see §1 of the script
  for the exact grouping).
- **Annotated object:** `Overall_shinanxi_New.rds` — a Seurat (v4) object,
  **67,742 cells in 8 γδ T subclusters** (including a Naive-like root and a GZMK⁺CSMD1⁺
  population). At **688 MB** it exceeds GitHub's 100 MB file limit, so it is distributed as a
  **Release asset** rather than tracked in Git — see the [Releases](../../releases) page.

## Analysis pipeline — `Rtext_Psoriasis_CSMD1.R`

A single working script (≈16,700 lines) organised into 19 sections:

| # | Section | Tooling |
|---|---|---|
| 1 | Initialize & read 10x data | Seurat v4 |
| 2 | Doublet removal | DoubletFinder |
| 3 | Standard quality control | Seurat |
| 4 | Dimensionality reduction & clustering | Seurat (PCA/UMAP) |
| 5 | Batch-effect removal | FastMNN (SeuratWrappers / batchelor) |
| 6 | Cell-type identification | SingleR + manual markers |
| 7 | Cell-cycle analysis | Seurat |
| 8 | Basic map plotting | scCustomize, ggplot2 |
| 9 | Pathway enrichment | clusterProfiler, org.Hs.eg.db, GO/DOSE |
| 10 | Cell–cell communication | CellChat |
| 11 | CellChat multi-group comparison | CellChat |
| 12 | CSMD1 subgroups via public databases | CellChat |
| 13 | Pseudotime | Monocle 2 |
| 14 | Developmental potential | CytoTRACE |
| 15 | Developmental potential | CytoTRACE 2 |
| 16 | Switch-gene analysis | GeneSwitches |
| 17 | Regulon / TF activity | SCENIC |
| 18 | Metabolite-mediated communication | MEBOCOST (Python) |
| 19 | In-silico perturbation | CellOracle (Python) |

## Dependencies

**R (Seurat v4 stack):** Seurat, SeuratWrappers, batchelor, harmony, DoubletFinder,
SingleR, celldex, clusterProfiler, org.Hs.eg.db, GO.db, DOSE, CellChat, monocle,
scCustomize, patchwork, tidyverse, pheatmap, RColorBrewer.

**Python:** MEBOCOST (§18), CellOracle (§19), CytoTRACE2 (§15) — run these chunks in a
Python environment, not in R.

## Notes

- The script keeps the **original hard-coded paths** (e.g. `setwd("F:\\03_Rtest\\...")`,
  `Read10X(...)`); adapt them to your own filesystem before running.
- It is a **mixed R + Python** working record — R for §1–17, Python for §18–19. Execute each
  block in the appropriate interpreter.
- Seurat is pinned to v4 (the script installs into a local `SeuratV4` libpath).

## Citation

To be added upon publication.
