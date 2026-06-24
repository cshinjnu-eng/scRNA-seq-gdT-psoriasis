# scRNA-seq-gdT-psoriasis

Single-cell RNA-seq analysis of peripheral γδ T cells in psoriasis.

Analysis code and processed data accompanying the manuscript (under revision).

## Data availability

- **Input:** six scRNA-seq libraries of sorted peripheral-blood γδ T cells —
  `NPBMC_gdt1–3` (control) and `PBMC_gdt1–3` (psoriasis).
- **Processed object:** `Overall_shinanxi_New.rds` (Seurat v4; 67,742 cells across 8 γδ T
  subclusters), available from the [Releases](../../releases) page.

## Key parameters

| Step                    | Setting                                                         |
| ----------------------- | --------------------------------------------------------------- |
| Species                 | Human                                                           |
| Quality control         | `nFeature_RNA` 200–8,000; `percent.MT` < 10%; `percent.HB` < 3% |
| Doublet removal         | DoubletFinder (expected doublet rate 0.008)                     |
| Feature selection       | 2,000 highly variable genes (vst)                               |
| Batch integration       | FastMNN                                                         |
| Clustering              | Louvain, resolution 0.5                                         |
| Dimensionality          | 15 principal components (UMAP)                                  |
| Cell-type annotation    | SingleR + manual marker curation                                |
| Cell–cell communication | CellChat (min. 10 cells per group)                              |

## Analysis pipeline — `Rtext_Psoriasis.R`

| #   | Section                               | Tooling                                |
| --- | ------------------------------------- | -------------------------------------- |
| 1   | Initialize & read count matrices      | Seurat v4                              |
| 2   | Doublet removal                       | DoubletFinder                          |
| 3   | Standard quality control              | Seurat                                 |
| 4   | Dimensionality reduction & clustering | Seurat (PCA/UMAP)                      |
| 5   | Batch-effect removal                  | FastMNN (SeuratWrappers / batchelor)   |
| 6   | Cell-type identification              | SingleR + manual markers               |
| 7   | Cell-cycle analysis                   | Seurat                                 |
| 8   | Map plotting                          | scCustomize, ggplot2                   |
| 9   | Pathway enrichment                    | clusterProfiler, org.Hs.eg.db, GO/DOSE |
| 10  | Cell–cell communication               | CellChat                               |
| 11  | CellChat multi-group comparison       | CellChat                               |
| 12  | Subgroup analysis via public datasets | CellChat                               |
| 13  | Pseudotime                            | Monocle 2                              |
| 14  | Developmental potential               | CytoTRACE                              |
| 15  | Switch-gene analysis                  | GeneSwitches                           |
| 16  | Regulon / TF activity                 | SCENIC                                 |
| 17  | Metabolite-mediated communication     | MEBOCOST (Python)                      |
| 18  | In-silico perturbation                | CellOracle (Python)                    |

## Software

- **R (Seurat v4):** Seurat, SeuratWrappers, batchelor, harmony, DoubletFinder, SingleR,
  celldex, clusterProfiler, org.Hs.eg.db, CellChat, monocle, scCustomize, tidyverse.
- **Python:** MEBOCOST (§17), CellOracle (§18).

File paths in the script reflect the original analysis environment and should be adapted
before running.

## Citation

Manuscript under revision; citation details will be added upon publication.
