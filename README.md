# NK Cells in Tumor Archetype: Analysis Pipeline

This repository contains the analysis pipeline for investigating NK (Natural Killer) cells within tumor archetypes. The analysis focuses on identifying distinct NK cell subtypes, their markers, and their functional relevance in tumor microenvironments.

---

## Project Overview

The project aims to:
1. Extract and characterize NK cell markers from bulk data.
2. Refine clustering of NK subtypes and identify unique transcriptional signatures.
3. Compare NK subtypes to other immune cells to identify distinct features.
4. Visualize overlaps in transcriptional markers across NK subtypes using UpSet plots.

---

## Workflow Summary

### 1. **Data Preprocessing**
- Load the CITE-seq dataset and preprocess the Seurat object.
- Filter for specific time points (e.g., pre-injection) and normalize RNA and ADT assays using SCTransform and CLR normalization methods.
- Remove contaminating and proliferating NK cells for cleaner downstream analysis.

### 2. **Marker Discovery**
- **Global NK Markers:** Identify markers for NK cells compared to all immune cells (`Table1_All_Markers_NK`).
- **Subtype-Specific Markers:** Identify markers for individual NK subtypes (NK1, NK2, NK3):
  - `Table5`: Markers for NK1 vs. NK2 + NK3.
  - `Table6`: Markers for NK2 vs. NK1 + NK3.
  - `Table7`: Markers for NK3 vs. NK1 + NK2.
- Save marker tables in Excel format for downstream analysis.

### 3. **NK Subtype Refinement**
- Isolate and re-cluster NK cells using weighted nearest neighbor (WNN) analysis with RNA (PCA) and ADT (apca) modalities.
- Generate UMAP visualizations to assess clustering outcomes and remove residual contaminations.
- Assign biologically relevant cluster names (`NK1`, `NK2`, `NK3`) and save pre-processed objects.

### 4. **Comparative Analysis**
- Compare NK subtypes against all non-NK immune cells:
  - `Table2`: NK1 vs. all non-NK cells.
  - `Table3`: NK2 vs. all non-NK cells.
  - `Table4`: NK3 vs. all non-NK cells.
- Use the Seurat `FindMarkers` function with a Wilcoxon test for statistical evaluation.

### 5. **Visualization**
- Generate an UpSet plot to highlight shared and distinct markers across NK subtypes using the top 40 markers for `Table5`, `Table6`, and `Table7`.

---

## Parameters and Settings

### Global Parameters

| Parameter                | Value                       | Description                           |
|--------------------------|-----------------------------|---------------------------------------|
| `SEED`                   | `42`                       | Random seed for reproducibility       |
| `VARIABLE_FEATURES_MAXNB`| `2000`                     | Max variable features for PCA         |
| `FINDMARKERS_METHOD`     | `"wilcox"`                 | Method for marker identification      |
| `FINDMARKERS_MINPCT`     | `0.2`                      | Minimum fraction of cells for testing |
| `FINDMARKERS_LOGFC_THR`  | `0.25`                     | Log2 fold-change threshold            |
| `FINDMARKERS_PVAL_THR`   | `0.001`                    | P-value threshold for significance    |

### Clustering Settings

| Parameter               | Value      | Description                          |
|-------------------------|------------|--------------------------------------|
| `FINDCLUSTERS_RESOLUTION`| `0.5`     | Resolution for clustering            |
| `FINDCLUSTERS_DIMS`     | `1:20`     | PCA dimensions for clustering        |
| `FINDNEIGHBORS_K`       | `30`       | Number of neighbors for WNN graph    |

---

## Outputs

1. **Marker Tables**:
   - Ordered marker tables (Excel format) for all comparisons:
     - NK vs. immune cells (`Table1`)
     - Subtype-specific markers (`Table5`, `Table6`, `Table7`)
     - NK subtypes vs. non-NK cells (`Table2`, `Table3`, `Table4`).

2. **Visualization**:
   - UMAP plots showing NK subtype clusters.
   - UpSet plots of overlapping markers across NK subtypes.

3. **Pre-Processed Objects**:
   - Saved Seurat objects for pre-processed and refined NK cells:
     - `Seurat_NKonly_MetaNK_Preprocessed.rds`.



## How to Use

1. Clone the repository and install dependencies.
2. Place the input datasets (`SAMPLE_ALL_CITEseq`, `SAMPLE_NK_CITEseq`) in the appropriate directory.
3. Run scripts sequentially for:
   - Preprocessing (`preprocessing.R`).
   - Marker discovery (`markers_discovery.R`).
   - Visualization (`visualization.R`).

4. Results are saved to the `PATH_ANALYSIS_OUTPUT` directory.

---

This repository provides a structured pipeline for the analysis of NK cells in tumor archetypes, facilitating the identification of biologically relevant subtypes and their signatures. For any issues or questions, please contact the repository maintainers.

