#Dependencies
library(SeuratDisk)

reference <- LoadH5Seurat(file.path(PATH_EXPERIMENT_REFERENCE , SAMPLE_ALL_CITEseq))
reference@assays[["ADT"]] = NULL #Remove ADT to free some memory

#Subset and re-Normalize
reference = SetIdent(reference,  value= "celltype.l1")
reference = subset(reference, subset= time== "0")
reference = SCTransform(reference, assay="SCT")           
        
#Find the markers
Table1_All_Markers_NK = FindMarkers(reference, ident.1 = "NK", only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE )

# Order the data frame by avg_log2FC in descending order
# Include rownames as a column and order by avg_log2FC
Table1_All_Markers_NK <- Table1_All_Markers_NK %>%
  tibble::rownames_to_column(var = "Gene") %>%
  arrange(desc(avg_log2FC))

#Have a quick look at it
Table1_All_Markers_NK

# Save the ordered data frame as an Excel file
output_file <- file.path(PATH_ANALYSIS_OUTPUT, "Ordered_Table1_Markers_NKvsAllImmuneCells.xlsx")
write.xlsx(Table1_All_Markers_NK, output_file, rowNames = FALSE)


