#Refernence loading
reference <- LoadH5Seurat(file.path(PATH_EXPERIMENT_REFERENCE , SAMPLE_ALL_CITEseq))
reference@assays[["ADT"]] = NULL #Remove ADT to free some memory

#Subset and re-Normalize
reference = SetIdent(reference,  value= "celltype.l1")
reference = subset(reference, subset= time== "0")
reference = SCTransform(reference, assay="SCT")      

#Exctract the cell names of all non NK immune cells
reference_not_NK = subset(reference, ident= "NK", invert = TRUE)
Cell_Names_reference_not_NK = Cells(reference_not_NK)

#Extract NK1, NK2 and NK3 cell names
Cell_Names_NK1 = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Cell_Names", "NK1_Cell_Names.rds"))
Cell_Names_NK2 = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Cell_Names", "NK2_Cell_Names.rds"))
Cell_Names_NK3 = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Cell_Names",  "NK3_Cell_Names.rds"))


#Table 2
#Find the markers
Table2_NK1vsAllnonNK = FindMarkers(reference, ident.1 =  Cell_Names_NK1 ,  ident.2 = Cell_Names_reference_not_NK, only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE )

# Order the data frame by avg_log2FC in descending order
# Include rownames as a column and order by avg_log2FC
Table2_NK1vsAllnonNK <- Table2_NK1vsAllnonNK %>%
  tibble::rownames_to_column(var = "Gene") %>%
  arrange(desc(avg_log2FC))


#Have a quick look at it
Table2_NK1vsAllnonNK

# Save the ordered data frame as an Excel file
output_file <- file.path(PATH_ANALYSIS_OUTPUT, "Ordered_Table2_Markers_Table2_NK1vsAllnonNK.xlsx")
write.xlsx(Table2_NK1vsAllnonNK, output_file, rowNames = FALSE)



#Table 3
#Find the markers
Table3_NK2vsAllnonNK = FindMarkers(reference, ident.1 =  Cell_Names_NK2 ,  ident.2 = Cell_Names_reference_not_NK, only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE )

# Order the data frame by avg_log2FC in descending order
# Include rownames as a column and order by avg_log2FC
Table3_NK2vsAllnonNK <- Table3_NK2vsAllnonNK %>%
  tibble::rownames_to_column(var = "Gene") %>%
  arrange(desc(avg_log2FC))

#Have a quick look at it
Table3_NK2vsAllnonNK

# Save the ordered data frame as an Excel file
output_file <- file.path(PATH_ANALYSIS_OUTPUT, "Ordered_Table3_Markers_NK2vsAllnonNK.xlsx")
write.xlsx(Table3_NK2vsAllnonNK, output_file, rowNames = FALSE)


#Table 4
#Find the markers
Table4_NK3vsAllnonNK = FindMarkers(reference, ident.1 =  Cell_Names_NK3 ,  ident.2 = Cell_Names_reference_not_NK, only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE )

# Order the data frame by avg_log2FC in descending order
# Include rownames as a column and order by avg_log2FC
Table4_NK3vsAllnonNK <- Table4_NK3vsAllnonNK %>%
  tibble::rownames_to_column(var = "Gene") %>%
  arrange(desc(avg_log2FC))

#Have a quick look at it
Table4_NK3vsAllnonNK

# Save the ordered data frame as an Excel file
output_file <- file.path(PATH_ANALYSIS_OUTPUT, "Ordered_Table4_Markers_NK3vsAllnonNK.xlsx")
write.xlsx(Table4_NK3vsAllnonNK, output_file, rowNames = FALSE)


