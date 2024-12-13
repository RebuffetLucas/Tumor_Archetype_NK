#Read the Seurat Object
Seurat_NK <- readRDS(file.path(PATH_EXPERIMENT_REFERENCE , SAMPLE_NK_CITEseq))


#Remove and Prolif recluster 
Seurat_NK = subset(Seurat_NK, idents= "NK Proliferating", invert= TRUE) #Remove the prolif NK cells
Seurat_NK = subset(Seurat_NK, subset= time== "0" ) # Keep only before clinical trial injection

#Normalize RNA and ADT
DefaultAssay(object = Seurat_NK) <- "SCT"
Seurat_NK = SCTransform(Seurat_NK, assay="SCT")

DefaultAssay(object = Seurat_NK) <- "ADT"
Seurat_NK = NormalizeData(Seurat_NK, normalization.method = 'CLR', margin = 2)

#Rerun wnn
Seurat_NK2 <- FindMultiModalNeighbors( Seurat_NK, reduction.list = list("pca", "apca"),   dims.list = list(1:50, 1:50), modality.weight.name = "SCT.weight" )
Seurat_NK2 <- RunUMAP(Seurat_NK2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Seurat_NK2 <- FindClusters(Seurat_NK2, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)

DimPlot(Seurat_NK2, reduction = "wnn.umap")

Seurat_NK3 = subset(Seurat_NK2, ident= "3", invert = TRUE ) #Remove µcluster of remaining contaminations

levels(Seurat_NK3$seurat_clusters) = c("NK1", "NK3", "NK2", "NK_Prolif")
Seurat_NK3 = SetIdent(Seurat_NK3, value = "seurat_clusters")


saveRDS(Seurat_NK3, file.path(PATH_ANALYSIS_OUTPUT, "Seurat_NKonly_MetaNK_Preprocessed.rds")) #Save the pre-processed Object
#Table 5
Table5_All_Markers_NK = FindMarkers(Seurat_NK3, ident.1 = "NK1", only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE )

Table5_All_Markers_NK <- Table5_All_Markers_NK %>%
  tibble::rownames_to_column(var = "Gene") %>%
  arrange(desc(avg_log2FC))

output_file <- file.path(PATH_ANALYSIS_OUTPUT, "Ordered_Table5_NK1vsNK2_NK3.xlsx")
write.xlsx(Table5_All_Markers_NK, output_file, rowNames = FALSE)


#Table 6
Table6_All_Markers_NK = FindMarkers(Seurat_NK3, ident.1 = "NK2", only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE )

Table6_All_Markers_NK <- Table6_All_Markers_NK %>%
  tibble::rownames_to_column(var = "Gene") %>%
  arrange(desc(avg_log2FC))

output_file <- file.path(PATH_ANALYSIS_OUTPUT, "Ordered_Table6_NK2vsNK1_NK3.xlsx")
write.xlsx(Table6_All_Markers_NK, output_file, rowNames = FALSE)

#Table 5
Table7_All_Markers_NK = FindMarkers(Seurat_NK3, ident.1 = "NK3", only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE )

Table7_All_Markers_NK <- Table7_All_Markers_NK %>%
  tibble::rownames_to_column(var = "Gene") %>%
  arrange(desc(avg_log2FC))

output_file <- file.path(PATH_ANALYSIS_OUTPUT, "Ordered_Table7_NK3vsNK1_NK2.xlsx")
write.xlsx(Table7_All_Markers_NK, output_file, rowNames = FALSE)

#Saving the cell identities
for (clust_subs in c("NK1" , "NK2" , "NK3")){
  SUBSET_OBJECT = subset(Seurat_NK3 , ident = clust_subs)
  output_file <- file.path(PATH_ANALYSIS_OUTPUT, "Cell_Names", paste0(clust_subs, "_Cell_Names.rds"))
  saveRDS(Cells(SUBSET_OBJECT), output_file)
}
