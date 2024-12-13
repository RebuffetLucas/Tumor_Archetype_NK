#Draw upset plot of some tables
## Load the dataframes
# Define the folder containing the Excel files
folder_path <- PATH_ANALYSIS_OUTPUT  # Replace with your actual folder path

# Get all Excel file names in the folder
excel_files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Create a list to store the data frames
excel_list <- list()

# Loop through each file and load it into the list
for (file in excel_files) {
  # Extract the base file name (without the extension)
  file_name <- basename(file)
  
  # Extract the table number using regex (assumes "Tablex" pattern in file name)
  table_name <- sub(".*(Table\\d+).*", "\\1", file_name)
  
  # Read the Excel file and store it in the list with the extracted name
  excel_list[[table_name]] <- read_xlsx(file)
}

# Print the names of the list to verify
print(names(excel_list))


# Select the top genes from each table
top_genes <- lapply(list_to_apply, function(table_name) {
  excel_list[[table_name]] %>% 
    slice_head(n = NUMBER_MARKERS_UPSET) %>% 
    pull(Gene)
})
names(top_genes) <- list_to_apply

# Create a binary matrix indicating the presence of each gene in the lists
all_genes <- unique(unlist(top_genes))
gene_matrix <- data.frame(
  Gene = all_genes,
  do.call(cbind, lapply(list_to_apply, function(table_name) {
    as.integer(all_genes %in% top_genes[[table_name]])
  })))
colnames(gene_matrix)[-1] <- list_to_apply

# Convert to a data frame suitable for UpSetR
upset_data <- gene_matrix
row.names(upset_data) <- upset_data$Gene
upset_data <- upset_data[, -1] # Remove the "Gene" column

# Draw the UpSet plot
print(upset(upset_data, 
            nsets = length(list_to_apply), # Number of sets dynamically set
            nintersects = NA, # Include all intersections
            order.by = "freq", 
            main.bar.color = "black", 
            sets.bar.color = "darkblue", 
            text.scale = c(1.5, 1.5, 1, 1, 1.5, 1),
            mainbar.y.label= paste("Top", NUMBER_MARKERS_UPSET ,  "Markers")))

