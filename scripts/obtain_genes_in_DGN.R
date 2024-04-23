# Check the results in the gene annotation databases

# Obtain the gene names from SKAT results
get_gene_from_SKAT <- function(skat_results_path, pvalue_threshold) {
  gene_result <- read.csv(skat_results_path, sep = "\t")
  gene_result <- subset(gene_result, pvalue < pvalue_threshold)
  genes <- gene_result$Gene_name
  return(genes)
}

# Get intersection in Disgenet
get_intersect_dgn <- function(genes, disgenet_path) {
  # Load the dataset for Disease gene associations from Disgenet
  g2d_disease <- read.csv(disgenet_path, sep = "\t")
  g2d_disease_genes <- g2d_disease$Gene
  common_genes_dgn <- intersect(genes, g2d_disease_genes)
  return(common_genes_dgn)
}

# Function to get genes from the input list that are not present in Disgenet database
get_diff_dgn <- function(genes, disgenet_path) {
  # Load the dataset for Disease gene associations from Disgenet
  g2d_disease <- read.csv(disgenet_path, sep = "\t")
  g2d_disease_genes <- g2d_disease$Gene

  # Find genes from the input list that are not present in Disgenet database
  diff_genes <- setdiff(genes, g2d_disease_genes)
  return(diff_genes)
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 3) {
  cat("Usage: Rscript script.R skat_results_path disgenet_path pvalue_threshold\n")
  quit(status = 1)
}

# Parse command-line arguments
skat_results_path <- args[1]
disgenet_path <- args[2]
pvalue_threshold <- as.numeric(args[3])

# Obtain significant genes from SKAT results
genes <- get_gene_from_SKAT(skat_results_path, pvalue_threshold)

# Get genes from SKAT results that are not present in Disgenet database
diff_genes <- get_diff_dgn(genes, disgenet_path)

# Get common genes in Disgenet
common_genes_dgn <- get_intersect_dgn(genes, disgenet_path)

# Calculate the ratio of significant genes in Disgenet over total significant genes
validation_ratio_dgn <- length(common_genes_dgn) / length(genes)

# Calculate the ratio of significant genes not in Disgenet over total significant genes
validation_ratio_diff <- length(diff_genes) / length(genes)

# Prepare output text
output_text <- paste0("P-value threshold: ", pvalue_threshold, "\n",
                      "Number of genes in SKAT results: ", length(genes), "\n",
                      "Number of common genes in Disgenet: ", length(common_genes_dgn), "\n",
                      "Number of significant genes not in Disgenet: ", length(diff_genes), "\n",
                      "Validation ratio as Number of genes in DGN/Number of significant genes: ", validation_ratio_dgn, "\n",
                      "Validation ratio as Number of significant genes not in DGN/Number of significant genes: ", validation_ratio_diff, "\n\n")

# Write output to a text file
file_name <- basename(skat_results_path)  # Extract file name without extension
pvalue_threshold_str <- gsub("\\.", "_", as.character(pvalue_threshold))  # Replace dot with underscore for valid filename
output_file <- paste0("output_plots/output_results_", pvalue_threshold_str, "_", gsub("\\..*$", "", file_name), ".txt")  # Construct output file path
# Open the file for writing
file_conn <- file(output_file, "w")
# Write the output text to the file
writeLines(output_text, file_conn)
# Close the file connection
close(file_conn)
cat("Output written to:", output_file, "\n")
