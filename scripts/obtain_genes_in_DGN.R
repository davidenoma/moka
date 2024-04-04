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

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 3) {
  cat("Usage: Rscript script.R skat_results_path disgenet_path pvalue_threshold\n")
  quit(status = 1)
}
skat_results_path <- args[1]
disgenet_path <- args[2]
pvalue_threshold <- as.numeric(args[3])
genes <- get_gene_from_SKAT(skat_results_path, pvalue_threshold)
common_genes_dgn <- get_intersect_dgn(genes, disgenet_path)
cat("P-value threshold: ",pvalue_threshold, "\n")
cat("Number of common genes in Disgenet: " , length(common_genes_dgn), "\n")
cat("Validation ratio as Number of genes in DGN/Number of significant genes: ",length(common_genes_dgn)/length(genes), "\n")
cat("\n")


# Prepare output text
output_text <- paste("P-value threshold: ", pvalue_threshold, "\n",
                     "Number of common genes in Disgenet: ", length(common_genes_dgn), "\n",
                     "Validation ratio as Number of genes in DGN/Number of significant genes: ",
                     length(common_genes_dgn)/length(genes), "\n\n")

# Write output to a text file
file_name <- basename(skat_results_path)  # Extract file name without extension
output_file <- paste0("output_plots/output_results_", gsub("\\..*$", "", file_name), ".txt")  # Construct output file path
writeLines(output_text, output_file)

cat("Output written to:", output_file, "\n")