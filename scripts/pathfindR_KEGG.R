library(pathfindR)

# Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if there are arguments passed
if (length(args) < 2) {
  cat("Usage: Rscript script_name.R <file_path> <pvalue_threshold>\n")
  quit(status = 1)
}

# Assign values to arguments
file_path <- args[1]
pvalue_threshold <- as.numeric(args[2])
cat("P-value threshold:", pvalue_threshold)

# Read data from the specified file
gene_data <- read.table(
  file_path,
  sep = "\t",  # Specify tab as the separator
  header = TRUE,  # Use the first row as column names
  col.names = c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
)

# Filter data based on the p-value threshold
gene_and_p_values <- subset(gene_data, pvalue < pvalue_threshold)
gene_and_p_values <- gene_and_p_values[, c("Gene_name", "pvalue")]
colnames(gene_and_p_values) <- c("GENE", "p_value")


# Specify the PNG file for the KEGG pathway plot
kegg_png_file_path <- paste0("output_plots/kegg_pathway_", gsub("\\..*$", "", basename(file_path)), ".png")
# Generate KEGG pathway plot
png(filename = kegg_png_file_path, width = 1500, height = 1000,units = "px",pointsize = 22,res = 150)
# Run pathfindR with the filtered data
enrichment_res <- run_pathfindR(gene_and_p_values)
dev.off()
