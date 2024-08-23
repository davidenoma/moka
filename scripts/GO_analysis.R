library(gprofiler2)
library(ggplot2)  

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

# Specify the PNG file for the GO pathway plot
GO_png_file_path <- paste0("output_plots/GO_", gsub("\\..*$", "", basename(file_path)), ".png")

# Generate GO pathway plot
png(filename = GO_png_file_path, width = 1700, height = 1800, units = "px", res = 150)

# Perform GO enrichment analysis
gost_res <- gost(
  query = gene_and_p_values$GENE,
  organism = "hsapiens",
  sources = c("GO:MF", "GO:BP", "GO:CC"),
  highlight = TRUE,
  evcodes = TRUE
)

# Generate GO enrichment plot
p <- gostplot(
  gostres = gost_res,
  pal = c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618"),
  interactive = FALSE,
  capped = TRUE
)

# Increase font size by modifying the plot object
p <- p + theme(
  text = element_text(size = 24),  # Adjust font size for all text
  axis.title = element_text(size = 20),  # Font size for axis titles
  axis.text = element_text(size = 20)  # Font size for axis labels
)

# Sort gostres$result dataframe by p_value column from lowest to highest
sorted_gostres <- gost_res$result[order(gost_res$result$p_value), ]

# Get the lowest 10 rows
lowest_10 <- sorted_gostres[1:10, ]

# Publish GO enrichment plot with the lowest 10 terms highlighted
publish_gostplot(p, highlight_terms = lowest_10$term_id)

dev.off()
