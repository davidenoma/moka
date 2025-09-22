library(pathfindR)
library(ggplot2)  # Required for customizing the plot
library(Cairo)
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
cat("P-value threshold:", pvalue_threshold, "\n")

# Read data from the specified file
gene_data <- read.table(
  file_path,
  sep = "\t",  # Specify tab as the separator
  header = TRUE,  # Use the first row as column names
  col.names = c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
)

# Adjust p-values using the Benjamini-Hochberg method
gene_data$adjusted_pvalue <- p.adjust(gene_data$pvalue, method = "fdr",n=length(gene_data$pvalue))

# Filter data based on the p-value threshold
gene_and_p_values <- subset(gene_data, adjusted_pvalue < pvalue_threshold)
gene_and_p_values <- gene_and_p_values[, c("Gene_name", "adjusted_pvalue")]
colnames(gene_and_p_values) <- c("GENE", "p_value")


if (nrow(gene_and_p_values) == 0) {
  cat("No genes passed the p-value threshold after adjustment. Exiting.\n")
  quit(status = 0)
}


# Run pathfindR with the filtered data
enrichment_res <- run_pathfindR(gene_and_p_values)

# Specify the PNG file for the KEGG pathway plot
kegg_png_file_path <- paste0("output_plots/kegg_pathway_", gsub("\\..*$", "", basename(file_path)), ".png")

# Open the PNG device
CairoPNG(filename = kegg_png_file_path, width = 1500, height = 1000, units = "px", pointsize = 22, res = 150)


# Create the enrichment chart and customize font size
enrichment_plot <- enrichment_chart(enrichment_res, top_terms = 10) +
  theme(
    text = element_text(size = 15),          # Increase font size for all text
    axis.title = element_text(size = 20),    # Font size for axis titles
    axis.text.y = element_text(size = 15),  # Increase font size of the pathway names
    axis.text.x = element_text(size = 20),  # Adjust x-axis font size if needed
    plot.title = element_text(size = 24, face = "bold"),  # Font size for the plot title
    legend.text = element_text(size = 16),   # Font size for legend text
    legend.title = element_text(size = 18)   # Font size for legend title
  )

# Render the plot
print(enrichment_plot)

# Close the PNG device to save the plot
dev.off()
