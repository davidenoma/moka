# Load necessary libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman")
if (!requireNamespace("Cairo", quietly = TRUE)) install.packages("Cairo")
library(ggplot2)
library(qqman)
library(Cairo)

# Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if there are arguments passed
if (length(args) == 0) {
  cat("Usage: Rscript script_name.R <file_path> \n")
  quit(status = 1)
}

# Assign values to arguments
file_path <- args[1]

# Read the gene data from the tab-delimited file
gene_data <- read.table(
  file_path,
  sep = "\t",  
  header = TRUE,  
  col.names = c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
)

# Calculate the p-value threshold using Bonferroni correction
num_rows <- nrow(gene_data)
pvalue_threshold <- 0.05 / num_rows
cat("P-value threshold:", pvalue_threshold, "\n")

# Prepare data for the Manhattan plot
selected_gene_data <- gene_data[, c("Gene_name", "Gene_chromosome", "Region_start", "pvalue")]
colnames(selected_gene_data) <- c("SNP", "CHR", "BP", "P")

# Remove the first data point (row)
selected_gene_data <- selected_gene_data[-1,]

# Create output directory if it doesn't exist
if (!dir.exists("output_plots")) dir.create("output_plots")
file_name <- basename(file_path)  
png_file_path <- paste0("output_plots/manhattan_", sub("\\.[^.]*$", "", file_name), ".png")

# Create the Manhattan plot using CairoPNG
CairoPNG(filename = png_file_path, width = 2940, height = 1782, units = "px", pointsize = 22, res = 150)
manhattan(
  selected_gene_data,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "P",
  col = c("grey", "skyblue", "pink"),
  annotatePval = pvalue_threshold,
  annotateTop = FALSE,
  genomewideline = -log10(pvalue_threshold),
  suggestiveline = -log10(10e-5),
  logp = TRUE,
  main = "Manhattan Plot",
  cex.axis = 0.7,  # Adjust axis labels to fit better
  cex.lab = 1,    # Adjust axis titles
  cex.main = 1.2,  # Adjust main title size
  xaxt = "n",      # Remove x-axis to manually control the space
  yaxt = "n"       # Remove y-axis to manually control the space
)

# Add the custom axis to remove space and adjust tick marks
axis(1, at = seq(1, max(selected_gene_data$CHR), by = 1), labels = as.character(seq(1, max(selected_gene_data$CHR), by = 1)), las = 2, cex.axis = 0.7)  # X-axis without gaps
axis(2, at = seq(0, 10, by = 2), las = 2, cex.axis = 0.7)  # Y-axis

# Draw a border around the plot manually
rect(0, 0, 21, 10, border = "black", lwd = 2)

dev.off()
cat("Manhattan plot saved to:", png_file_path, "\n")
