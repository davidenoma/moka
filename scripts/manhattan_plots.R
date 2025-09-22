# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install and load required packages
# packages <- c("ggplot2", "qqman")
# for (pkg in packages) {
#   if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#     install.packages(pkg)
#     library(pkg, character.only = TRUE)
#   }
# }
library(ggplot2)
library(qqman)

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
  header = TRUE
  # col.names = c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
)

# Calculate the p-value threshold using Bonferroni correction
num_rows <- nrow(gene_data)
pvalue_threshold <- 0.05 / num_rows
cat("P-value threshold:", pvalue_threshold, "\n")

# Prepare data for the Manhattan plot
selected_gene_data <- gene_data[, c("Gene_name", "Gene_chromosome", "Region_start", "pvalue")]
colnames(selected_gene_data) <- c("SNP", "CHR", "BP", "P")

# Ensure CHR is numeric if possible (qqman expects numeric chromosomes)
# If your chromosomes are like "X", "Y", "MT", convert or drop accordingly.
suppressWarnings({
  selected_gene_data$CHR <- as.numeric(selected_gene_data$CHR)
})

# Ensure p-values are finite and non-zero
selected_gene_data <- selected_gene_data[!is.na(selected_gene_data$P) & selected_gene_data$P > 0, ]

# Check if there are valid p-values left
if (nrow(selected_gene_data) == 0) {
  cat("Error: No valid p-values found after filtering.\n")
  quit(status = 1)
}

# Create output directory if it doesn't exist
if (!dir.exists("output_plots")) dir.create("output_plots")
file_name <- basename(file_path)
png_file_path <- paste0("output_plots/manhattan_", sub("\\.[^.]*$", "", file_name), ".png")

# Get top SNPs by smallest p-values (top 10)
top_snps <- head(selected_gene_data[order(selected_gene_data$P), ], 10)
# Threshold for annotation = largest p among the selected top SNPs
top_pvalue <- max(top_snps$P, na.rm = TRUE)

# Open a standard PNG device (no Cairo)
png(filename = png_file_path, width = 2940, height = 1782, units = "px", pointsize = 20, res = 250)

manhattan(
  selected_gene_data,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "P",
  col = c("grey", "skyblue", "pink"),
  annotatePval = top_pvalue,
  highlight = top_snps$SNP,
  annotateTop = FALSE,
  genomewideline = -log10(pvalue_threshold),
  suggestiveline = FALSE,
  logp = TRUE,
  cex.main = 3.0,
  cex = 1.0
)

dev.off()
cat("Manhattan plot saved to:", png_file_path, "\n")
