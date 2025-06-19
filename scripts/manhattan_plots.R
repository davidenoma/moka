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
  # col.names = c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
)

# Calculate the p-value threshold using Bonferroni correction
num_rows <- nrow(gene_data)
pvalue_threshold <- 0.05 / num_rows
cat("P-value threshold:", pvalue_threshold, "\n")

# Prepare data for the Manhattan plot
selected_gene_data <- gene_data[, c("Gene_name", "Gene_chromosome", "Region_start", "pvalue")]
colnames(selected_gene_data) <- c("SNP", "CHR", "BP", "P")

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

top_snps <- head(selected_gene_data[order(selected_gene_data$P), ], 10)
top_idx <- which.max(top_snps$P)
# Extract the p-value for that SNP
top_pvalue <- top_snps$P[top_idx]

# Create the Manhattan plot using CairoPNG
CairoPNG(filename = png_file_path, width = 2940, height = 1782, units = "px", pointsize = 20, res = 250)
# quartz()
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
  cex.main= 3.0,
  cex = 1.0
)

 # with(subset(selected_gene_data, P < top_pvalue), textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 1.5 ))
# ── NEW lines: draw larger gene labels ───────────────────────────────

## ─────────────────────────────────

# Wait for the user to press enter so the plot stays open
# readline(prompt = "Press [enter] to close the plot window.")
# Annotate the top 10 SNPs with smallest p-values
# with(top_snps, text(BP, -log10(P), labels = SNP, pos = 3, cex = 0.8))
# dev.print()
dev.off()
cat("Manhattan plot saved to:", png_file_path, "\n")
