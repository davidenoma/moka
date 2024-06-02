# Load the necessary libraries
library(ggplot2)
library(qqman)


# create_manhattan_plot <- function(data, threshold) {
#   ggplot(data, aes(x = as.factor(Gene_chromosome), y = neg_log_pvalue, label = ifelse(neg_log_pvalue > threshold, Gene_name, ""))) +
#     geom_point(aes(color = as.factor(Gene_chromosome)), size = 2.5, alpha = 0.8, position = position_dodge(width = 0.5)) +
#     scale_color_manual(values = rep(c("grey", "skyblue","pink"), 22 )) +
#     geom_text(vjust = -0.5, hjust = 1, size = 4) +
#     geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
#     labs(
#       x = "Chromosome",
#       y = "-log10(pvalue)",
#       title = "Manhattan Plot",
#       subtitle = "P-values for Genes",
#       caption = paste("p-value_threshold: ", 10^-(threshold))
#     ) +
#
#     theme_classic() +
#     theme(
#       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#       panel.grid.minor = element_blank()
#     ) +
#     theme(
#       legend.position="none",
#       panel.border = element_blank(),
#       panel.grid.major.x = element_blank(),
#       panel.grid.minor.x = element_blank()
#     )+
#   scale_x_discrete(labels = 1:22)
# }


# Ensure the "pvalue" column is numeric
# gene_data$pvalue <- as.numeric(gene_data$pvalue)
# Calculate the negative logarithm (base 10) of the p-values
# gene_data$neg_log_pvalue <- -log10(gene_data$pvalue)

# Example of usage
# Set the p-value threshold for labeling genes
# pvalue_threshold <- 0.05/1735
# pvalue_threshold_EU <- 0.05/2414

# create_manhattan_plot(gene_data, pvalue_threshold_EU)

# genes <- read.csv("/Users/davidenoma/Desktop/genes_test.txt", header = FALSE, sep="\t", col.names = c("SNP", "CHR", "BP", "P"))
# Read the gene data from the tab-delimited file

# genes <- read.csv("/Users/davidenoma/Desktop/genes_test.txt",header = FALSE,sep="\t",col.names = )


# Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if there are arguments passed
if (length(args) == 0) {
  cat("Usage: Rscript script_name.R <file_path> <pvalue_threshold>\n")
  quit(status = 1)
}

# Assign values to arguments
file_path <- args[1]
pvalue_threshold <- as.numeric(args[2])
cat("P-value threshold:",pvalue_threshold)
gene_data <- read.table(
  file_path,
  sep = "\t",  # Specify tab as the separator
  header = TRUE,  # Use the first row as column names
  col.names = c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
)
selected_gene_data <- gene_data[, c("Gene_name", "Gene_chromosome", "Region_start", "pvalue")]
colnames(selected_gene_data) <- c("SNP", "CHR", "BP", "P")

# Save the Manhattan plot
file_name <- basename(file_path)  # Extract file name without extension
png_file_path <- paste0("output_plots/manhattan_", gsub("\\..*$", "", file_name), ".png")  # Construct PNG file path
# Create the directory if it does not exist
dir.create(dirname(png_file_path), showWarnings = FALSE, recursive = TRUE)
png(filename = png_file_path, width = 2940, height = 1782, units = "px", pointsize = 22,res = 150)
# png(filename = png_file_path, width = 2000, height = 1600, units = "px", pointsize = 16)
# Call the function to create Manhattan plot
manhattan(selected_gene_data, chr = "CHR", bp = "BP", snp = "SNP",
          p = "P", col = c("grey", "skyblue", "pink"),
          annotatePval = pvalue_threshold, annotateTop = FALSE,
          genomewideline = -log10(pvalue_threshold), suggestiveline = FALSE,
          logp = TRUE)

# Close the PNG device
dev.off()
