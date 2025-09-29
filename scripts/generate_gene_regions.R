# Install and load rtracklayer if not available
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("rtracklayer", quietly = TRUE))
    BiocManager::install("rtracklayer")

library(rtracklayer)

# Function to load GFF and extract gene features
get_gene_features <- function(gff_path) {
  gff_url <- "https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz"
  if (!file.exists(gff_path)) {
    dir.create(dirname(gff_path), showWarnings = FALSE)
    download.file(gff_url, destfile = gff_path, mode = "wb")
  }
  gff <- readGFF(gff_path)
  gene_features <- subset(gff, type == "gene" & seqid %in% paste0("chr", 1:22))
  gene_features <- gene_features[complete.cases(gene_features$Name), ]
  return(gene_features)
}

# Function to create gene regions with flanking
create_gene_regions <- function(gene_features, flank_size) {
  gene_name <- gene_features$Name
  gene_start <- gene_features$start
  gene_end <- gene_features$end
  gene_chromosome <- sub("^chr", "", gene_features$seqid)
  region_start <- as.integer(pmax(1, gene_start - flank_size))
  region_end <- as.integer(gene_end + flank_size)
  gene_regions <- data.frame(
    GeneName = gene_name,
    Start = region_start,
    Stop = region_end,
    Chromosome = gene_chromosome,
    stringsAsFactors = FALSE
  )
  return(gene_regions)
}

# Main script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript generate_gene_regions.R <flank_size>")
}
flank_size <- as.integer(args[1])
cat("Using flank size:", flank_size, "bp\n")

gff_path <- "helper/Homo_sapiens.GRCh38.115.gff3"
gene_features <- get_gene_features(gff_path)
gene_regions <- create_gene_regions(gene_features, flank_size)

output_file <- paste0("helper/gene_regions_", as.character(flank_size), "bp.csv")
write.csv(gene_regions, output_file, row.names = FALSE)
