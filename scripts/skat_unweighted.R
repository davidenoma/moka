# Load required R libraries
library(SKAT)

# Define the function to perform the SKAT test
perform_skat_test <- function(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, result_file, is_binary = TRUE) {
  tryCatch({
    genotype_prefix <- paste0(genotype_path, genotype_prefix)
    genotype_bed <- paste0(genotype_prefix, ".bed")

    # Read the SNP list file
    snp_list <- read.table(paste0(genotype_path, "snp_list_", chr, ".txt"), header = FALSE, col.names = "SNP")

    system(paste(
      "plink --bfile",
      genotype_prefix,
      "--extract",
      paste0(genotype_path, "snp_list_", chr, ".txt"),
      "--make-bed",
      "--silent",
      "--allow-no-sex",
      "--out",
      paste0(genotype_prefix, "_SKAT")
    ))

    # Prepend "SET" to the SNP IDs
    snp_list <- paste0("SET\t", snp_list$SNP)

    # Write the modified SNP list to the output file
    write.table(snp_list, paste0(genotype_prefix, ".setid"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    Generate_SSD_SetID(genotype_bed, paste0(genotype_prefix, ".bim"), paste0(genotype_prefix, ".fam"), paste0(genotype_prefix, ".setid"), paste0(genotype_prefix, ".ssd"), paste0(genotype_prefix, ".info"))

    # Perform SKAT test based on phenotype type (binary or continuous)
    genotype_fam <- paste0(genotype_prefix, ".fam")
    FAM <- Read_Plink_FAM(genotype_fam, Is.binary = is_binary)
    y <- FAM$Phenotype  # This is the outcome variable

    # Set the null model based on outcome type
    if (is_binary) {
      obj <- SKAT_Null_Model(y ~ 1, out_type = "D")  # Binary outcome
    } else {
      obj <- SKAT_Null_Model(y ~ 1, out_type = "C")  # Continuous outcome
    }

    genotype_ssd <- paste0(genotype_prefix, ".ssd")
    SSD.INFO <- Open_SSD(genotype_ssd, paste0(genotype_prefix, ".info"))
    id <- 1
    SetID <- SSD.INFO$SetInfo$SetID[id]

    # Use 'gene_snps' as the subset of SNPs for SKAT test
    skat_test <- SKAT.SSD.OneSet(SSD.INFO, SetID, obj, kernel = "linear.weighted", obj.SNPWeight = NULL)

    # Append SKAT test results to result file
    ss2 <- c(gene_name, gene_chromosome, region_start, region_end, toString(skat_test$Q), toString(skat_test$p.value))
    write(ss2, file = result_file, append = TRUE, ncol = 6, sep = '\t')

    # Clean up temporary files
    file.remove(paste0(genotype_prefix, "_SKAT.bed"))
    file.remove(paste0(genotype_prefix, "_SKAT.bim"))
    file.remove(paste0(genotype_prefix, "_SKAT.fam"))
    file.remove(paste0(genotype_prefix, "_SKAT.log"))
    file.remove(paste0(genotype_prefix, ".setid"))
    file.remove(paste0(genotype_prefix, ".ssd"))
    file.remove(paste0(genotype_prefix, ".info"))

  }, error = function(e) {
    cat("Error in SKAT analysis for gene:", gene_name, "\n")
    cat("Error message:", e$message, "\n")
  })
}

# Perform SKAT for a specific chromosome
snvs_and_skat_chr <- function(genotype_prefix, gene_regions_file, genotype_path, chr, result_folder, is_binary = TRUE) {

  gene_regions <- read.csv(gene_regions_file, header = TRUE)
  system(paste0("mkdir -p ", result_folder))
  result_file <- file.path(result_folder, paste0(genotype_prefix, "_result_chr_", chr, ".txt"))

  ss <- c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
  write(ss, file = result_file, ncol = 6, sep = '\t')

  gene_count <- 0
  gene_regions <- subset(gene_regions, Chromosome == chr)

  snp_data <- read.csv(paste0(genotype_path, genotype_prefix, ".bim"), sep = "\t", header = FALSE, col.names = c("chr", "snp_id", "pos_cm", "loc", "alt", "ref"))
  snp_data <- subset(snp_data, chr == chr)

  for (i in 1:nrow(gene_regions)) {
    gene_name <- gene_regions$GeneName[i]
    region_start <- gene_regions$Start[i]
    region_end <- gene_regions$Stop[i]
    gene_chromosome <- gene_regions$Chromosome[i]

    gene_snps <- subset(snp_data, loc >= region_start & loc <= region_end)

    if (nrow(gene_snps) > 0) {
      write.table(gene_snps$snp_id, paste0(genotype_path, "snp_list_", chr, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      print(paste("Gene number: ", gene_count))

      # Call the SKAT test function with selected SNPs
      perform_skat_test(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, result_file, is_binary)
      gene_count <- gene_count + 1
      file.remove(paste0(genotype_path, "snp_list_", chr, ".txt"))
    } else {
      print("Not found for gene region")
    }
  }
  print(paste("Done with SKAT analysis for chromosome", chr))
}

# Helper functions to handle SKAT files for each chromosome
prepare_SKAT_files_per_chr <- function(genotype_path, genotype_prefix)  {
  old_file_name_bed <- paste0(genotype_path, genotype_prefix, ".bed")
  new_file_name_bed <- paste0(genotype_path, genotype_prefix, "_", chr, ".bed")
  system(paste("cp ", old_file_name_bed, new_file_name_bed))

  old_file_name_bim <- paste0(genotype_path, genotype_prefix, ".bim")
  new_file_name_bim <- paste0(genotype_path, genotype_prefix, "_", chr, ".bim")
  system(paste("cp ", old_file_name_bim, new_file_name_bim))

  old_file_name_fam <- paste0(genotype_path, genotype_prefix, ".fam")
  new_file_name_fam <- paste0(genotype_path, genotype_prefix, "_", chr, ".fam")
  system(paste("cp ", old_file_name_fam, new_file_name_fam))
}

remove_SKAT_files_per_chr <- function(genotype_path, genotype_prefix)  {
  old_file_name_bed <- paste0(genotype_path, genotype_prefix, ".bed")
  system(paste("rm ", old_file_name_bed))

  old_file_name_bim <- paste0(genotype_path, genotype_prefix, ".bim")
  system(paste("rm ", old_file_name_bim))

  old_file_name_fam <- paste0(genotype_path, genotype_prefix, ".fam")
  system(paste("rm ", old_file_name_fam))

  old_file_name_log <- paste0(genotype_path, genotype_prefix, ".log")
  system(paste("rm ", old_file_name_log))

  old_file_name_info <- paste0(genotype_path, genotype_prefix, ".info")
  system(paste("rm ", old_file_name_info))
}

# Combine SKAT results from multiple chromosomes
combine_skat_results <- function(genotype_prefix, result_folder) {
  output_file <- file.path(result_folder, paste0(genotype_prefix, "_combined_association.tsv"))

  if (!dir.exists(result_folder)) {
    stop("Result folder does not exist.")
  }

  write.table(matrix(c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue"), nrow=1),
              file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  result_files <- list.files(result_folder, pattern = paste0(genotype_prefix, "_result_chr_"), full.names = TRUE)
  for (file in result_files) {
    if (file != output_file) {
      write.table(read.table(file, header = TRUE, sep = "\t"),
                  file = output_file, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = "\t")
    }
  }

  cat("SKAT association results combined successfully into:", output_file, "\n")
}

# Main execution starts here
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript skat_analysis.R <genotype_prefix> <gene_regions_file> <genotype_path> <chromosome> <output_folder> [is_binary]\n")
  quit(status = 1)
}

# Extract command-line arguments
genotype_prefix <- args[1]
gene_regions_file <- args[2]
genotype_path <- args[3]
chr <- args[4]
output_folder <- args[5]
is_binary <- ifelse(length(args) > 5, as.logical(args[6]), TRUE)

# Perform SKAT analysis for the given chromosome
prepare_SKAT_files_per_chr(genotype_path, genotype_prefix)
print("Done creating SKAT genotype")
genotype_prefix <- paste0(genotype_prefix, "_", chr)
snvs_and_skat_chr(genotype_prefix, gene_regions_file, genotype_path, chr, output_folder, is_binary)

# Combine results from all chromosomes
combine_skat_results(genotype_prefix, output_folder)

# Clean up SKAT files after analysis
remove_SKAT_files_per_chr(genotype_path, genotype_prefix)

cat("Done with SKAT analysis\n")
