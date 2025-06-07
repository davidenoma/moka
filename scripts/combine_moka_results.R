# File: `scripts/combine_moka_results.R`
combine_skat_results <- function(genotype_prefix, weights_type, result_folder) {
  if (weights_type == "" || is.na(weights_type)) {
    output_file <- paste0(genotype_prefix, "_combined_association.tsv")
    file_pattern <- file.path(result_folder, paste0("*", genotype_prefix, "*"))
  } else {
    output_file <- paste0(genotype_prefix, "_", weights_type, "_combined_association.tsv")
    file_pattern <- file.path(result_folder, paste0("*", genotype_prefix, "*", weights_type, "*"))
  }

  output_path <- file.path(result_folder, output_file)

  # Write header
  system(paste('echo "Gene_name\tGene_chromosome\tRegion_start\tRegion_end\tQ_test\tpvalue" >', output_path))

  # Loop through matching files and combine them
  system(paste('for file in', file_pattern, '; do if [ "$file" !=', output_path, ']; then tail -n +2 "$file" >>', output_path, '; fi; done'))

  cat("SKAT association results combined successfully.\n")
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Expecting at least two arguments; if weights type is missing, set to empty string
if (length(args) < 2) {
  cat("Usage: Rscript combine_results.R <genotype_prefix> <result_folder> [weights_type]\n")
  quit(status = 1)
}

genotype_prefix <- args[1]
if (length(args) >= 3) {
  weights_type <- args[2]
  result_folder <- args[3]
} else {
  # If weights type is not provided, use empty string and adjust parameters accordingly
  weights_type <- ""
  result_folder <- args[2]
}

combine_skat_results(genotype_prefix, weights_type, result_folder)