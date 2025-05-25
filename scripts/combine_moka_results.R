# Define the function combine_skat_results
combine_skat_results <- function(genotype_prefix, weights_type, result_folder) {
  # Define the output file name with weights_type included
  output_file <- paste0(genotype_prefix, "_", weights_type, "_combined_association.tsv")

  # Execute shell commands to combine SKAT files using pattern matching
  system(paste('echo "Gene_name\tGene_chromosome\tRegion_start\tRegion_end\tQ_test\tpvalue" >', file.path(result_folder, output_file)))
  system(paste('for file in', file.path(result_folder, paste0("*", genotype_prefix, "*", weights_type, "*")), '; do if [ "$file" !=', file.path(result_folder, output_file), ']; then tail -n +2 "$file" >>', file.path(result_folder, output_file), '; fi; done'))

  # Optional: Print message indicating completion
  cat("SKAT association results combined successfully.\n")
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 3) {
  cat("Usage: Rscript combine_results.R <genotype_prefix> <weights_type> <result_folder>\n")
  quit(status = 1)
}

# Assign command-line arguments to variables
genotype_prefix <- args[1]
weights_type <- args[2]
result_folder <- args[3]

# Call the function combine_skat_results
combine_skat_results(genotype_prefix, weights_type, result_folder)
