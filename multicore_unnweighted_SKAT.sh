#!/bin/bash

# Check if genotype prefix, base directory, and genotype folder are provided
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
  echo "Usage: $0 <genotype_prefix> <base_pipeline_directory> <genotype_folder>"
  exit 1
fi

# Set genotype prefix, base directory, and genotype folder from command-line arguments
genotype_prefix=$1
base_dir=$2
genotype_folder=$3

# Create output directory if it doesn't exist
output_dir="output_association"
mkdir -p $output_dir

# Get the number of available cores
num_cores=$(nproc)

# Print out the number of cores
echo "Using $num_cores cores for parallel processing."

# Run Rscript in parallel for chromosomes 1 to 22 using the genotype folder as input
seq 1 22 | parallel -j $num_cores --verbose Rscript ${base_dir}/scripts/skat_unweighted.R $genotype_prefix ${base_dir}/helper/gene_regions.csv $genotype_folder {} $output_dir/

# Confirm that all jobs are complete
echo "All parallel jobs are complete. Proceeding to merge results."

# Merge all results into a single file
echo -e "Gene_name\tGene_chromosome\tRegion_start\tRegion_end\tQ_test\tpvalue" > "${genotype_prefix}_assoc_test.tsv"
for file in "${output_dir}/${genotype_prefix}"*; do
  if [ "$file" != "${genotype_prefix}_assoc_test.tsv" ]; then
    tail -n +2 "$file" >> "${genotype_prefix}_assoc_test.tsv"
  fi
done

echo "Merging complete."
