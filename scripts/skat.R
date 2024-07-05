# Load required R libraries
library(SKAT)
library(parallel)

# Define the function to perform the SKAT test
perform_skat_test <- function(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, result_folder, result_file) {
  tryCatch({
    # Rest of the function remains the same
    genotype_prefix <- paste0(genotype_path, genotype_prefix)
    genotype_bed <- paste0(genotype_prefix, ".bed")

    # Read the SNP list file
    snp_list <- read.table(paste0(genotype_path,"snp_list_", weights_type, chr,".txt"), header = FALSE, col.names = "SNP")

    system(paste(
      "plink --bfile",
      genotype_prefix,
      "--covar",
      paste0(genotype_prefix,".cov"),
      "--extract",
      paste0(genotype_path,"snp_list_", weights_type, chr,".txt"),
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
    # Perform SKAT test
    genotype_fam <- paste0(genotype_prefix, ".fam")
    genotype_cov <- paste0(genotype_prefix,".cov")
    # FAM <- Read_Plink_FAM(genotype_fam, Is.binary = TRUE)
    FAM_cov <- Read_Plink_FAM_Cov(genotype_fam, genotype_cov, Is.binary = FALSE,cov_header=TRUE)
    age <- FAM_cov$Age
    sex <- FAM_cov$Sex
    y <- FAM_cov$Phenotype
    obj <- SKAT_Null_Model(y ~ age + sex, out_type = "D")
    genotype_ssd <- paste0(genotype_prefix, ".ssd")
    SSD.INFO <- Open_SSD(genotype_ssd, paste0(genotype_prefix, ".info"))
    id <- 1
    Z <- Get_Genotypes_SSD(SSD.INFO, id)

    # Use 'gene_snps' as the subset of SNPs for SKAT test
    skat_test <- SKAT(Z, obj, kernel = "linear.weighted", weights = gene_snps$Weight)

    # Append SKAT test results to result file
    ss2 <- c(gene_name, gene_chromosome, region_start, region_end, toString(skat_test$Q), toString(skat_test$p.value))
    write(ss2, file = result_file, append = TRUE, ncol = 6, sep = '\t')


    # Clean up temporary files
    file.remove(paste0(genotype_prefix, "_SKAT.bed"), paste0(genotype_prefix, "_SKAT.bim"), paste0(genotype_prefix, "_SKAT.fam"), paste0(genotype_prefix, ".setid"), paste0(genotype_prefix, ".ssd"), paste0(genotype_prefix, ".info"))

    }, error = function(e) {
    cat("Error in SKAT analysis for gene:", gene_name, "\n")
    cat("Error message:", e$message, "\n")
  })
}
# Define a function to extract weights for SNVs and perform SKAT for a specific chromosome
extract_weights_for_snvs_and_skat_chr <- function(genotype_prefix, gene_regions_file, weights_file, genotype_path, weights_type, result_folder, chr) {
  # Read the weights vector file
  weights_vector <- read.csv(weights_file, sep = ",", header = TRUE)
  colnames(weights_vector) <- c("SNP", "Chr", "Pos", "Weight")

  # Read the gene regions files
  gene_regions <- read.csv(gene_regions_file, header = TRUE)

  weights_vector <- subset(weights_vector, Chr == chr)

  # Create results folder
  system("mkdir -p result_folder")
  result_file <- file.path(result_folder, paste0(genotype_prefix, "_", weights_type, "_result_chr_", chr, ".txt"))

  ss <- c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
  write(ss, file = result_file, ncol = 6, sep = '\t')

  gene_count <- 0
  gene_regions <- subset(gene_regions, Chromosome == chr)
  # Iterate through gene regions
  for (i in 1:nrow(gene_regions)) {
    gene_name <- gene_regions$GeneName[i]
    region_start <- gene_regions$Start[i]
    region_end <- gene_regions$Stop[i]
    gene_chromosome <- gene_regions$Chromosome[i]

    # Select SNPs within the gene region
    gene_snps <- subset(weights_vector, Pos >= region_start & Pos <= region_end)

    # Perform SKAT test if there are SNPs in the gene region
    if (nrow(gene_snps) > 0) {
      # Create a temporary SNP list file
      write.table(gene_snps$SNP, paste0(genotype_path,"snp_list_", weights_type, chr,".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      print(paste("Gene number: ", gene_count))
      # Call the SKAT test function with selected SNPs
      perform_skat_test(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, result_folder, result_file)
      gene_count <- gene_count + 1
      # Remove the temporary SNP list file
      file.remove(paste0(genotype_path,"snp_list_", weights_type, chr,".txt"))
    } else {
      print("Not found for gene region")
    }
  }
  print(paste("Done with SKAT analysis for chromosome", chr))
}

prepare_SKAT_files_per_chr <- function(genotype_path, genotype_prefix)  {
  old_file_name_bed <- paste0(genotype_path,genotype_prefix,".bed")
  new_file_name_bed <- paste0(genotype_path,genotype_prefix,"_",chr,".bed")
  system(paste("cp ",old_file_name_bed,new_file_name_bed))

  old_file_name_bim <- paste0(genotype_path,genotype_prefix,".bim")
  new_file_name_bim <- paste0(genotype_path,genotype_prefix,"_",chr,".bim")
  system(paste("cp ",old_file_name_bim,new_file_name_bim))

  old_file_name_fam <- paste0(genotype_path,genotype_prefix,".fam")
  new_file_name_fam <- paste0(genotype_path,genotype_prefix,"_",chr,".fam")
  system(paste("cp ",old_file_name_fam,new_file_name_fam))

  old_file_name_cov <- paste0(genotype_path,genotype_prefix,".cov")
  new_file_name_cov <- paste0(genotype_path,genotype_prefix,"_",chr,".cov")
  system(paste("cp ",old_file_name_cov,new_file_name_cov))

}

remove_SKAT_files_per_chr <- function(genotype_path, genotype_prefix)  {
  old_file_name_bed <- paste0(genotype_path,genotype_prefix,".bed")
  new_file_name_bed <- paste0(genotype_path,genotype_prefix,"_",chr,".bed")
  system(paste("rm ",old_file_name_bed,new_file_name_bed))

  old_file_name_bim <- paste0(genotype_path,genotype_prefix,".bim")
  new_file_name_bim <- paste0(genotype_path,genotype_prefix,"_",chr,".bim")
  system(paste("rm ",old_file_name_bim,new_file_name_bim))

  old_file_name_fam <- paste0(genotype_path,genotype_prefix,".fam")
  new_file_name_fam <- paste0(genotype_path,genotype_prefix,"_",chr,".fam")
  system(paste("rm ",old_file_name_fam,new_file_name_fam))

  old_file_name_cov <- paste0(genotype_path,genotype_prefix,".cov")
  new_file_name_cov <- paste0(genotype_path,genotype_prefix,"_",chr,".cov")
  system(paste("rm ",old_file_name_cov,new_file_name_cov))


}

combine_skat_results <- function(genotype_prefix, result_folder) {
  # Define the output file name
  output_file <- paste0(genotype_prefix, "_combined_association.tsv")

  # Execute shell commands to combine SKAT files
  system(paste('echo -e "Gene_name\tGene_chromosome\tRegion_start\tRegion_end\tQ_test\tpvalue" >', file.path(result_folder, output_file)))
  system(paste('for file in', file.path(result_folder, paste0("*", genotype_prefix, "*")), '; do if [ "$file" !=', file.path(result_folder, output_file), ']; then tail -n +2 "$file" >>', file.path(result_folder, output_file), '; fi; done'))

  # Optional: Print message indicating completion
  cat("SKAT association results combined successfully.\n")
}



# Required files
#Plink formatted bim
#Plink formatted bed
#Plink formatted fam

#Gene regions file

#Weight files

# EXAMPLE USAGE
#Rscript skat.R genotype_prefix gene_regions.csv combined.csv "EDAS/EUR_SKAT/" "phylo" 10



# CODE STARTS RUNNING FROM HERE
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6 ) {
  cat("Usage: Rscript skat_analysis.R <genotype_prefix> <gene_regions_file> <weights_file> <genotype_path> <weights_type> <chromosome> \n")
  quit(status = 1)
}

# Extract command-line arguments
genotype_prefix <- args[1]
gene_regions_file <- args[2]
weights_file <- args[3]
genotype_path <- args[4]
weights_type <- args[5]
chr <- args[6]

chr <- gsub("chr", "", chr)

#work on code for result folder
result_folder <- "result_folder"
prepare_SKAT_files_per_chr(genotype_path, genotype_prefix)
print("done creating skat genotype")
genotype_prefix <- paste0(genotype_prefix,"_",chr)
extract_weights_for_snvs_and_skat_chr(genotype_prefix, gene_regions_file, weights_file, genotype_path, weights_type, result_folder, chr)

#rm skatfiles.

cat("Done with SKAT analysis\n")
remove_SKAT_files_per_chr(genotype_path,genotype_prefix)
