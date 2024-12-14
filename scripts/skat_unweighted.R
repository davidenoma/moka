# Load required R libraries
library(SKAT)

# Define the function to perform the SKAT test
perform_skat_test <- function(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, result_file) {
  tryCatch({
    # Rest of the function remains the same
    genotype_prefix <- paste0(genotype_path, genotype_prefix)
    genotype_bed <- paste0(genotype_prefix, ".bed")
    
    # Read the SNP list file
    snp_list <- read.table(paste0(genotype_path,"snp_list_", chr,".txt"), header = FALSE, col.names = "SNP")
    
    system(paste(
      "plink --bfile", 
      genotype_prefix, 
      "--extract", 
      paste0(genotype_path,"snp_list_", chr,".txt"), 
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
    FAM <- Read_Plink_FAM(genotype_fam, Is.binary = TRUE) 
    y <- FAM$Phenotype
    
    obj <- SKAT_Null_Model(y ~ 1, out_type = "D")
    genotype_ssd <- paste0(genotype_prefix, ".ssd")
    SSD.INFO <- Open_SSD(genotype_ssd, paste0(genotype_prefix, ".info"))
    id <- 1
    # Z <- Get_Genotypes_SSD(SSD.INFO, id)
    SetID <- SSD.INFO$SetInfo$SetID[id]
    # Use 'gene_snps' as the subset of SNPs for SKAT test
    skat_test <- SKAT.SSD.OneSet(SSD.INFO,SetID,obj, kernel = "linear.weighted", obj.SNPWeight=NULL)
    
    # Append SKAT test results to result file
    ss2 <- c(gene_name, gene_chromosome, region_start, region_end, toString(skat_test$Q), toString(skat_test$p.value))
    write(ss2, file = result_file, append = TRUE, ncol = 6, sep = '\t')
    
    
    # Clean up temporary files
    # file.remove(paste0(genotype_prefix, "_SKAT.bed"), paste0(genotype_prefix, "_SKAT.bim"), paste0(genotype_prefix, "_SKAT.fam"), paste0(genotype_prefix, ".setid"), paste0(genotype_prefix, ".ssd"), paste0(genotype_prefix, ".info"))
    file.remove(paste0(genotype_prefix, "_SKAT.bed"))
    file.remove (paste0(genotype_prefix, "_SKAT.bim"))
    file.remove(paste0(genotype_prefix, "_SKAT.fam"))
    file.remove(paste0(genotype_prefix, "_SKAT.log"))
    file.remove ( paste0(genotype_prefix, ".setid"))
    file.remove(paste0(genotype_prefix, ".ssd"))
    file.remove(paste0(genotype_prefix, ".info"))
    file.remove (paste0(genotype_prefix, ".ssd_LOG.txt"))
    file.remove (paste0(genotype_prefix, ".info.TEMP.txt"))
    
  }, error = function(e) {
    cat("Error in SKAT analysis for gene:", gene_name, "\n")
    cat("Error message:", e$message, "\n")
  })
}


# Perform SKAT for a specific chromosome
snvs_and_skat_chr <- function(genotype_prefix, gene_regions_file, genotype_path, chr, result_folder) {

  # Read the gene regions files
  gene_regions <- read.csv(gene_regions_file, header = TRUE)
  
  # Create results folder
  system(paste0("mkdir -p ",result_folder))
  result_file <- file.path(result_folder, paste0(genotype_prefix, "_result_chr_", chr, ".txt"))
  
  ss <- c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue")
  write(ss, file = result_file, ncol = 6, sep = '\t')
  
  gene_count <- 0 
  gene_regions <- subset(gene_regions, Chromosome == chr)
  
  #Read in the bim file
  snp_data <- read.csv(paste0(genotype_path, genotype_prefix,".bim"),sep = "\t", header = FALSE, col.names = c("chr","snp_id","pos_cm","loc","alt","ref"))

  snp_data <- subset(snp_data,chr == chr)
  
  
  # Iterate through gene regions
  for (i in 1:nrow(gene_regions)) {
    gene_name <- gene_regions$GeneName[i]
    region_start <- gene_regions$Start[i]
    region_end <- gene_regions$Stop[i]
    gene_chromosome <- gene_regions$Chromosome[i]
    
    # Select SNPs within the gene region from the snp data or bim file.
    gene_snps <- subset(snp_data, loc >= region_start & loc <= region_end)
    
    # Perform SKAT test if there are SNPs in the gene region
    if (nrow(gene_snps) > 0) {
      # Create a temporary SNP list file
      write.table(gene_snps$snp_id, paste0(genotype_path,"snp_list_", chr,".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      print(paste("Gene number: ", gene_count))
      # Call the SKAT test function with selected SNPs
      perform_skat_test(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, result_file)
      gene_count <- gene_count + 1
      # Remove the temporary SNP list file
      file.remove(paste0(genotype_path,"snp_list_", chr,".txt"))
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
  
}

remove_SKAT_files_per_chr <- function(genotype_path, genotype_prefix)  {
  old_file_name_bed <- paste0(genotype_path,genotype_prefix,".bed")
  new_file_name_bed <- paste0(genotype_path,genotype_prefix,"_",chr,".bed")
  system(paste("rm ",old_file_name_bed,new_file_name_bed))
  
  old_file_name_bim <- paste0(genotype_path,genotype_prefix,".bim")
  new_file_name_bim <- paste0(genotype_path,genotype_prefix,"_",chr,".bim")
  system(paste("cp ",old_file_name_bim,new_file_name_bim))
  
  old_file_name_fam <- paste0(genotype_path,genotype_prefix,".fam")
  new_file_name_fam <- paste0(genotype_path,genotype_prefix,"_",chr,".fam")
  system(paste("cp ",old_file_name_fam,new_file_name_fam))
  
}
combine_skat_results <- function(genotype_prefix, result_folder) {
  # Define the output file name
  output_file <- file.path(result_folder, paste0(genotype_prefix, weights_type,"_combined_association.tsv"))
  
  # Check if the result folder exists
  if (!dir.exists(result_folder)) {
    stop("Result folder does not exist.")
  }
  
  # Write the header to the combined results file
  write.table(matrix(c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue"), nrow=1), 
              file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # Combine SKAT results files
  result_files <- list.files(result_folder, pattern = paste0(genotype_prefix, "_", weights_type, "_result_chr_"), full.names = TRUE)
  for (file in result_files) {
    if (file != output_file) {
      # Append the content of each result file, skipping the header
      write.table(read.table(file, header = TRUE, sep = "\t"), 
                  file = output_file, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = "\t")
    }
  }
  
  # Optional: Print message indicating completion
  cat("SKAT association results combined successfully into:", output_file, "\n")
}

# CODE STARTS RUNNING FROM HERE
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 ) {
  cat("Usage: Rscript skat_analysis.R <genotype_prefix> <gene_regions_file> <genotype_path> <chromosome> <output_folder> \n")
  quit(status = 1)
}




# Extract command-line arguments
genotype_prefix <- args[1]
gene_regions_file <- args[2]
genotype_path <- args[3]
chr <- args[4]
output_folder <- args[5]
chr <- gsub("chr", "", chr)
#work on code for result folder
result_folder <- output_folder





prepare_SKAT_files_per_chr(genotype_path, genotype_prefix)
print("done creating skat genotype")
genotype_prefix <- paste0(genotype_prefix,"_",chr)
snvs_and_skat_chr(genotype_prefix, gene_regions_file, genotype_path, chr, result_folder)
# combine_skat_results(genotype_prefix, result_folder)
#rm skatfiles.


cat("Done with SKAT analysis\n")
#TODO 
