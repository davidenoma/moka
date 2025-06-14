# Load required R libraries
library(SKAT)
library(parallel)
# import
# Define the function to perform the SKAT test
perform_skat_test <- function(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, result_folder, result_file, is_binary = TRUE) {
  tryCatch({
    genotype_prefix <- paste0(genotype_path, genotype_prefix)
    genotype_bed <- paste0(genotype_prefix, ".bed")

    # Read the SNP list file
    snp_list <- read.table(paste0(genotype_path,"snp_list_", weights_type, chr,".txt"), header = FALSE, col.names = "SNP")

    system(paste(
      "plink --bfile",
      genotype_prefix,
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

    # Perform SKAT test based on phenotype type (binary or continuous)
    genotype_fam <- paste0(genotype_prefix, ".fam")
    FAM <- Read_Plink_FAM(genotype_fam, Is.binary = is_binary)
    y <- FAM$Phenotype

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
    skat_test <- SKAT.SSD.OneSet(SSD.INFO, SetID, obj, kernel = "linear.weighted", weights = gene_snps$Weight)

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
perform_skat_test_decomposition <- function(
  gene_name, gene_chromosome, region_start, region_end, gene_snps,
  genotype_prefix, genotype_path, result_folder, result_file, D, is_binary = TRUE
) {
  genotype_prefix <- paste0(genotype_path, genotype_prefix)

  tryCatch({
    message("Starting SKAT for gene: ", gene_name)

    # snp_list_file <- file.path(genotype_path, paste0("snp_list_", weights_type, chr, ".txt"))
    prefix_skat <- paste0(genotype_prefix, "_SKAT")

    # ----- Step 0: Create subset binary PLINK files -----
    system(paste(
      "plink --bfile", genotype_prefix,
      "--extract", paste0(genotype_path,"snp_list_", weights_type, chr,".txt"),
      "--make-bed --allow-no-sex",
      "--out", prefix_skat
    ))

    # ----- Step 1: Generate .raw file for dosage matrix -----
    system(paste(
      "plink --bfile", prefix_skat,
      "--recode A --out", prefix_skat,
      "--allow-no-sex "
    ))
    print('raw file generated')



    # ----- Step 3: Read genotype matrix -----
    raw_file <- paste0(prefix_skat, ".raw")
    if (!file.exists(raw_file)) stop("No .raw file found.") else (print("raw file found"))
# Read the genotype data
  geno_df <- read.table(raw_file, header = TRUE, sep = " ")

  all_cols <- colnames(geno_df)
  snp_indices <- which(sapply(all_cols, function(cn) {
    any(startsWith(cn, paste0(gene_snps$SNP, "_")))
  }))

  if (length(snp_indices) == 0) {
    warning("No SNPs matched for gene: ", gene_name)
    return(NULL)
  }

    genotype_matrix <- as.matrix(geno_df[, snp_indices])

    # ----- Step 4: Read phenotype -----
    # fam_file <- file.path(genotype_path, paste0(genotype_prefix, ".fam"))
    fam <- read.table(paste0(genotype_prefix, ".fam"), header = FALSE)
    Y <- fam$V6
    Y <- scale(Y, center = TRUE, scale = FALSE)
    #     # Step 2: Generate GRM using GCTA
    # system(paste(
    #   "gcta64 --bfile", prefix_skat,
    #   "--make-grm-bin --out", prefix_skat
    # ))
    # ----- Step 2.5: Generate GCTA-compatible .pheno file from .fam -----
  pheno_file <- paste0(prefix_skat, ".pheno")
  fam <- read.table(paste0(genotype_prefix, ".fam"), header = FALSE)
  fam_pheno <- fam[, c(1, 2, 6)]  # FID, IID, PHENO
  write.table(fam_pheno, file = pheno_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  # Step 3: Estimate h² using GCTA REML
#   system(paste(
#     "gcta64 --grm", prefix_skat,
#     "--pheno", pheno_file, "--thread-num", 22,
#     "--reml --out", prefix_skat
#   ))
#   h2_file <- paste0(prefix_skat, ".hsq")
# # # Read the hsq file with header set to TRUE and fill parameter for safety
# h2_data <- read.table(h2_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# # # Extract the row for 'V(G)/Vp'
# h2_line <- h2_data[h2_data$Source == "V(G)/Vp", ]
# # # Convert the Variance column value to numeric
# h2 <- as.numeric(h2_line$Variance)
# cat(h2)



    cat('Reading GRM\n')
    X <- genotype_matrix
    # cat(str(X), "\n")
    # G <- grm
    # # G <- read_grm_plink(prefix_skat)
    #
    # # cat(str(G))
    # cat('Decomposing GRM\n')
    # # ----- Step 8: Decorrelate using GRM -----
    # eig <- eigen(G, symmetric = TRUE)
    #
    # U <- unname(eig$vectors)
    # S <- eig$values
    # D <- U %*% diag(1 / sqrt(h2 * S + 1))
    # cat(str(D), "\n")
    cat('Done decomposition\n')
    Y_star <- D %*% Y
    cat(str(Y_star), "\n")
    X_star <- D %*% X
    cat(str(X_star), "\n")
    # intercept <- D %*% rep(1, length(Y_star))
    # cat(str(intercept), "\n")

    # ----- Step 9: Run SKAT -----
    obj <- SKAT_Null_Model(Y_star ~ 1,out_type="C")
    skat_result <- SKAT(X_star, obj, kernel = "linear.weighted", weights = gene_snps$Weight)
    # cat(skat_result$Q, skat_result$p.value, "\n")
    # ----- Step 10: Write result -----
    ss2 <- c(gene_name, gene_chromosome, region_start, region_end,
             toString(skat_result$Q), toString(skat_result$p.value))
    write(ss2, file = result_file, append = TRUE, ncol = 6, sep = "\t")
    # cat('Written association result for gene:')
    # ----- Step 11: Clean up -----
# Remove SKAT temporary files, excluding the .rel file
# Remove SKAT temporary files, excluding .rel and GRM files created by PLINK (e.g., .grm, .grm.bin, .grm.id, .grm.N.bin)
files_to_remove <- Sys.glob(file.path(genotype_path, paste0(prefix_skat, ".*")))
files_to_remove <- files_to_remove[!grepl("\\.(rel|grm|grm\\.bin|grm\\.id|grm\\.N\\.bin)$", files_to_remove)]
unlink(files_to_remove)

raw_files <- Sys.glob(file.path(genotype_path, paste0(genotype_prefix, "*.raw")))
unlink(raw_files)
  }, error = function(e) {
    cat("Error in SKAT for gene:", gene_name, "\n")
    cat("Message:", e$message, "\n")
     cat("Full Error Details:\n")
  print(e)
  })
}





# Define a function to extract weights for SNVs and perform SKAT for a specific chromosome
extract_weights_for_snvs_and_skat_chr <- function(genotype_prefix, gene_regions_file, weights_file, genotype_path, weights_type, result_folder, chr, D,is_binary = TRUE) {
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
      write.table(gene_snps$SNP, paste0(genotype_path, "snp_list_", weights_type, chr, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      print(paste("Gene number: ", gene_count))
      # Call the SKAT test function with selected SNPs
          if (spectral_decorrelated) {
        perform_skat_test_decomposition(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, genotype_path, result_folder, result_file,D, is_binary)
      } else {
        perform_skat_test(gene_name, gene_chromosome, region_start, region_end, gene_snps, genotype_prefix, result_folder, result_file, is_binary)
      }
        gene_count <- gene_count + 1
      file.remove(paste0(genotype_path, "snp_list_", weights_type, chr, ".txt"))
    } else {
      print("Not found for gene region")
    }
  }

  print(paste("Done with SKAT analysis for chromosome", chr))
}
read_grm <- function(prefix) {
      grm_bin <- paste0(prefix, ".grm.bin")
      grm_id <- paste0(prefix, ".grm.id")
      grm_vals <- readBin(grm_bin, what = "numeric", n = 1e9, size = 4)
      ids <- read.table(grm_id)
      n <- nrow(ids)
      G <- matrix(0, n, n)
      k <- 1
      for (i in 1:n) {
        for (j in 1:i) {
          G[i, j] <- grm_vals[k]
          G[j, i] <- grm_vals[k]
          k <- k + 1
        }
      }
      return(G)
    }
# PLINK-based GRM generation function
read_grm_plink <- function(prefix) {
  system(paste(
    "plink --bfile", prefix, "--allow-no-sex",
    "--make-rel square --out", prefix
  ))
  grm_file <- paste0(prefix, ".rel")
  G <- as.matrix(read.table(grm_file, header = FALSE))
  return(G)
}

make_grm_bin_plink <- function(prefix){

    system(paste(
      "plink --bfile", prefix,
      "--make-grm-bin --out", prefix,
      "--allow-no-sex "
    ))
    print("GRM generated")
}
# Helper function to prepare files for SKAT per chromosome
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

# Remove SKAT files per chromosome
remove_SKAT_files_per_chr <- function(genotype_path, genotype_prefix)  {
  old_file_name_bed <- paste0(genotype_path, genotype_prefix, ".bed")
  system(paste("rm ", old_file_name_bed))

  old_file_name_bim <- paste0(genotype_path, genotype_prefix, ".bim")
  system(paste("rm ", old_file_name_bim))

  old_file_name_fam <- paste0(genotype_path, genotype_prefix, ".fam")
  system(paste("rm ", old_file_name_fam))

  new_file_name_log <- paste0(genotype_path, genotype_prefix, "_", "SKAT", ".log")
  system(paste("rm ", new_file_name_log))
}

# Combine SKAT results
combine_skat_results <- function(genotype_prefix, result_folder, weights_type) {
  output_file <- file.path(result_folder, paste0(genotype_prefix, weights_type, "_combined_association.tsv"))

  if (!dir.exists(result_folder)) {
    stop("Result folder does not exist.")
  }

  write.table(matrix(c("Gene_name", "Gene_chromosome", "Region_start", "Region_end", "Q_test", "pvalue"), nrow=1),
              file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  result_files <- list.files(result_folder, pattern = paste0(genotype_prefix, "_", weights_type, "_result_chr_"), full.names = TRUE)
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
if (length(args) < 6) {
  cat("Usage: Rscript skat_analysis.R <genotype_prefix> <gene_regions_file> <weights_file> <genotype_path> <weights_type> <chromosome> [is_binary] [spectral_decorrelated]\n")
  quit(status = 1)
}

genotype_prefix <- args[1]
gene_regions_file <- args[2]
weights_file <- args[3]
genotype_path <- args[4]
weights_type <- args[5]
chr <- args[6]
is_binary <- ifelse(length(args) >= 7, as.logical(args[7]), TRUE)
spectral_decorrelated <- ifelse(length(args) >= 8, as.logical(args[8]), TRUE)
# script_path <- commandArgs(trailingOnly = TRUE)[9]  # Assume the 9th argument is the python script
   if (spectral_decorrelated) {
# First check if h2 file exists
h2_file <- paste0(genotype_path, genotype_prefix, "_h2.txt")
if (!file.exists(h2_file)) {
  # Run FastLMM and save h2
  h2_raw <- system(
    paste("python", "scripts/estimate_h2_fastlmm.py", "--snp_prefix", paste0(genotype_path, genotype_prefix)),
    intern = TRUE
  )
  h2_numeric <- as.numeric(trimws(h2_raw))
  h2 <- h2_numeric[which.max(!is.na(h2_numeric))]  # Take last non-NA
  # Write h2 to file
  write(sprintf("%.15f", h2), file = h2_file)
  cat(sprintf("Calculated and saved h²: %.15f\n", h2))
} else {
  # Read existing h2
  h2 <- as.numeric(readLines(h2_file))
  cat(sprintf("Loaded existing h²: %.15f\n", h2))
}

# Check for GRM file
grm_file <- paste0(genotype_path, genotype_prefix, ".grm.bin")
if (!file.exists(grm_file)) {
  # Generate GRM if it doesn't exist
  cat("Generating GRM...\n")
  make_grm_bin_plink(paste0(genotype_path, genotype_prefix))
  grm <- read_grm(paste0(genotype_path, genotype_prefix))
  G <- grm
    # G <- read_grm_plink(prefix_skat)

    # cat(str(G))
    cat('Decomposing GRM\n')
    # ----- Step 8: Decorrelate using GRM -----
    eig <- eigen(G, symmetric = TRUE)

    U <- unname(eig$vectors)
    S <- eig$values
    D <- U %*% diag(1 / sqrt(h2 * S + 1))
}  else{
  cat("Loading existing GRM...\n")
  grm <- read_grm(paste0(genotype_path, genotype_prefix))
  G <- grm
    # G <- read_grm_plink(prefix_skat)

    # cat(str(G))
    cat('Decomposing GRM\n')
    # ----- Step 8: Decorrelate using GRM -----
    eig <- eigen(G, symmetric = TRUE)

    U <- unname(eig$vectors)
    S <- eig$values
    D <- U %*% diag(1 / sqrt(h2 * S + 1))
  # grm <- as.matrix(read.table(grm_file, header = FALSE))
}
   }
chr <- gsub("chr", "", chr)
result_folder <- "result_folder"
prepare_SKAT_files_per_chr(genotype_path, genotype_prefix)
print("done creating skat genotype")
genotype_prefix <- paste0(genotype_prefix, "_", chr)
#assign GRM TO d
extract_weights_for_snvs_and_skat_chr(genotype_prefix, gene_regions_file, weights_file, genotype_path, weights_type, result_folder, chr,D, is_binary)

#combine_skat_results(genotype_prefix, result_folder, weights_type)
remove_SKAT_files_per_chr(genotype_path, genotype_prefix)

cat("Done with SKAT analysis\n")
