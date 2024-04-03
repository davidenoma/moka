# 🌉 Bridged Association Studies (BAS) Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


## Introduction
Bridged Association Studies (BAS) aims to implement a Snakemake pipeline to automate data bridge kernel-based association tests. This pipeline offers flexibility for various types of association studies with different bridge weights.

## 🚀 Usage
To run the BAS pipeline:

1. Install **Snakemake** and any required dependencies.
   - [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
     ```bash
     conda install -n base -c conda-forge mamba
     mamba create -c conda-forge -c bioconda -n snakemake snakemake
     mamba activate snakemake
     snakemake --help
     ```
2. Install R dependencies 
3. Configure the pipeline parameters in the `config.yaml` file.
3. Execute the pipeline using the command:


## 📚 Rules
### Rule: association_test
- **Input:** Preprocessed data and weight files.
- **Output:** Results of association tests.

snakemake --cores <num_cores>

### Rule: merge_results
- **Input:** Individual association test results.
- **Output:** Merged association test results.

snakemake --cores all -R --until merge_skat_results

### Rule: annotate_results
- **Input:** Merged association test results.
- **Output:** Annotated association test results.

snakemake --cores all -R --until disgenet_annotation

### Rule: visualize_results
- **Input:** Annotated association test results.
- **Output:** Visual representations of association test results.

snakemake --cores all -R --until manhattan_plots

## Dependencies

### Software

- **Snakemake**
- **R**
- **Python**
- **PLINK**
- **Rscript**

### R Packages

- **manhattan:** R package for creating manhattan plots, commonly used in genome-wide association studies (GWAS).
- **SKAT:** R package for SKAT (Sequence Kernel Association Test) which is a powerful gene-based association test.
- **PARALLEL:** R package for parallel computing capabilities in R.
- **QQMAN:** R package for creating QQ (Quantile-Quantile) plots, commonly used in GWAS to assess whether observed p-values deviate from the expected distribution under the null hypothesis.
- **GGPLOT:** R package for creating highly customizable plots and graphics.

### Other

- **PLINK:** PLINK is a widely used software toolset for genome-wide association studies (GWAS) and analysis of DNA sequencing data.
- **Data Files:** Plink genotyped Bim, Bed & Fam files, bridge weights file (format is specified), gene regions file, DisGeNET reference file 


## 📋 Configuration
- **genotype_prefix:** Prefix for genotype data files.
- **weights_type:** Type of bridge weights to be used such as EQTL, Imaging, Conservation, Representation etc.
- **genotype_file_path:** Path to genotype data files.
- **weight_file:** Path to weight files used for association tests.
- **disgenet_reference_file:** Disease database specific weights from https://disgenet.org [For gene disease associations only!]
- **Plink:** Path to plink installation
- **sample size:** Sample size of GWAS Cohort
  
## 📖 Additional Information
For more information on the BAS pipeline and its usage, refer to the documentation provided in the repository or contact the project maintainers.
