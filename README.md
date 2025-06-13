

# 🌉 Multi-omics bridged SNP-set kernel association test (MOKA) Pipeline
<div align="center">
  <img src="https://github.com/user-attachments/assets/6d7d5099-aac4-44e2-a3d2-eacf6921a395" alt="image" width="500">
</div>

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

## Introduction
Multi-Omics bridged SNP-set Kernel Association test (MOKA) aims to implement a Snakemake pipeline to automate data bridge kernel-based association tests. This pipeline offers flexibility for various types of association studies with different bridge weights.

## 🚀 Usage
To run the BAS pipeline:

1.Minimal data inputs
   1. GWAS genotype files in plink format
   2. Variant specific weights for each SNP ('SNP_ID, CHROMOSOME, POSITION, WEIGHT)
2. Install **Snakemake** 
   - [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
     ```bash
     conda install -n base -c conda-forge mamba
     mamba create -c conda-forge -c bioconda -n snakemake
     mamba activate snakemake
     snakemake --help
     ```
3. Install Plink, Python & R (Rscript configure)  [!important Check Dependencies section](https://github.com/davidenoma/moka/blob/main/README.md#dependencies)
4. Download and install moka
```bash
git clone https://github.com/davidenoma/moka
cd moka
```
5. Configure the pipeline parameters in the `config.yaml` file.
6. Execute the pipeline using the software:


## 📚 Rules
### Rule: moka association_test
- **Input:** Preprocessed genotype data and weight files.
- **Output:** Results of association tests.
```bash
snakemake --cores <num_cores>
```

   If you do not have all the dependencies with Python and R you can get it configured on conda, utilize with:

   ```bash
   snakemake --cores <num_cores> --use-conda
   ```
However, some R packages are not available to best to be installed R package manager.

### Rule: merge_results
- **Input:** Individual association test results.
- **Output:** Merged association test results.
- 
```bash
snakemake --cores 1 merge_moka_results
```


### Rule: annotate_results
- **Input:** Merged association test results.
- **Output:** Annotated association test results with DisGeNet database

```bash
snakemake --cores 1 disgenet_annotation_005
```

### Rule: visualize_results
- **Input:** Merged association test results.
- **Output:** Manhattan plots with visual representations of association test results.

```bash
snakemake --cores 1 manhattan_plots
```

### Rule: go_analysis
- **Input:** Merged association test results.
- **Output:** GO analysis results.

```bash
snakemake --cores 1 go_analysis
```

### Rule: kegg_pathway_analysis
- **Input:** Merged association test results.
- **Output:** KEGG pathway analysis results.

```bash
snakemake --cores 1 kegg_pathway_analysis
```
### Rule:  Skat test with linear kernel
- **Input:** Genotype
- **Output:** results for association mapping, folder: output_association/
```bash
snakemake --cores 22 skat

```
## Dependencies

### Software
**They must be configured on your path**
- **Snakemake (8.0.1+)**
- **R(4.2.0+)**
- **Python (3.9+)**
- **PLINK (1.9+)**: [https://www.cog-genomics.org/plink/1.9/]
- **Rscript**

### Python Packages
- **FaST-LMM**  Factored Spectrally Transformed Linear Mixed Models, is a program for performing genome-wide association studies (GWAS) on datasets of all sizes 
- **PySnpTools**  PySnpTools is a library for reading and manipulating genetic data.
  
```Python
pip install pysnptools fastlmm
```
### R Packages

- **manhattan:** R package for creating manhattan plots, commonly used in genome-wide association studies (GWAS).
- **SKAT:** R package for SKAT (Sequence Kernel Association Test) which is a powerful gene-based association test.
- **QQMAN:** R package for creating QQ (Quantile-Quantile) plots, commonly used in GWAS to assess whether observed p-values deviate from the expected distribution under the null hypothesis.
- **GGPLOT:** R package for creating highly customizable plots and graphics.
- **gprofiler2:** R package for gene set enrichment analysis (GO analysis).
- **pathfindR:** R package for pathway analysis, including KEGG pathway analysis.

Installation steps:

```R
install.packages(c("BiocManager","SKAT","ggplot2"))
BiocManager::install(c( "gprofiler2", "pathfindR","manhattan","qqman"))
```

### Other Software
- **Parallel:** Linux Parallel GNU : https://www.gnu.org/software/parallel/
```bash
apt install parallel #linux or WSL windows
brew install parallel #macos
```
### Input file format
- **Data Files:** Plink https://www.cog-genomics.org/plink/1.9/  format genotyped BIM, BED & FAM files [!required]
- Multi-omics **Bridge weights.csv** file (SNP_ID,Chromosome,Position,Weight) [!required for moka]
- Gene regions file provied in GRCh38 or hg38. (Genome Research Consortium Human Build 38)
- DisGeNET gene disease database reference file ( If disease external validation needed)

### Liftover protocol 
You much lift over to GRCh38 format check here: Liftover GWAS: [https://github.com/davidenoma/LiftOver] 

## 📋 Configuration
- **genotype_prefix:** Prefix for genotype data files.
- **weights_type:** Text string for type of bridge weights to be used e.g. "eqtl", "imaging"
- **genotype_file_path:** Path to genotype data files.
- **weight_file:** Path to weight files used for association tests.
- **disgenet_reference_file:** External disease database specific gene-disease associations from https://disgenet.org [For gene disease associations only!]
- **spectral decomposition:** Flag for performation decomposition and transformation of genotype and phenotype, default: TRUE
- **is_binary:** Flag for binary/ quantitative trait, default: TRUE 
- **Plink:** Path to plink installation e.g. "~/software/plink"
  
## 📖 Additional Information
For more information on the MOKA pipeline and its usage, refer to the documentation provided in the repository or contact the project maintainers.
david.enoma@ucalgary.ca 

