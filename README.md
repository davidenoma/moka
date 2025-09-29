[![Snakemake](https://img.shields.io/badge/snakemake-8.15.2-brightgreen.svg)](https://snakemake.github.io)

<div align="center">
  <img src="https://github.com/user-attachments/assets/6d7d5099-aac4-44e2-a3d2-eacf6921a395" alt="image" width="500">
</div>

# 🌉 Multi-omics bridged Kernel Association test (MOKA) Pipeline
MOKA implements a Snakemake pipeline to automate data bridge kernel-based association tests. 
This pipeline offers flexibility of GWAS analysis & visualizations with different multi-omics variant specific weights.
<div align="center"> 
<img width="700" height="700" alt="moka-figure" src="https://github.com/user-attachments/assets/253b74bd-435e-4405-97b2-cf63ad73a3ef" />

</div>
Publication available at: https://www.medrxiv.org/content/10.1101/2025.07.06.25330974v1 

## 🚀 Usage
To run the moka pipeline:

1.Minimal data inputs

-GWAS genotype files in PLINK format (bed, bim & fam)

-Variant specific weights for each SNP ('SNP_ID, CHROMOSOME, POSITION, WEIGHT)

2. Install **Snakemake**
   
   - [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
     ```bash
      conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
      conda activate snakemake
      snakemake --help
     ```
3. Download and install moka

```bash
git clone https://github.com/davidenoma/moka
cd moka
```

4. Configure the pipeline parameters in the `config.yaml` file.

## 📚 Rules

### Rule: moka association_test
- **Input:** Preprocessed genotype data and weight files.
- **Output:** Results of association tests.
```bash
snakemake --cores <num_cores> --use-conda
```
    <num_cores> are the number of cores to use e.g. 8 

### Note
The initial run would take some time as the software installs the core dependencies requires from workflow/envs/moka.yaml.


### Rule: merge_results
- **Input:** Individual association test results.
- **Output:** Merged association test results.
```bash
snakemake --cores 1 merge_moka_results --use-conda
```

### Rule: visualize_results
- **Input:** Merged association test results.
- **Output:** Manhattan plots with visual representations of association test results.

```bash
snakemake --cores 1 manhattan_plots --use-conda
```

### Rule: go_analysis
- **Input:** Merged association test results.
- **Output:** GO analysis results.

```bash
snakemake --cores 1 go_analysis --use-conda
```

### Rule: kegg_pathway_analysis
- **Input:** Merged association test results.
- **Output:** KEGG pathway analysis results.

```bash
snakemake --cores 1 kegg_pathway_analysis --use-conda
```

### Rule: annotate_results
- **Input:** Merged association test results.
- **Output:** Annotated association test results with DisGeNet database

```bash
snakemake --cores 1 disgenet_annotation_005 --use-conda
```

### Rule: generate_gene_regions
- **Documentation:** Generates gene region files with specified flanking size from GFF3 annotation for gene-based association testing. Uses `config.flank_size` from the config file.

**How to execute:**
```bash
snakemake --cores 1 generate_gene_regions --use-conda
```

### Dependencies
All required Python and R packages, as well as other software dependencies, are specified in the provided `envs.yaml` files in the `workflow/envs/` directory. These environments are automatically created and managed by Snakemake when you use the `--use-conda` flag. You do not need to install packages manually.

For more details, see the environment YAML files in `workflow/envs/`.


### Input file format
- **Data Files:** Plink https://www.cog-genomics.org/plink/1.9/  format genotyped BIM, BED & FAM files [!required]
- Multi-omics **Bridge weights.csv** file (SNP_ID,Chromosome,Position,Weight) [!required for moka]
- Gene regions file (provided in GRCh38 or hg38) is generated from Ensembl GFF3 annotation files, e.g., Homo_sapiens.GRCh38.115.gff3.gz from Ensembl release 115.
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

## Publication reference 
MOKA: A pipeline for multi-omics bridged SNP-set kernel association test
https://www.medrxiv.org/content/10.1101/2025.07.06.25330974v1

## 🐳 Docker

You can run the MOKA pipeline using the official Docker image:

```bash
docker pull davidenoma/moka-gwas
```

For usage instructions and examples, see [https://hub.docker.com/r/davidenoma/moka-gwas](https://hub.docker.com/r/davidenoma/moka-gwas).
