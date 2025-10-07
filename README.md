[![Snakemake](https://img.shields.io/badge/snakemake-8.15.2-brightgreen.svg)](https://snakemake.github.io)

<div align="center">
  <img src="https://github.com/user-attachments/assets/6d7d5099-aac4-44e2-a3d2-eacf6921a395" alt="image" width="500">
</div>


# 🌉 Multi-omics bridged Kernel Association test (MOKA) Pipeline
MOKA implements a Snakemake pipeline to automate data bridge kernel-based association tests.  
This pipeline offers flexibility of GWAS analysis & visualizations with different multi-omics variant-specific weights.



<div align="center"> 
<img width="700" height="700" alt="moka-figure" src="https://github.com/user-attachments/assets/e9be1f99-1000-4cdc-a7c7-4acf13890bc0" />
</div




## 🚀 Usage
### 1. Prepare data inputs
- **GWAS genotype files** in **PLINK format** (`.bed`, `.bim`, `.fam`) [!required]
- **Variant-specific weights** in CSV format (SNP_ID,Chromosome,Position,Weight) [!required for moka]
- **Gene regions file** (provided in GRCh38 or hg38) , generated from Ensembl GFF3 annotations (e.g., Homo_sapiens.GRCh38.115.gff3.gz, Ensembl Release 115).
- **DisGeNET gene disease database** reference file ( If disease external validation needed)

### 2. Install **Snakemake**
Follow the [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):
```bash
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
conda activate snakemake
snakemake --help
```

### 3. Download and install **moka**
```bash
git clone https://github.com/davidenoma/moka
cd moka
```

### 4. Configure the pipeline parameters in the `config.yaml` file
This configuration controls paths to inputs and analysis settings:
- ***genotype_prefix:*** Prefix for the genotype data files (without file extensions).
- ***weights_type:*** Text descriptor indicating the source or type of biologically informed functional weights used in the analysis.
- ***genotype_file_path:*** Directory path to the genotype data files.
- ***weight_file:*** File path to the functional weight file applied in the association tests.
- ***disgenet_reference_file:*** Path to the DisGeNET reference file containing disease-specific gene–disease associations https://disgenet.org [For gene disease associations only!].
- ***spectral decomposition:*** Boolean flag to enable spectral (eigenvalue) decomposition and transformation of genotype and phenotype matrices. Default: TRUE.
- ***is_binary:*** Boolean flag specifying the phenotype type (TRUE for binary traits; FALSE for quantitative traits). Default: TRUE. 
- ***Plink:*** Path to the PLINK executable, e.g."~/software/plink".

### 5. Running MOKA
We provide a demo example with configuration located at ./config/config.yaml. The pipeline executes the following steps:

**Step 1. Kernel-based association testing**    
  - Integrates GWAS genotype data with the provided weights.
  - For the weight file, supports diverse data sources derived SNP-level weights. 
  - Performs SNP-set kernel-based association tests to model the **joint effect of multiple variants**.
  - Execute one chromosome at a time.
  - Optionally applies **decorrelation** to account for population structure or relatedness.
  
***Note:*** The initial run would take some time as the software installs the core dependencies requires from workflow/envs/moka.yaml. 

**Input:** Preprocessed genotype data (./genotype_data/test_geno.fam test_geno.bim test_geno.bed ) and weight files (./weights/test_geno_weights.csv).

**Output:** Results of association tests under ./result_folder/
  
```bash
snakemake --cores <num_cores> --use-conda
```
    <num_cores> are the number of cores to use e.g. 8 

**Step 2. Merge results from all chromosomes to a single file**  
**Input:** Individual association test results.

**Output:** Merged association test results under ./result_folder/
  
```bash
snakemake --cores 1 merge_moka_results --use-conda
```

**Step 3. Functional enrichment analysis**  
- Significant genes from association testing are assessed for:  
  - **KEGG pathway enrichment**  
  - **Gene Ontology (GO) enrichment**  
- Helps identify biological processes and pathways underlying GWAS signals.

**Input:** Merged association test results.

**Output:** GO analysis results under ./output_plots/

```bash
snakemake --cores 1 go_analysis --use-conda
```

**Output:** KEGG pathway analysis results under ./output_plots/

```bash
snakemake --cores 1 kegg_pathway_analysis --use-conda
```

**Step 4. Visualization**  
- Generates publication-ready visual summaries: **Manhattan plots**  

**Input:** Merged association test results.

**Output:** Manhattan plots with visual representations of association test results under ./output_plots/
  
```bash
snakemake --cores 1 manhattan_plots --use-conda
```

**Step 5. External validation**  
- Cross-references significant genes against the **DisGeNET database**.  
- Reports a **validation ratio** of overlapping associated genes, strengthening interpretation of GWAS findings.  

**Input:** Merged association test results.

**Output:** Annotated association test results with DisGeNet database under ./output_plots/

```bash
snakemake --cores 1 disgenet_annotation_005 --use-conda
```

**Additional function: generate_gene_regions**  
- Generates gene region files with specified flanking size from GFF3 annotation for gene-based association testing. Uses `config.flank_size` from the config file.

**How to execute:**
```bash
snakemake --cores 1 generate_gene_regions --use-conda
```

### Dependencies
All required Python and R packages, as well as other software dependencies, are specified in the provided `envs.yaml` files in the `workflow/envs/` directory. These environments are automatically created and managed by Snakemake when you use the `--use-conda` flag. You do not need to install packages manually.

For more details, see the environment YAML files in `workflow/envs/`.

### Liftover protocol 
You much lift over to GRCh38 format check here: Liftover GWAS: [https://github.com/davidenoma/LiftOver] 


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
