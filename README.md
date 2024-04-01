# 🌉 Bridged Association Studies (BAS) Pipeline

## Introduction
Bridged Association Studies (BAS) aims to implement a Snakemake pipeline to automate data bridge kernel-based association tests. This pipeline offers flexibility for various types of association studies beyond traditional eQTL analysis.

## 🚀 Usage
To run the BAS pipeline:
1. Install **Snakemake** and any required dependencies.
2. Configure the pipeline parameters in the `config.yaml` file.
3. Execute the pipeline using the command:

snakemake --cores <num_cores>


## 📚 Rules
### Rule: obtain_weights
- **Input:** Genotype data and relevant metadata.
- **Output:** Weight files for association tests.

### Rule: preprocess_data
- **Input:** Raw data files.
- **Output:** Preprocessed data ready for association testing.

### Rule: association_test
- **Input:** Preprocessed data and weight files.
- **Output:** Results of association tests.

### Rule: merge_results
- **Input:** Individual association test results.
- **Output:** Merged association test results.

### Rule: annotate_results
- **Input:** Merged association test results.
- **Output:** Annotated association test results.

### Rule: visualize_results
- **Input:** Annotated association test results.
- **Output:** Visual representations of association test results.

## 🛠️ Dependencies
- **Software:** Snakemake, Python, R
- **Data Files:** Genotype data, metadata, additional resources required for association testing

## 📋 Configuration
- **genotype_prefix:** Prefix for genotype data files.
- **analysis_type:** Type of association test to perform (e.g., eQTL, GWAS).
- **genotype_file_path:** Path to genotype data files.
- **weight_file:** Path to weight files used for association tests.
- **metadata_file:** Path to metadata files required for preprocessing.
- **output_folder:** Output directory for results.

## 📖 Additional Information
For more information on the BAS pipeline and its usage, refer to the documentation provided in the repository or contact the project maintainers.
