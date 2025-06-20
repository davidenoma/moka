import os
import sys
import pandas as pd
from pysnptools.snpreader import Bed
from pysnptools.snpreader import Pheno
from fastlmm.association import single_snp
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Perform MLMA using FaST-LMM with PLINK SNP data.")
    parser.add_argument("--snp_prefix", required=True, help="Absolute path to the PLINK-formatted SNP prefix file")

    # Parse the arguments
    args = parser.parse_args()
    snp_data_loc = args.snp_prefix

    # Load PLINK Bed file
    bedfile = Bed(snp_data_loc)

    # Load the phenotype data from the .fam file
    phenotype_file = snp_data_loc + ".fam"
    phenotypes = pd.read_csv(phenotype_file, sep=' ', header=None)
    phenotypes.columns = ["FID", "IID", "father", "mother", "sex", "phenotype"]  # Adjust column names if necessary

    # Create Pheno object
    pheno_dict = {
        'header': ['phenotype'],
        'vals': phenotypes[['phenotype']].values,
        'iid': phenotypes[['FID', 'IID']].astype(str)
    }
    pheno = Pheno(pheno_dict)

    # Perform MLMA using FaST-LMM
    # print("Starting MLMA analysis with FaST-LMM... to get HERITABILITY ")
    # Perform analysis and extract only null h2
    results_df = single_snp(test_snps=bedfile, pheno=pheno, leave_out_one_chrom=False)
    # h2 = results_df['Nullh2'].iloc[0]  # Get first value since null h2 is constant
    h2 = results_df['Nullh2'].iloc[0]
    if pd.isna(h2):
        h2 = 0.5
    print(f"{h2:.15f}")
if __name__ == "__main__":
    main()
