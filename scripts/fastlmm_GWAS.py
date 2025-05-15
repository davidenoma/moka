import argparse
import pandas as pd
from pysnptools.snpreader import Bed
from pysnptools.snpreader import Pheno
from fastlmm.association import single_snp

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snp_prefix", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    bedfile = Bed(args.snp_prefix)
    phenotypes = pd.read_csv(args.snp_prefix + ".fam", sep=' ', header=None)
    phenotypes.columns = ["FID", "IID", "father", "mother", "sex", "phenotype"]

    pheno_dict = {
        'header': ['phenotype'],
        'vals': phenotypes[['phenotype']].values,
        'iid': phenotypes[['FID', 'IID']].astype(str)
    }
    pheno = Pheno(pheno_dict)

    results = single_snp(test_snps=bedfile, pheno=pheno)
    h2_est = results["h2"][0]
    with open(args.out, "w") as f:
        f.write(str(h2_est) + "\n")

if __name__ == "__main__":
    main()
