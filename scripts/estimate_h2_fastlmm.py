#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from fastlmm.association import heritability_spatial_correction

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geno", required=True)
    parser.add_argument("--pheno", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    X = np.load(args.geno)
    Y = pd.read_csv(args.pheno, header=None).values.flatten()

    result = heritability_spatial_correction.h2(Y, X)
    h2 = result["h2"]
    with open(args.out, "w") as f:
        f.write(f"{h2:.6f}\n")

if __name__ == "__main__":
    main()
