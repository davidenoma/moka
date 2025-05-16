import argparse
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed, Pheno, SnpData
from fastlmm.association import heritability_spatial_correction
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os

def run_pca(snp_prefix, out_file, n_components=2):
    print(f"📥 Reading genotype matrix from: {snp_prefix}")
    bed = Bed(snp_prefix, count_A1=False)
    geno = bed.read().val
    iids = bed.iid

    print(f"✅ Loaded genotype matrix: {geno.shape}")
    print("⚙️  Running PCA on standardized genotypes...")
    X_std = StandardScaler().fit_transform(geno)
    pcs = PCA(n_components=n_components).fit_transform(X_std)

    df = pd.DataFrame(iids, columns=["FID", "IID"])
    for i in range(n_components):
        df[f"PC{i+1}"] = pcs[:, i]

    df.to_csv(out_file, sep="\t", index=False)
    print(f"✅ PCA spatial coordinates saved to: {out_file}")

def load_spatial_coords(coord_file):
    df = pd.read_csv(coord_file, sep="\t")
    coords = df[["PC1", "PC2"]].values
    spatial_iid = df[["FID", "IID"]].values
    return coords, spatial_iid

def estimate_heritability(snp_prefix, coord_file, out_file):
    bed = Bed(snp_prefix, count_A1=False)

    # Load phenotype from FAM
    fam_file = snp_prefix + ".fam"
    fam = pd.read_csv(fam_file, sep=" ", header=None)
    fam.columns = ["FID", "IID", "father", "mother", "sex", "phenotype"]
    pheno_dict = {
        'header': ['phenotype'],
        'vals': fam[['phenotype']].values,
        'iid': fam[['FID', 'IID']].astype(str)
    }
    pheno = Pheno(pheno_dict)
    alpha_grid = [int(v) for v in np.logspace(np.log10(100), np.log10(1e10), 10)]
    # coords, spatial_iid = load_spatial_coords(coord_file)
    df = pd.read_csv(coord_file, sep="\t")
    coords = df[["PC1", "PC2"]].values
    spatial_iid = df[["FID", "IID"]]
    print("📊 Running heritability_spatial_correction...")
    df = heritability_spatial_correction(
        bed,
        spatial_coor=coords,
        spatial_iid=spatial_iid,
        alpha_list=alpha_grid,
        alpha_power=2,
        pheno=pheno,
        jackknife_count=100,
        permute_plus_count=10,
        permute_times_count=10,
        seed=42,
        allow_gxe2=True
    )

    df.to_csv(out_file, sep="\t", index=False)
    print(f"✅ Heritability results written to {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Estimate heritability with spatial correction using PCA-based spatial coordinates.")
    parser.add_argument("--snp_prefix", required=True, help="Prefix for PLINK .bed/.bim/.fam files")
    parser.add_argument("--pca_out", default=None, help="Path to save PCA coordinates (optional)")
    parser.add_argument("--out", required=True, help="Output file for heritability results")
    args = parser.parse_args()

    # Define default PCA coord output
    pca_out_file = args.pca_out or f"{args.snp_prefix}_spatial_coords.tsv"

    # Step 1: Run PCA to get coordinates
    if os.path.exists(pca_out_file):
        print("PCA coordinates file exists, loading...")
        coords_df = pd.read_csv(pca_out_file, sep="\t")
    else:
        print("PCA coordinates file not found, running PCA...")
        run_pca(args.snp_prefix, pca_out_file)


    # Step 2: Run heritability estimation
    estimate_heritability(args.snp_prefix, pca_out_file, args.out)

if __name__ == "__main__":
    main()
