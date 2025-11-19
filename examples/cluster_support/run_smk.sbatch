#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --cpus-per-task=22
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/smk_%j.out
#SBATCH --error=logs/smk_%j.err
# #SBATCH --partition=compute          # uncomment & set if your cluster requires a partition

set -euo pipefail

# Make sure logs dir exists
mkdir -p logs

# Pin micromamba in the batch environment (don’t rely on ~/.bashrc)
export MAMBA_EXE="$HOME/moka/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export SNAKE_CONDA_PREFIX="$HOME/.snakemake/conda"

# Use the submission directory as the project root (default SLURM behavior, but explicit here)
cd "$SLURM_SUBMIT_DIR"

# Run the wrapper; pass additional Snakemake args by editing this line or use --wrap variant
bash ./run_smk.sh
