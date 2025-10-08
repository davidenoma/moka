#!/usr/bin/env bash
# Wrapper to run Snakemake on ARC using a user-space micromamba + fresh conda
# Submit with:  sbatch run_smk.sbatch   OR run directly:  bash run_smk.sh [extra snakemake args]
set -euo pipefail

# --- Configurable knobs (override via environment) ---
: "${MAMBA_EXE:=$HOME/moka/bin/micromamba}"
: "${MAMBA_ROOT_PREFIX:=$HOME/micromamba}"
: "${SNAKE_CONDA_PREFIX:=$HOME/.snakemake/conda}"

# Determine threads: prefer SLURM_CPUS_PER_TASK, else default to 22
J="${SLURM_CPUS_PER_TASK:-22}"

# Ensure micromamba hook is loaded in non-interactive shells
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"

# Activate the driver env that includes snakemake + a new conda (>=24.7.1)
if ! micromamba activate snk; then
  echo "[ERROR] Failed to activate 'snk' environment. Did you run the install steps in the README?"
  exit 1
fi

# Ensure conda-prefix exists to avoid surprises on first run
mkdir -p "$SNAKE_CONDA_PREFIX"

echo "[INFO] Using $(snakemake --version) with $(conda --version)"
echo "[INFO] Threads: $J  |  Conda prefix: $SNAKE_CONDA_PREFIX"
echo "[INFO] Working dir: $(pwd)"

# Run Snakemake (forward any extra args provided to this script)
snakemake -j"$J" --use-conda --conda-prefix "$SNAKE_CONDA_PREFIX" "$@"
