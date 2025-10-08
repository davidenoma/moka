# README — Run Snakemake on ARC with your own Conda (no root, no cluster changes)

This guide shows you how to run:

```bash
snakemake -c 22 --use-conda
```

on a cluster where the **system conda is too old** (e.g., 23.9.0) and you **can’t** update it.
We’ll install a *self-contained* micromamba/conda in your **$HOME**, then make Snakemake use **that** fresh conda (≥ 24.7.1).

---

## Why these steps?

- **Snakemake shells out to `conda` internally** and requires **conda ≥ 24.7.1**.
  If your cluster’s `conda` is older, Snakemake fails even if you use `--use-conda`.
- **Micromamba** is a fast, zero-admin package manager you can put in `$HOME`.
  We create a clean env **`snk`** that contains **snakemake** *and* a **new `conda`** binary.
  When you activate `snk`, Snakemake finds that new `conda` first and everything “just works”.
- We keep **all per-rule environments** created by Snakemake in a predictable folder (e.g., `~/.snakemake/conda`) so they’re cached and reusable.

---

## TL;DR (copy–paste)

> Do this once per account (safe on ARC; no root needed):

```bash
# 0) Make micromamba available in THIS shell
export MAMBA_EXE="$HOME/moka/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$($MAMBA_EXE shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"

# 1) Create a clean env that bundles snakemake + a NEW conda
micromamba create -n snk -y -c conda-forge -c bioconda python=3.11 snakemake conda

# 2) Activate it (adds the new 'conda' to PATH)
micromamba activate snk

# 3) (Optional, set defaults once)
mkdir -p ~/.config/snakemake
cat > ~/.config/snakemake/config.yaml <<'YAML'
conda-prefix: /home/'$USER'/.snakemake/conda
YAML

# 4) Run your workflow
snakemake -c 22 --use-conda
```

If you **haven’t** downloaded micromamba yet, see Step A below.

---

## Step A — One-time micromamba install (no root)

Pick a folder (e.g., your project or `$HOME/moka`) and run:

```bash
# download & unpack micromamba into the current dir
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# initialize (writes a small block into ~/.bashrc)
./bin/micromamba shell init -s bash ~/micromamba
source ~/.bashrc
```

**Why:**

- Puts a tiny `micromamba` binary under your control.
- Sets `MAMBA_ROOT_PREFIX=~/micromamba` (where your user-space envs live).
- Keeps you isolated from the system conda.

> If your current shell isn’t login-style and didn’t pick up `~/.bashrc`, either `source ~/.bashrc` or run `exec $SHELL -l`.

---

## Step B — Create the “driver” environment (`snk`)

```bash
# ensure micromamba is active *in this shell*
export MAMBA_EXE="$HOME/moka/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$($MAMBA_EXE shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"

micromamba create -n snk -y -c conda-forge -c bioconda python=3.11 snakemake conda
micromamba activate snk
```

**Why include `conda` in this env?**  
Recent Snakemake deprecates alternative frontends and shells out to **`conda`**. By installing a **new** `conda` *inside* `snk`, Snakemake finds that first (≥ 24.7.1) and stops hitting the cluster’s outdated one.

**Sanity checks:**
```bash
conda --version   # must be 24.7.1 or newer
snakemake --version
```

---

## Step C — Configure where Snakemake stores envs (optional but recommended)

```bash
mkdir -p ~/.config/snakemake
cat > ~/.config/snakemake/config.yaml <<'YAML'
# Where Snakemake will create/restore per-rule environments from your env YAMLs
# Use $HOME for portability; you can point this to $SCRATCH for faster I/O on some clusters.
conda-prefix: /home/'$USER'/.snakemake/conda
YAML
```

Now you can simply run:

```bash
micromamba activate snk
snakemake -c 22 --use-conda
```

**Why:**

- Centralizes cache of per-rule envs in a stable location.
- Saves resolve time across reruns and projects.

> Tip: For project-local caches, use `conda-prefix: .snakemake/conda` in a **project** config file instead.

---

## Day-to-day usage

From any new shell:

```bash
# Load your user-space micromamba
export MAMBA_EXE="$HOME/moka/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$($MAMBA_EXE shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"

# Activate the driver env
micromamba activate snk

# Run your pipeline
snakemake -c 22 --use-conda
```

Common flags you might add:

- `--rerun-incomplete` to rebuild incompletes.
- `--use-singularity` if your rules use containers.
- `--resources mem_mb=...` to respect memory limits.
- `--conda-cleanup-pkgs` to drop solver caches when space is tight.

---

## Wrapper + SLURM batch

### `run_smk.sh` (call directly or via `sbatch --wrap "bash run_smk.sh"`)

```bash
#!/usr/bin/env bash
set -euo pipefail

# --- Configurable knobs (override via environment) ---
: "${MAMBA_EXE:=$HOME/moka/bin/micromamba}"
: "${MAMBA_ROOT_PREFIX:=$HOME/micromamba}"
: "${SNAKE_CONDA_PREFIX:=$HOME/.snakemake/conda}"

# Resolve threads from SLURM_CPUS_PER_TASK if available, else default
J="${SLURM_CPUS_PER_TASK:-22}"

# Ensure micromamba hook is loaded in non-interactive shells
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"

# Activate the driver env that includes snakemake + a new conda
micromamba activate snk

# Create conda-prefix if missing
mkdir -p "$SNAKE_CONDA_PREFIX"

# Run Snakemake (forward any extra args)
snakemake -c "$J" --use-conda --conda-prefix "$SNAKE_CONDA_PREFIX" "$@"
```

### `run_smk.sbatch` (submit with `sbatch run_smk.sbatch`)

```bash
#!/bin/bash
#SBATCH --c ob-name=snakemake
#SBATCH --cpus-per-task=22
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/smk_%j.out
#SBATCH --error=logs/smk_%j.err
# #SBATCH --partition=compute          # uncomment & set if your cluster requires a partition

set -euo pipefail

# Make sure logs dir exists
mkdir -p logs

# Optional: pin micromamba in the batch environment (no reliance on ~/.bashrc)
export MAMBA_EXE="$HOME/moka/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export SNAKE_CONDA_PREFIX="$HOME/.snakemake/conda"

# Use the submission directory as the project root
cd "$SLURM_SUBMIT_DIR"

# Run the wrapper; pass any Snakemake args after a --
bash ./run_smk.sh
```

Usage examples:

```bash
# local interactive
bash run_smk.sh -k --rerun-incomplete

# submit to SLURM (uses #SBATCH cpus/time/mem)
sbatch run_smk.sbatch
```

---

## FAQ / Troubleshooting

**“Conda must be version 24.7.1 or later…”**  
You’re still picking up the system `conda`. Fix:
1. `micromamba activate snk`
2. `conda --version` (verify ≥ 24.7.1)
3. Run `snakemake -c 22 --use-conda` again.

If it persists on compute nodes, ensure your batch job initializes micromamba (the wrapper and sbatch do this explicitly).

**`snakemake: command not found`**  
You deactivated base (where snakemake lived) but didn’t install snakemake elsewhere. Activate `snk` or (once) run:
```bash
micromamba create -n snk -y -c conda-forge -c bioconda python=3.11 snakemake conda
micromamba activate snk
```

**`micromamba: command not found`**  
Run it by path or re-init the hook:
```bash
~/moka/bin/micromamba --version

export MAMBA_EXE="$HOME/moka/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$($MAMBA_EXE shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
```

**Creating per-rule envs is slow on shared storage**  
- Use a fast prefix (e.g., `$SCRATCH/.snakemake/conda`) in your config.
- Micromamba/conda will cache packages; subsequent runs should be much faster.

**Can I still pass `--conda-frontend mamba`?**  
Snakemake has deprecated alternative frontends. Rely on the **new `conda`** you installed in `snk`. That’s exactly why we included `conda` in the driver env.

**How do I update later?**
```bash
micromamba activate snk
micromamba update -n snk snakemake conda
```

---

## Optional: Per-project Snakemake profile

Put this in `.snakemake/config.yaml` under your project to keep everything local:

```yaml
# project/.snakemake/config.yaml
conda-prefix: .snakemake/conda
latency-wait: 120
printshellcmds: true
rerun-incomplete: true
```

Run:
```bash
micromamba activate snk
snakemake -c 22 --use-conda
```

---

## Clean up & disk hygiene

- Remove old cached envs:
  ```bash
  rm -rf ~/.snakemake/conda/*  # or your project’s .snakemake/conda
  ```
- Clear micromamba package cache:
  ```bash
  micromamba clean -a -y
  ```

---

## Why not just update system conda?

On shared clusters you typically don’t have permission, and changing central tooling can break other users’ workflows. This user-space approach is **safe, reversible, and reproducible**.
