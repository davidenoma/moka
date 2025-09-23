# syntax=docker/dockerfile:1.7
FROM continuumio/miniconda3:latest
SHELL ["/bin/bash", "-lc"]

# --- knobs ---
ARG REPO_URL="https://github.com/davidenoma/moka.git"
ARG REPO_REF="main"
ARG ENV_FILE="docker_environment.yml"

# Use only conda-forge + bioconda; remove defaults to avoid TOS
RUN conda update -n base -y conda \
 && conda config --system --remove channels defaults || true \
 && conda config --system --add channels conda-forge \
 && conda config --system --add channels bioconda \
 && conda config --system --set channel_priority strict \
 && conda install -n base -y conda-libmamba-solver \
 && conda config --system --set solver libmamba

# ---- base runtime; your YAMLs define most deps ----
COPY ${ENV_FILE} /tmp/base_env.yml
# strip any explicit 'defaults' entries in your YAML (belt & suspenders)
RUN sed -i '/^\s*-\s*defaults\s*$/d' /tmp/base_env.yml || true
RUN if [ -s /tmp/base_env.yml ]; then \
      echo ">> updating base from ${ENV_FILE}"; \
      conda env update -n base -f /tmp/base_env.yml || true; \
    fi \
 && conda install -n base -y snakemake git wget unzip parallel pyyaml \
 && conda clean -a -y

# ---- repo ----
WORKDIR /home/conda
RUN mkdir -p moka
WORKDIR /home/conda/moka
RUN git clone --depth 1 --branch "${REPO_REF}" "${REPO_URL}" .

# ---- merge ALL rule envs into base (so no --use-conda at runtime) ----
RUN set -eux; \
    shopt -s nullglob; \
    for f in workflow/envs/*.yml workflow/envs/*.yaml; do \
      sed -i '/^\s*-\s*defaults\s*$/d' "$f" || true; \
      echo ">> merging $f into base"; \
      conda env update -n base -f "$f" || true; \
    done; \
    conda clean -a -y

# ---- install any pip: extras from those YAMLs into base ----
RUN python - <<'PY'
import glob, sys, yaml, subprocess, tempfile
pip_pkgs=set()
for p in glob.glob('workflow/envs/*.y*ml'):
    try:
        y=yaml.safe_load(open(p)) or {}
        for d in y.get('dependencies', []):
            if isinstance(d, dict) and 'pip' in d:
                pip_pkgs.update(d['pip'] or [])
    except Exception as e:
        print("WARN:", p, e, file=sys.stderr)
if pip_pkgs:
    with tempfile.NamedTemporaryFile('w', delete=False) as fh:
        fh.write('\n'.join(sorted(pip_pkgs))+'\n'); name=fh.name
    subprocess.check_call(["bash","-lc",f"conda run -n base pip install -r {name}"])
else:
    print("No pip extras found.")
PY

# ---- default entrypoint: run snakemake in base env ----
RUN printf '%s\n' '#!/usr/bin/env bash' \
                  'set -euo pipefail' \
                  'exec conda run -n base snakemake "$@"' \
   > /usr/local/bin/snakemake-entry.sh \
 && chmod +x /usr/local/bin/snakemake-entry.sh

ENTRYPOINT ["/usr/local/bin/snakemake-entry.sh"]
CMD ["--help"]
