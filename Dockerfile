# syntax=docker/dockerfile:1.7
FROM mambaorg/micromamba:1.5.3
SHELL ["/bin/bash", "-lc"]

# ---- build-time args ----
ARG REPO_URL="https://github.com/davidenoma/moka.git"
ARG REPO_REF="docker_config"
ARG PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
ARG ENV_FILE="docker_environment.yml"

# micromamba root + PATH (include user's ~/.local/bin for plink)
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH="/opt/conda/bin:/usr/local/bin:/home/mambauser/.local/bin:${PATH}"

# ---- base runtime (snakemake, mamba, etc.) ----
COPY ${ENV_FILE} /tmp/base_env.yml
RUN micromamba install -y -n base -f /tmp/base_env.yml -c conda-forge -c bioconda \
 && micromamba clean -a -y
# --- 2) Add a real 'conda' + 'mamba' into the SAME base env and expose on PATH
# ---- base runtime (snakemake, etc.) ----
COPY ${ENV_FILE} /tmp/base_env.yml
RUN micromamba install -y -n base -f /tmp/base_env.yml -c conda-forge -c bioconda \
 && micromamba install -y -n base conda mamba -c conda-forge \
 && micromamba clean -a -y \
 && micromamba run -n base conda info --json >/dev/null

# Ensure we actually have wget+unzip available in the build
RUN micromamba install -y -n base wget unzip -c conda-forge \
 && micromamba clean -a -y

# ---- working directory: /home/mambauser/moka (no chown, no root) ----
WORKDIR /home/mambauser
WORKDIR ./moka

# ---- clone your repo into the current dir (.) ----
RUN micromamba run -n base git clone --depth 1 --branch "${REPO_REF}" "${REPO_URL}" .

# ---- PLINK 1.9 (download into /opt/conda/bin so PATH sees it) ----
RUN micromamba run -n base wget -q -O /tmp/plink.zip "${PLINK_URL}" \
 && micromamba run -n base unzip -j /tmp/plink.zip -d /opt/conda/bin \
 && chmod +x /opt/conda/bin/plink \
 && micromamba run -n base plink --version \
 && rm -f /tmp/plink.zip



# ---- OPTIONAL: pre-install rule envs into the same env (no --use-conda) ----
# 1) conda sections from workflow/envs/*.yml
RUN set -eux; \
    if compgen -G 'workflow/envs/*.y*ml' > /dev/null; then \
      for f in workflow/envs/*.y*ml; do \
        echo "Installing conda deps from $f"; \
        micromamba install -y -n base -f "$f" -c conda-forge -c bioconda || true; \
      done; \
      micromamba clean -a -y; \
    else \
      echo 'No workflow/envs/*.y*ml found; skipping.'; \
    fi

# 2) any pip: entries inside those YAMLs
RUN micromamba run -n base python - <<'PY'
import glob, sys
try:
    import yaml
except Exception:
    print("PyYAML not found; skipping pip extraction.", file=sys.stderr); sys.exit(0)
pkgs=set()
for p in glob.glob('workflow/envs/*.y*ml'):
    try:
        with open(p) as fh:
            y=yaml.safe_load(fh) or {}
        for d in y.get('dependencies', []):
            if isinstance(d, dict) and 'pip' in d:
                pkgs.update(d['pip'])
    except Exception as e:
        print('WARN:', p, e, file=sys.stderr)
open('/tmp/pip.txt','w').write('\n'.join(sorted(pkgs))+'\n')
print('pip packages:', len(pkgs))
PY
RUN if [ -s /tmp/pip.txt ]; then \
      micromamba run -n base pip install -r /tmp/pip.txt; \
    fi

# ---- default: run Snakemake with preloaded env ----
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "snakemake"]
CMD ["--help"]
