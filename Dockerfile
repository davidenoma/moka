# syntax=docker/dockerfile:1.7
FROM mambaorg/micromamba:1.5.3

# --- build-time knobs ---
ARG REPO_URL="https://github.com/davidenoma/moka.git"
ARG REPO_REF="main"
# Pick the current Linux x86_64 zip from https://www.cog-genomics.org/plink/1.9/
ARG PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
ARG ENV_FILE="docker_environment.yml"

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH="/usr/local/bin:$PATH"
WORKDIR /opt

# base runtime (snakemake, mamba, wget, unzip, etc.)
COPY ${ENV_FILE} /tmp/base_env.yml
RUN micromamba install -y -n base -f /tmp/base_env.yml -c conda-forge -c bioconda \
 && micromamba clean -a -y

# clone your repo into the image (public HTTPS; use SSH + --ssh default if private)
RUN micromamba run -n base git clone --depth 1 --branch "${REPO_REF}" "${REPO_URL}" /opt/moka
WORKDIR /opt/moka

# ---- PLINK 1.9 (wget + unzip) ----
RUN micromamba run -n base bash -lc "\
    echo 'Downloading PLINK from ${PLINK_URL}'; \
    wget -q -O /tmp/plink.zip ${PLINK_URL} && \
    unzip -j /tmp/plink.zip -d /usr/local/bin && \
    chmod +x /usr/local/bin/plink && \
    /usr/local/bin/plink --version || true && \
    rm -f /tmp/plink.zip"

# ---- OPTIONAL BUT RECOMMENDED: pre-install all rule envs so no --use-conda ----
# 1) conda sections from every workflow/envs/*.yml
RUN set -eux; \
    for f in $(find workflow/envs -maxdepth 1 -type f \( -name "*.yml" -o -name "*.yaml" \)); do \
      echo 'Installing conda deps from' "$f"; \
      micromamba install -y -n base -f "$f" -c conda-forge -c bioconda || true; \
    done; \
    micromamba clean -a -y

# 2) any pip: entries inside those YAMLs
RUN micromamba run -n base python - <<'PY'\n\
import yaml,glob\n\
pkgs=set()\n\
for p in glob.glob('workflow/envs/*.y*ml'):\n\
    try:\n\
        y=yaml.safe_load(open(p)) or {}\n\
        for d in y.get('dependencies',[]):\n\
            if isinstance(d, dict) and 'pip' in d:\n\
                pkgs.update(d['pip'])\n\
    except Exception as e:\n\
        print('WARN:', p, e)\n\
open('/tmp/pip.txt','w').write('\\n'.join(sorted(pkgs))+'\\n')\n\
PY
RUN if [ -s /tmp/pip.txt ]; then \
      echo 'pip installing extras...'; \
      micromamba run -n base pip install -r /tmp/pip.txt; \
    fi

# default entrypoint: snakemake (no --use-conda needed)
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "snakemake"]
CMD ["--help"]
