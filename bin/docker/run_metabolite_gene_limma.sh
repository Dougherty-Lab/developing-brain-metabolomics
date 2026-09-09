#!/bin/bash
# run_metabolite_gene_limma.sh
# ---------------------------------------------------------------------------
# Submit metabolite-gene limma-voom loop to SLURM via Singularity.
#
# Inputs:  4.1_metabolite_gene_limma.R
# Outputs: results/gene-metabolite/parquet-log2_na/*.parquet
#
# Usage: sbatch src/gene-metabolite/run_metabolite_gene_limma.sh
#        (run from the repo root, or the script auto-detects it)
# ---------------------------------------------------------------------------

#SBATCH --job-name=metab_gene_limma
#SBATCH --output=logs/metab_gene_limma_%j.out
#SBATCH --error=logs/metab_gene_limma_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL --mail-user=YOUR_EMAIL@institution.edu  # <-- set to your email

# ---- auto-detect repo root from script location -----------------------------
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_PATH}/../.." && pwd)"

# Verify we found the right place
if [ ! -d "${REPO}/.git" ]; then
  echo "ERROR: could not find repo root (no .git at ${REPO})" >&2
  exit 1
fi

SIF=${REPO}/bin/docker/developing-brain-metabolomics_1.0.5.sif
SCRIPT_DIR=${REPO}/src/gene-metabolite

# ---- log dir -----------------------------------------------------------------
mkdir -p ${SCRIPT_DIR}/logs

echo "Job ID     : ${SLURM_JOB_ID}"
echo "Start time : $(date)"
echo "Node       : $(hostname)"
echo "Repo       : ${REPO}"
echo "SIF        : ${SIF}"
echo "Script dir : ${SCRIPT_DIR}"
echo "---"

# ---- run ---------------------------------------------------------------------
cd ${SCRIPT_DIR}

singularity exec \
  --bind ${REPO}:${REPO} \
  ${SIF} \
  Rscript src/gene-metabolite/4.1_metabolite_gene_limma.R

EXIT_CODE=$?
echo "---"
echo "End time  : $(date)"
echo "Exit code : ${EXIT_CODE}"
exit ${EXIT_CODE}

# ---- AI assistance disclosure ------------------------------------------------
# Code in this script was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.
