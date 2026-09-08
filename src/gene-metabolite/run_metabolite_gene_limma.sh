#!/bin/bash
# run_metabolite_gene_limma.sh
# ---------------------------------------------------------------------------
# Submit metabolite-gene limma-voom loop to HTCF via Singularity.
# Run from any directory:
#   sbatch run_metabolite_gene_limma.sh
#
# METHOD is set inside metabolite_gene_limma.R.
# To run both methods, change METHOD in the R script between submissions --
# outputs go to separate parquet/ vs parquet-<method>/ dirs 
# ---------------------------------------------------------------------------

#SBATCH --job-name=metab_gene_limma
#SBATCH --output=logs/metab_gene_limma_%j.out
#SBATCH --error=logs/metab_gene_limma_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL --mail-user=sneha.chaturvedi@wustl.edu

# ---- paths -------------------------------------------------------------------
REPO=/scratch/jdlab/sneha/developing-brain-metabolomics
SIF=${REPO}/bin/docker/developing-brain-metabolomics_1.0.4.sif
SCRIPT_DIR=${REPO}/src/gene-metabolite

# ---- log dir (relative to SCRIPT_DIR, where output paths resolve) -----------
mkdir -p ${SCRIPT_DIR}/logs

echo "Job ID     : ${SLURM_JOB_ID}"
echo "Start time : $(date)"
echo "Node       : $(hostname)"
echo "SIF        : ${SIF}"
echo "Script dir : ${SCRIPT_DIR}"
echo "---"

# ---- run ---------------------------------------------------------------------
cd ${SCRIPT_DIR}

singularity exec \
  --bind ${REPO}:${REPO} \
  ${SIF} \
  Rscript metabolite_gene_limma.R

EXIT_CODE=$?
echo "---"
echo "End time  : $(date)"
echo "Exit code : ${EXIT_CODE}"
exit ${EXIT_CODE}
