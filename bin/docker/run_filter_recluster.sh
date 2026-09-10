#!/bin/bash
#SBATCH --job-name=filter_recluster
#SBATCH --mem=400GB
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=filter_recluster_%j.log
#SBATCH --error=filter_recluster_%j.err

# configure spack
. /ref/jdlab/software/spack/share/spack/setup-env.sh

# load singularity
spack load singularityce@3.8.0%gcc@8.5.0

# run filter_recluster.R inside container
singularity exec \
   -B /scratch/jdlab/emma/developing-brain-metabolomics \
   developing-brain-metabolomics_1.0.5.sif \
   Rscript /scratch/jdlab/emma/developing-brain-metabolomics/src/gene-hormone/filter_recluster.R
