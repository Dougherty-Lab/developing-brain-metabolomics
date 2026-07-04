#!/bin/bash
#SBATCH --job-name=filter_recluster
#SBATCH --mem=400GB
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --output=filter_recluster_%j.log
#SBATCH --error=filter_recluster_%j.err

# configure spack
. /ref/jdlab/software/spack/share/spack/setup-env.sh

# load singularity
spack load singularityce

# run filter_recluster.R inside container
singularity exec \
   -B /scratch/jdlab/sneha/developing-brain-metabolomics \
   developing-brain-metabolomics_1.0.4.sif \
   Rscript /scratch/jdlab/sneha/developing-brain-metabolomics/src/gene-hormone/filter_recluster.R
