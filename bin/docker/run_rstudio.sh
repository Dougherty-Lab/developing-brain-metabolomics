#!/bin/bash
#SBATCH --job-name=rstudio
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00


# configure spack
. /ref/jdlab/software/spack/share/spack/setup-env.sh

# load singularity
spack load singularityce@3.8.0%gcc@8.5.0

# print ssh address
host=$(hostname)

echo -e "
        SSH tunnel to $host
        ssh $USER@login.htcf.wustl.edu -N -L 8787:$host:8787

        Open in browser: http://localhost:8787
"
# execute container
singularity exec \
   -B /tmp:/var/lib/rstudio-server \
   -B /tmp:/var/run/rstudio-server \
   -B /scratch/jdlab/sneha/developing-brain-metabolomics \
   developing-brain-metabolomics_1.0.5.sif \
   rserver \
     --server-user=$USER

# expected output: INFO:    Converting SIF file to temporary sandbox...
