#!/bin/bash

#SBATCH --time=15:00:00  # right now, 15 hours. change depending on time expectation to run
#SBATCH --mem-per-cpu=10G
#SBATCH -J your_jobname.out
#SBATCH -o your_jobname.out

ml miniconda

# until HTCF updates and spack is available, this works
source activate /scratch/mblab/chasem/rnaseq_pipeline/conda_envs/nextflow

mkdir tmp

nextflow run /path/to/brentlab_rnaseq_nf/main.nf -params-file /path/to/your_params.json
