#!/bin/bash

#SBATCH --mem-per-cpu=3G
#SBATCH -n 30
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=chasem@wustl.edu

module load mpich/3.3-python-2.7.15
module load r

# note: this is currently correct, don't change it. Eventually it will be moved
# into a 'group' directory or singularity container
libpath=/scratch/mblab/chasem/rnaseq_pipeline/experiments/zev_induction/r_library_for_mpi_deseq
export R_LIBS=$libpath

# on HTCF, this shouldn't need to be changed
LD_PRELOAD=/usr/lib/openmpi/lib/libmpi.so

# add #SBATCH --array=1-4%4, or however many files you want to run, to the top,
# and uncomment this line. change the path to the correct lookup path.
#read dds_input < <(sed -n ${SLURM_ARRAY_TASK_ID}p lookup.txt)

dds_input=input.rds

# for now, you'll need to copy the deseq_de.R into your path, or to your local
# directory and make it executable. Eventually this will be in some sort of
# package
mpiexec -usize $SLURM_NTASKS -np 1 Rscript ./deseq_de.R -d $dds_input -l $libpath
