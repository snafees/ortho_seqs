#!/bin/bash
#SBATCH --job-name=data_process
#SBATCH --nodes=1
#SBATCH -o data_proc_output_%j.out
#SBATCH -e data_proc_errors_%j.err

module load anaconda

conda activate ortho_seq

python /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/code/data_process.py
