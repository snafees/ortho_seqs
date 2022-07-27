#!/bin/bash
#SBATCH --job-name=split_amyloid_job
#SBATCH --nodes=1

module load anaconda

conda activate ortho_seq

python /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/code/split_amyloid.py
