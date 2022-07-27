#!/bin/bash
#SBATCH --job-name=gk_all_second_job
#SBATCH --nodes=1
#SBATCH -o --gkall2_%j.out
#SBATCH -e --gkall2_%j.err

module load anaconda

conda activate ortho_seq

ortho_seq orthogonal-polynomial /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/data/amyloid_gk_nostop.csv --molecule protein --poly_order second --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/gk_all_second --pheno_name nscore

