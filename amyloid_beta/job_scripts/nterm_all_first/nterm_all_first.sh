#!/bin/bash
#SBATCH --job-name=nterm_all_first_job
#SBATCH --nodes=1

module load anaconda

conda activate ortho_seq

ortho_seq orthogonal-polynomial /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/data/amyloid_nterm_nostop.csv --molecule protein --poly_order first --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/nterm_all_first --pheno_name nscore
