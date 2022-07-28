#!/bin/bash
#SBATCH --job-name=abeta_charge_first_job
#SBATCH --nodes=1

module load anaconda

conda activate ortho_seq

ortho_seq orthogonal-polynomial /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/data/complete_amyloid1.csv --molecule protein --poly_order first --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/abeta_charge_first --alphbt_input DE,RHK --pheno_name nscore
