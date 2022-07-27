#!/bin/bash
#SBATCH --job-name=gk_alberts_second_job
#SBATCH --nodes=1

module load anaconda

conda activate ortho_seq

ortho_seq orthogonal-polynomial /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/data/amyloid_gk_nostop.csv --molecule protein --poly_order second --alphbt_input alberts --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/gk_alberts_second --pheno_name nscore

