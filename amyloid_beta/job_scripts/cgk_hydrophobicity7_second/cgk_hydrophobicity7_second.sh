#!/bin/bash
#SBATCH --job-name=cgk_hydrophobicity7_second_job
#SBATCH --nodes=1

module load anaconda

conda activate ortho_seq

ortho_seq orthogonal-polynomial /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/data/amyloid_cgk_nostop.csv --molecule protein --poly_order second --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/cgk_hydrophobicity7_second --alphbt_input LIFWVM,CYA,TEGSQ --pheno_name nscore
