#!/bin/bash
#SBATCH --job-name=ngk_hydrophobicity7_second_job
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64G

module load anaconda

conda activate ortho_seq

export PYTHONUNBUFFERED=TRUE

ortho_seq orthogonal-polynomial /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/data/amyloid_ngk_nostop.csv --molecule protein --poly_order second --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/ngk_hydrophobicity7_second --alphbt_input LIFWVM,CYA,TEGSQ --pheno_name nscore
