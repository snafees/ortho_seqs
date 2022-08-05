#!/bin/bash
#SBATCH --job-name=cgk_polarity_second_job
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16G

module load anaconda

conda activate ortho_seq

export PYTHONUNBUFFERED=TRUE

ortho_seq orthogonal-polynomial /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/data/amyloid_cgk_nostop.csv --molecule protein --poly_order second --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/cgk_polarity_second --alphbt_input protein_pnp --pheno_name nscore
