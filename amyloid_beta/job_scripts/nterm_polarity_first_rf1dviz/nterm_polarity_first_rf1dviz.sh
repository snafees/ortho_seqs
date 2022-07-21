#!/bin/bash
#SBATCH --job-name=nterm_polarity_first_rf1dviz_job
#SBATCH --nodes=1

module load anaconda

conda activate ortho_seq

ortho_seq rf1d-viz /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/nterm_polarity_first/amyloid_nterm_nostop_regressions.npz --alphbt_input 0,z --molecule protein --phenotype nscore --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/nterm_polarity_first_rf1dviz --action ALL

