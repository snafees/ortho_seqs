#!/bin/bash
#SBATCH --job-name=ngk_polarity_first_rf1dviz_job
#SBATCH --nodes=1

module load anaconda

conda activate ortho_seq

ortho_seq rf1d-viz /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/ngk_polarity_first/amyloid_full_nterm_regressions.npz --molecule protein --alphbt_input polar,nonpolar --out_dir /hpc/projects/data_lg/olivia.yoo/ortho_seqs/amyloid_beta/server_ortho_seq_results/ngk_polarity_first_rf1dviz --phenotype nscore --action ALL

