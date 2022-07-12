#!/bin/bash
#BSUB -J gk_all_first_job
#BSUB -n 1
#BSUB -e %J.err
#BSUB -o %J.out

module load python

conda activate ortho_seqs

cd /mnt/ibm_lg/olivia.yoo/ortho_seqs

ortho_seq orthogonal-polynomial /mnt/ibm_lg/olivia.yoo/ortho_seqs/amyloid_beta/data/amyloid_gk_nostop.csv --molecule protein --poly_order first --out_dir /mnt/ibm_lg/olivia.yoo/ortho_seqs/amyloid_beta/ortho_seq_results/gk_all_first --pheno_name nscore
