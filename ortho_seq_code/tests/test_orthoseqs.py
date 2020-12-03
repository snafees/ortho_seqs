
from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
import numpy as np
import os

def test_nucleotide_first_order(test_data_dir, nucleotide_expected_outputs_dir):
    seqs_filename = os.path.join(test_data_dir, )
     #def orthogonal_polynomial(filename, pheno_file, molecule, sites, dm, pop_size, poly_order, precomputed, out_dir)


def test_protein_first_order(test_data_dir, protein_expected_outputs_dir):
    seqs_filename = os.path.join(test_data_dir, 'protein_seqs_nopad.txt')
    pheno_filename = os.path.join(test_data_dir, 'protein_pheno_nopad.txt')
    molecule = 'protein'
    sites = 6
    dm = 20
    pop_size = 6
    poly_order = 'first'
    out_dir = '/tmp'

    orthogonal_polynomial(seqs_filename, pheno_filename, molecule, sites, dm, pop_size, poly_order, out_dir)

    basename = os.path.basename(seqs_filename)
    expected_mean = np.load(protein_expected_outputs_dir, f"{basename}_mean.npy")
    actual_mean = np.load(f"{out_dir}/{basename}_mean.npy")
    assert expected_mean == actual_mean
