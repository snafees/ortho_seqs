
from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
import numpy as np
import os

def test_nucleotide_first_order(test_data_dir, nucleotide_expected_outputs_dir):
    seqs_filename = os.path.join(test_data_dir, 'test_seqs_2sites_dna.txt')
    pheno_filename = os.path.join(test_data_dir, 'trait_test_seqs_2sites_dna.txt')
    molecule = 'DNA'
    sites = 2
    dm = 4
    pop_size = 10
    poly_order = 'first'
    out_dir = '/tmp'

    orthogonal_polynomial(seqs_filename, pheno_filename, molecule, sites, dm, pop_size, poly_order, out_dir)

    basename = os.path.basename(seqs_filename)

    expected_mean = np.load(nucleotide_expected_outputs_dir, f"{basename}_mean.npy")
    actual_mean = np.load(f"{out_dir}/{basename}_mean.npy")
    assert expected_mean == actual_mean

    expected_P = np.load(nucleotide_expected_outputs_dir, f"{basename}_P.npy")
    actual_P = np.load(f"{out_dir}/{basename}_P.npy")
    assert expected_P == actual_P

    expected_var = np.load(nucleotide_expected_outputs_dir, f"{basename}_var.npy")
    actual_var = np.load(f"{out_dir}/{basename}_var.npy")
    assert expected_var == actual_var

    expected_cov = np.load(nucleotide_expected_outputs_dir, f"{basename}_cov.npy")
    actual_cov = np.load(f"{out_dir}/{basename}_cov.npy")
    assert expected_cov == actual_cov

    expected_reg11 = np.load(nucleotide_expected_outputs_dir, f"{basename}_reg11.npy")
    actual_reg11 = np.load(f"{out_dir}/{basename}_reg11.npy")
    assert expected_reg11 == actual_reg11

    expected_Pa = np.load(nucleotide_expected_outputs_dir, f"{basename}_Pa.npy")
    actual_Pa = np.load(f"{out_dir}/{basename}_Pa.npy")
    assert expected_Pa == actual_Pa

    expected_P1i1 = np.load(nucleotide_expected_outputs_dir, f"{basename}_P1i1.npy")
    actual_P1i1 = np.load(f"{out_dir}/{basename}_P1i1.npy")
    assert expected_P1i1 == actual_P1i1

    expected_varP1i1 = np.load(nucleotide_expected_outputs_dir, f"{basename}_varP1i1.npy")
    actual_varP1i1 = np.load(f"{out_dir}/{basename}_varP1i1.npy")
    assert expected_varP1i1 == actual_varP1i1

    expected_cov11i1 = np.load(nucleotide_expected_outputs_dir, f"{basename}_cov11i1.npy")
    actual_cov11i1 = np.load(f"{out_dir}/{basename}_cov11i1.npy")
    assert expected_cov11i1 == actual_cov11i1

    expected_reg11i1 = np.load(nucleotide_expected_outputs_dir, f"{basename}_reg11i1.npy")
    actual_reg11i1 = np.load(f"{out_dir}/{basename}_reg11i1.npy")
    assert expected_reg11i1 == actual_reg11i1

    expected_Pa1i1 = np.load(nucleotide_expected_outputs_dir, f"{basename}_Pa1i1.npy")
    actual_Pa1i1 = np.load(f"{out_dir}/{basename}_Pa1i1.npy")
    assert expected_Pa1i1 == actual_Pa1i1

    expected_P1D = np.load(nucleotide_expected_outputs_dir, f"{basename}_P1D.npy")
    actual_P1D = np.load(f"{out_dir}/{basename}_P1D.npy")
    assert expected_P1D == actual_P1D

    expected_varP1D = np.load(nucleotide_expected_outputs_dir, f"{basename}_varP1D.npy")
    actual_varP1D = np.load(f"{out_dir}/{basename}_varP1D.npy")
    assert expected_varP1D == actual_varP1D

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

    expected_P = np.load(protein_expected_outputs_dir, f"{basename}_P.npy")
    actual_P = np.load(f"{out_dir}/{basename}_P.npy")
    assert expected_P == actual_P

    expected_var = np.load(protein_expected_outputs_dir, f"{basename}_var.npy")
    actual_var = np.load(f"{out_dir}/{basename}_var.npy")
    assert expected_var == actual_var

    expected_cov = np.load(protein_expected_outputs_dir, f"{basename}_cov.npy")
    actual_cov = np.load(f"{out_dir}/{basename}_cov.npy")
    assert expected_cov == actual_cov

    expected_reg11 = np.load(protein_expected_outputs_dir, f"{basename}_reg11.npy")
    actual_reg11 = np.load(f"{out_dir}/{basename}_reg11.npy")
    assert expected_reg11 == actual_reg11

    expected_Pa = np.load(protein_expected_outputs_dir, f"{basename}_Pa.npy")
    actual_Pa = np.load(f"{out_dir}/{basename}_Pa.npy")
    assert expected_Pa == actual_Pa

    expected_P1i1 = np.load(protein_expected_outputs_dir, f"{basename}_P1i1.npy")
    actual_P1i1 = np.load(f"{out_dir}/{basename}_P1i1.npy")
    assert expected_P1i1 == actual_P1i1

    expected_varP1i1 = np.load(protein_expected_outputs_dir, f"{basename}_varP1i1.npy")
    actual_varP1i1 = np.load(f"{out_dir}/{basename}_varP1i1.npy")
    assert expected_varP1i1 == actual_varP1i1

    expected_cov11i1 = np.load(protein_expected_outputs_dir, f"{basename}_cov11i1.npy")
    actual_cov11i1 = np.load(f"{out_dir}/{basename}_cov11i1.npy")
    assert expected_cov11i1 == actual_cov11i1

    expected_reg11i1 = np.load(protein_expected_outputs_dir, f"{basename}_reg11i1.npy")
    actual_reg11i1 = np.load(f"{out_dir}/{basename}_reg11i1.npy")
    assert expected_reg11i1 == actual_reg11i1

    expected_Pa1i1 = np.load(protein_expected_outputs_dir, f"{basename}_Pa1i1.npy")
    actual_Pa1i1 = np.load(f"{out_dir}/{basename}_Pa1i1.npy")
    assert expected_Pa1i1 == actual_Pa1i1

    expected_P1D = np.load(protein_expected_outputs_dir, f"{basename}_P1D.npy")
    actual_P1D = np.load(f"{out_dir}/{basename}_P1D.npy")
    assert expected_P1D == actual_P1D

    expected_varP1D = np.load(protein_expected_outputs_dir, f"{basename}_varP1D.npy")
    actual_varP1D = np.load(f"{out_dir}/{basename}_varP1D.npy")
    assert expected_varP1D == actual_varP1D
