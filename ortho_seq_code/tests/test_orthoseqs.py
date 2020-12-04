
from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
import numpy as np
import os

def test_nucleotide_first_order(nucleotide_two_sites_data_dir, nucleotide_expected_output_dir):
    test_data_dir = nucleotide_two_sites_data_dir
    seqs_filename = os.path.join(test_data_dir, 'test_seqs_2sites_dna.txt')
    pheno_filename = os.path.join(test_data_dir, 'trait_test_seqs_2sites_dna.txt')
    molecule = 'DNA'
    sites = 2
    dm = 4
    pop_size = 12
    poly_order = 'first'
    out_dir = '/tmp'

    orthogonal_polynomial(seqs_filename, pheno_filename, molecule, sites, dm, pop_size, poly_order, False, out_dir)

    basename = os.path.basename(seqs_filename)

    expected_mean = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_mean.npy"))
    actual_mean = np.load(f"{out_dir}/{basename}_mean.npy")
    assert expected_mean.all() == actual_mean.all()

    expected_P = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_P.npy"))
    actual_P = np.load(f"{out_dir}/{basename}_P.npy")
    assert expected_P.all() == actual_P.all()

    expected_var = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_var.npy"))
    actual_var = np.load(f"{out_dir}/{basename}_var.npy")
    assert expected_var.all() == actual_var.all()

    expected_cov = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_cov.npy"))
    actual_cov = np.load(f"{out_dir}/{basename}_cov.npy")
    assert expected_cov.all() == actual_cov.all()

    expected_reg11 = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_reg11.npy"))
    actual_reg11 = np.load(f"{out_dir}/{basename}_reg11.npy")
    assert expected_reg11.all() == actual_reg11.all()

    expected_Pa = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_Pa.npy"))
    actual_Pa = np.load(f"{out_dir}/{basename}_Pa.npy")
    assert expected_Pa.all() == actual_Pa.all()

    expected_P1i1 = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_P1i1.npy"))
    actual_P1i1 = np.load(f"{out_dir}/{basename}_P1i1.npy")
    assert expected_P1i1.all() == actual_P1i1.all()

    expected_varP1i1 = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_varP1i1.npy"))
    actual_varP1i1 = np.load(f"{out_dir}/{basename}_varP1i1.npy")
    assert expected_varP1i1.all() == actual_varP1i1.all()

    expected_cov11i1 = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_cov11i1.npy"))
    actual_cov11i1 = np.load(f"{out_dir}/{basename}_cov11i1.npy")
    assert expected_cov11i1.all() == actual_cov11i1.all()

    expected_reg11i1 = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_reg11i1.npy"))
    actual_reg11i1 = np.load(f"{out_dir}/{basename}_reg11i1.npy")
    assert expected_reg11i1.all() == actual_reg11i1.all()

    expected_Pa1i1 = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_Pa1i1.npy"))
    actual_Pa1i1 = np.load(f"{out_dir}/{basename}_Pa1i1.npy")
    assert expected_Pa1i1.all() == actual_Pa1i1.all()

    expected_P1D = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_P1D.npy"))
    actual_P1D = np.load(f"{out_dir}/{basename}_P1D.npy")
    assert expected_P1D.all() == actual_P1D.all()

    expected_varP1D = np.load(os.path.join(nucleotide_expected_output_dir, f"{basename}_varP1D.npy"))
    actual_varP1D = np.load(f"{out_dir}/{basename}_varP1D.npy")
    assert expected_varP1D.all() == actual_varP1D.all()

def test_protein_first_order(protein_data_dir, protein_nopad_expected_output_dir):
    test_data_dir = protein_data_dir
    protein_expected_output_dir = protein_nopad_expected_output_dir
    seqs_filename = os.path.join(test_data_dir, 'protein_seqs_nopad.txt')
    pheno_filename = os.path.join(test_data_dir, 'protein_pheno_nopad.txt')
    molecule = 'protein'
    sites = 6
    dm = 20
    pop_size = 6
    poly_order = 'first'
    out_dir = '/tmp'

    orthogonal_polynomial(seqs_filename, pheno_filename, molecule, sites, dm, pop_size, poly_order, False, out_dir)

    basename = os.path.basename(seqs_filename)

    expected_mean = np.load(os.path.join(protein_expected_output_dir, f"{basename}_mean.npy"))
    actual_mean = np.load(f"{out_dir}/{basename}_mean.npy")
    assert expected_mean.all() == actual_mean.all()

    expected_P = np.load(os.path.join(protein_expected_output_dir, f"{basename}_P.npy"))
    actual_P = np.load(f"{out_dir}/{basename}_P.npy")
    assert expected_P.all() == actual_P.all()

    expected_var = np.load(os.path.join(protein_expected_output_dir, f"{basename}_var.npy"))
    actual_var = np.load(f"{out_dir}/{basename}_var.npy")
    assert expected_var.all() == actual_var.all()

    expected_cov = np.load(os.path.join(protein_expected_output_dir, f"{basename}_cov.npy"))
    actual_cov = np.load(f"{out_dir}/{basename}_cov.npy")
    assert expected_cov.all() == actual_cov.all()

    expected_reg11 = np.load(os.path.join(protein_expected_output_dir, f"{basename}_reg11.npy"))
    actual_reg11 = np.load(f"{out_dir}/{basename}_reg11.npy")
    assert expected_reg11.all() == actual_reg11.all()

    expected_Pa = np.load(os.path.join(protein_expected_output_dir, f"{basename}_Pa.npy"))
    actual_Pa = np.load(f"{out_dir}/{basename}_Pa.npy")
    assert expected_Pa.all() == actual_Pa.all()

    expected_P1i1 = np.load(os.path.join(protein_expected_output_dir, f"{basename}_P1i1.npy"))
    actual_P1i1 = np.load(f"{out_dir}/{basename}_P1i1.npy")
    assert expected_P1i1.all() == actual_P1i1.all()

    expected_varP1i1 = np.load(os.path.join(protein_expected_output_dir, f"{basename}_varP1i1.npy"))
    actual_varP1i1 = np.load(f"{out_dir}/{basename}_varP1i1.npy")
    assert expected_varP1i1.all() == actual_varP1i1.all()

    expected_cov11i1 = np.load(os.path.join(protein_expected_output_dir, f"{basename}_cov11i1.npy"))
    actual_cov11i1 = np.load(f"{out_dir}/{basename}_cov11i1.npy")
    assert expected_cov11i1.all() == actual_cov11i1.all()

    expected_reg11i1 = np.load(os.path.join(protein_expected_output_dir, f"{basename}_reg11i1.npy"))
    actual_reg11i1 = np.load(f"{out_dir}/{basename}_reg11i1.npy")
    assert expected_reg11i1.all() == actual_reg11i1.all()

    expected_Pa1i1 = np.load(os.path.join(protein_expected_output_dir, f"{basename}_Pa1i1.npy"))
    actual_Pa1i1 = np.load(f"{out_dir}/{basename}_Pa1i1.npy")
    assert expected_Pa1i1.all() == actual_Pa1i1.all()

    expected_P1D = np.load(os.path.join(protein_expected_output_dir, f"{basename}_P1D.npy"))
    actual_P1D = np.load(f"{out_dir}/{basename}_P1D.npy")
    assert expected_P1D.all() == actual_P1D.all()

    expected_varP1D = np.load(os.path.join(protein_expected_output_dir, f"{basename}_varP1D.npy"))
    actual_varP1D = np.load(f"{out_dir}/{basename}_varP1D.npy")
    assert expected_varP1D.all() == actual_varP1D.all()
