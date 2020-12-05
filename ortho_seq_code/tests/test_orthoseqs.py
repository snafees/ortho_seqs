import numpy as np
import os

from click.testing import CliRunner

from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
from ortho_seq_code.cli import cli


def output_file_templates(basename):
    
    return [
        f"{basename}_mean.npy",
        f"{basename}_P.npy",
        f"{basename}_var.npy",
        f"{basename}_cov.npy",
        f"{basename}_reg11.npy",
        f"{basename}_Pa.npy",
        f"{basename}_P1i1.npy",
        f"{basename}_varP1i1.npy",
        f"{basename}_varP1i1.npy",
        f"{basename}_cov11i1.npy",
        f"{basename}_reg11i1.npy",
        f"{basename}_Pa1i1.npy",
        f"{basename}_P1D.npy",
        f"{basename}_varP1D.npy",
    ]


def test_cli(
        protein_seqs_no_padding,
        protein_pheno_no_padding
):
    molecule = 'protein'
    sites = 6
    dm = 20
    pop_size = 6
    poly_order = 'first'
    out_dir = '/tmp'

    runner = CliRunner()

    result = runner.invoke(
        cli, [
            protein_seqs_no_padding,
            "--pheno_file", protein_pheno_no_padding,
            "--molecule",  molecule,
            "--sites", sites,
            "--dm", dm,
            "--pop_size", pop_size,
            "--poly_order", poly_order,
            "--out_dir", out_dir,
        ]
    )

    assert result.exit_code == 0

    
def test_nucleotide_first_order(
        nucleotide_two_sites_data_dir,
        nucleotide_expected_output_dir,
        nucleotide_params_first_order
):

    orthogonal_polynomial(
        *nucleotide_params_first_order
    )

    basename = os.path.basename(nucleotide_params_first_order.seqs_filename)

    for t in output_file_templates(basename):
        expected = np.load(os.path.join(nucleotide_expected_output_dir, t))
        actual = np.load(os.path.join(nucleotide_params_first_order.out_dir, t))
        assert expected.all() == actual.all()


def test_protein_first_order(
        protein_data_dir,
        protein_nopad_expected_output_dir,
        protein_params_first_order
):

    orthogonal_polynomial(*protein_params_first_order)

    basename = os.path.basename(protein_params_first_order.seqs_filename)

    for t in output_file_templates(basename):
        expected = np.load(os.path.join(protein_nopad_expected_output_dir, t))
        actual = np.load(os.path.join(protein_params_first_order.out_dir, t))
        assert expected.all() == actual.all()
