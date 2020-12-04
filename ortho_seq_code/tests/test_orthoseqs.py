
from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
from ortho_seq_code.cli import cli
import numpy as np
import os

from click.testing import CliRunner


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



def test_protein_first_order(
        protein_seqs_no_padding,
        protein_pheno_no_padding,
        protein_expected_output_dir

):
    molecule = 'protein'
    sites = 6
    dm = 20
    pop_size = 6
    poly_order = 'first'
    out_dir = '/tmp'

    orthogonal_polynomial(
        protein_seqs_no_padding,
        protein_pheno_no_padding,
        molecule,
        sites,
        dm,
        pop_size,
        poly_order,
        False,
        out_dir
    )

    basename = os.path.basename(protein_seqs_no_padding)
    expected_mean = np.load(protein_expected_output_dir, f"{basename}_mean.npy")
    actual_mean = np.load(f"{out_dir}/{basename}_mean.npy")
    assert expected_mean == actual_mean
