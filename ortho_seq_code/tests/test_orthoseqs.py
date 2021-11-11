import numpy as np
from ortho_seq_code.orthogonal_polynomial import ortho_poly_command
import pandas as pd
import os
from matplotlib import pyplot as plt

from click.testing import CliRunner

from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
from ortho_seq_code.cli import cli
from ortho_seq_code.tests import orthoseqs_tst_utils as utils
from ortho_seq_code.utils import get_seq_info


def test_cli(protein_seqs_no_padding, protein_pheno_no_padding):
    molecule = "protein"
    poly_order = "first"
    out_dir = "/tmp"
    alphbt_input = None
    min_pct = 75

    runner = CliRunner()

    result = runner.invoke(
        ortho_poly_command,
        [
            protein_seqs_no_padding,
            "--pheno_file",
            protein_pheno_no_padding,
            "--molecule",
            molecule,
            "--poly_order",
            poly_order,
            "--out_dir",
            out_dir,
            "--alphbt_input",
            alphbt_input,
            "--min_pct",
            min_pct,
        ],
    )

    assert result.exit_code == 0


def test_cli(protein_seqs_padding, protein_pheno_padding):
    molecule = "protein"
    poly_order = "first"
    out_dir = "/tmp"
    alphbt_input = None
    min_pct = 75

    runner = CliRunner()

    result = runner.invoke(
        ortho_poly_command.main,
        [
            protein_seqs_padding,
            "--pheno_file",
            protein_pheno_padding,
            "--molecule",
            molecule,
            "--poly_order",
            poly_order,
            "--out_dir",
            out_dir,
            "--alphbt_input",
            alphbt_input,
            "--min_pct",
            min_pct,
        ],
    )

    assert result.stderr == 0


def test_cli_precomputed(
    protein_seqs_no_padding, protein_pheno_no_padding, protein_data_dir
):
    molecule = "protein"
    poly_order = "first"
    out_dir = protein_data_dir
    alphbt_input = None
    min_pct = 75

    runner = CliRunner()

    result = runner.invoke(
        ortho_poly_command,
        [
            protein_seqs_no_padding,
            "--pheno_file",
            protein_pheno_no_padding,
            "--molecule",
            molecule,
            "--poly_order",
            poly_order,
            "--out_dir",
            out_dir,
            "--precomputed",
            alphbt_input,
            "--alphbt_input",
            min_pct,
            "--min_pct",
        ],
    )

    assert result.stderr == 0


def assert_equality(expected_path, actual_path):
    assert os.path.exists(expected_path)
    assert os.path.exists(actual_path)
    obtained_arrays = np.load(actual_path)
    expected_arrays = np.load(expected_path)
    for key, obtained_array in obtained_arrays.items():
        expected_array = expected_arrays[key]
        np.testing.assert_array_equal(
            expected_array, obtained_array, "error at {}".format(key)
        )


def test_nucleotide_first_order(
    nucleotide_first_order_data_dir, nucleotide_params_first_order
):

    with utils.TempDirectory() as location:
        nucleotide_params_first_order = nucleotide_params_first_order._replace(
            out_dir=location
        )
        location += "0"  # Test file paths already exists, but location isn't updated with out_dir, this will update it
        orthogonal_polynomial(*nucleotide_params_first_order)
        basename = os.path.basename(nucleotide_params_first_order.seqs_filename)
        basename_pheno = os.path.basename(nucleotide_params_first_order.pheno_filename)
        expected_path = os.path.join(nucleotide_first_order_data_dir, basename + ".npz")
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            nucleotide_first_order_data_dir, basename_pheno + "_regressions.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_regressions.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            nucleotide_first_order_data_dir, basename_pheno + "_covs_with_F.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_covs_with_F.npz")
        expected_path = np.load(
            os.path.join(nucleotide_first_order_data_dir, basename_pheno + "_Fm.npy")
        )
        obtained_path = np.load(os.path.join(location, basename_pheno + "_Fm.npy"))
        np.testing.assert_array_equal(expected_path, obtained_path)


def test_nucleotide_second_order(
    nucleotide_second_order_data_dir, nucleotide_params_second_order
):
    with utils.TempDirectory() as location:
        nucleotide_params_second_order = nucleotide_params_second_order._replace(
            out_dir=location
        )
        location += "0"
        orthogonal_polynomial(*nucleotide_params_second_order)

        basename = os.path.basename(nucleotide_params_second_order.seqs_filename)
        basename_pheno = os.path.basename(nucleotide_params_second_order.pheno_filename)
        expected_path = os.path.join(
            nucleotide_second_order_data_dir, basename + ".npz"
        )
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            nucleotide_second_order_data_dir, basename_pheno + "_regressions.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_regressions.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            nucleotide_second_order_data_dir, basename_pheno + "_covs_with_F.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_covs_with_F.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = np.load(
            os.path.join(nucleotide_second_order_data_dir, basename_pheno + "_Fm.npy")
        )
        obtained_path = np.load(os.path.join(location, basename_pheno + "_Fm.npy"))
        np.testing.assert_array_equal(expected_path, obtained_path)


def test_protein_first_order(protein_data_dir, protein_params_first_order):

    with utils.TempDirectory() as location:
        protein_params_first_order = protein_params_first_order._replace(
            out_dir=location
        )
        location += "0"
        orthogonal_polynomial(*protein_params_first_order)

        basefile = os.path.abspath(protein_params_first_order.seqs_filename)
        assert get_seq_info(basefile, None, None)[:-4] == [18, 6, 6]

        basename = os.path.basename(protein_params_first_order.seqs_filename)

        basename_pheno = os.path.basename(protein_params_first_order.pheno_filename)
        expected_path = os.path.join(protein_data_dir, basename + ".npz")
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            protein_data_dir, basename_pheno + "_regressions.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_regressions.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            protein_data_dir, basename_pheno + "_covs_with_F.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_covs_with_F.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = np.load(
            os.path.join(protein_data_dir, basename_pheno + "_Fm.npy")
        )
        obtained_path = np.load(os.path.join(location, basename_pheno + "_Fm.npy"))
        np.testing.assert_array_equal(expected_path, obtained_path)


def test_protein_padded_first_order(
    protein_data_dir, protein_params_first_order_padded
):

    with utils.TempDirectory() as location:
        protein_params_first_order_padded = protein_params_first_order_padded._replace(
            out_dir=location
        )
        location += "0"
        orthogonal_polynomial(*protein_params_first_order_padded)

        basefile = os.path.abspath(protein_params_first_order_padded.seqs_filename)
        assert get_seq_info(basefile, None, None)[:-4] == [21, 6, 10]

        basename = os.path.basename(protein_params_first_order_padded.seqs_filename)
        basename_pheno = os.path.basename(
            protein_params_first_order_padded.pheno_filename
        )
        expected_path = os.path.join(protein_data_dir, basename + ".npz")
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            protein_data_dir, basename_pheno + "_regressions.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_regressions.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            protein_data_dir, basename_pheno + "_covs_with_F.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_covs_with_F.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = np.load(
            os.path.join(protein_data_dir, basename_pheno + "_Fm.npy")
        )
        obtained_path = np.load(os.path.join(location, basename_pheno + "_Fm.npy"))
        np.testing.assert_array_equal(expected_path, obtained_path)


def test_protein_paddded_custom_aa(protein_data_dir, protein_params_custom_aa):

    with utils.TempDirectory() as location:
        protein_params_custom_aa = protein_params_custom_aa._replace(out_dir=location)
        location += "0"
        orthogonal_polynomial(*protein_params_custom_aa)

        basefile = os.path.abspath(protein_params_custom_aa.seqs_filename)

        indices = [0, 1, 2, 5, 6]
        assert [get_seq_info(basefile, "ARSY", "protein")[x] for x in indices] == [
            6,
            6,
            10,
            ["A", "R", "S", "Y", "z", "n"],
            ["A", "R", "S", "Y", "z", "n"],
        ]

        basename = os.path.basename(protein_params_custom_aa.seqs_filename)
        basename_pheno = os.path.basename(protein_params_custom_aa.pheno_filename)
        expected_path = os.path.join(protein_data_dir, basename + ".npz")
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            protein_data_dir, basename_pheno + "_regressions.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_regressions.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            protein_data_dir, basename_pheno + "_covs_with_F.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_covs_with_F.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = np.load(
            os.path.join(protein_data_dir, basename_pheno + "_Fm.npy")
        )
        obtained_path = np.load(os.path.join(location, basename_pheno + "_Fm.npy"))
        np.testing.assert_array_equal(expected_path, obtained_path)


def test_protein_padded_custom_aa_2(protein_data_dir, protein_params_custom_aa_2):

    with utils.TempDirectory() as location:
        protein_params_custom_aa_2 = protein_params_custom_aa_2._replace(
            out_dir=location
        )
        location += "0"
        orthogonal_polynomial(*protein_params_custom_aa_2)

        basefile = os.path.abspath(protein_params_custom_aa_2.seqs_filename)

        indices = [0, 1, 2, 5, 6]
        assert [get_seq_info(basefile, "AR,SY", "protein")[x] for x in indices] == [
            4,
            6,
            10,
            ["0", "1", "2", "3"],
            ["AR", "SY", "CDEFGHIKLMNPQTVW", "n"],
        ]

        basename = os.path.basename(protein_params_custom_aa_2.seqs_filename)
        basename_pheno = os.path.basename(protein_params_custom_aa_2.pheno_filename)
        expected_path = os.path.join(protein_data_dir, basename + ".npz")
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            protein_data_dir, basename_pheno + "_regressions.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_regressions.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = os.path.join(
            protein_data_dir, basename_pheno + "_covs_with_F.npz"
        )
        obtained_path = os.path.join(location, basename_pheno + "_covs_with_F.npz")
        assert_equality(expected_path, obtained_path)
        expected_path = np.load(
            os.path.join(protein_data_dir, basename_pheno + "_Fm.npy")
        )
        obtained_path = np.load(os.path.join(location, basename_pheno + "_Fm.npy"))
        np.testing.assert_array_equal(expected_path, obtained_path)
