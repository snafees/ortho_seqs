import pytest
from pathlib import Path
from collections import namedtuple
import os


Params = namedtuple(
    "Params",
    [
        "seqs_filename",
        "pheno_filename",
        "molecule",
        "poly_order",
        "precomputed",
        "out_dir",
        "alphbt_input",
        "min_pct",
    ],
)


@pytest.fixture
def test_data_dir():
    return os.path.join(Path(__file__).parent, "tests", "data")


@pytest.fixture
def protein_data_dir(test_data_dir):
    return os.path.join(test_data_dir, "protein")


@pytest.fixture
def protein_seqs_no_padding(protein_data_dir):
    return os.path.join(protein_data_dir, "protein_seqs_nopad.txt")


@pytest.fixture
def protein_pheno_no_padding(protein_data_dir):
    return os.path.join(protein_data_dir, "protein_pheno_nopad.txt")


@pytest.fixture
def protein_seqs_padding(protein_data_dir):
    return os.path.join(protein_data_dir, "protein_seqs_padded.txt")


@pytest.fixture
def protein_pheno_padding(protein_data_dir):
    return os.path.join(protein_data_dir, "protein_pheno_padded.txt")


@pytest.fixture
def protein_nopad_expected_output_dir(protein_expected_output_dir):
    return os.path.join(protein_expected_output_dir, "protein_seqs_nopad")


@pytest.fixture
def protein_pad_expected_output_dir(protein_expected_output_dir):
    return os.path.join(protein_expected_output_dir, "protein_seqs_padded")


@pytest.fixture
def nucleotide_data_dir(test_data_dir):
    return os.path.join(test_data_dir, "nucleotide")


@pytest.fixture
def nucleotide_first_order_data_dir(nucleotide_data_dir):
    return os.path.join(nucleotide_data_dir, "first_order")


@pytest.fixture
def nucleotide_second_order_data_dir(nucleotide_data_dir):
    return os.path.join(nucleotide_data_dir, "second_order")


@pytest.fixture
def nucleotide_params_first_order(nucleotide_first_order_data_dir):
    seqs_filename = os.path.join(
        nucleotide_first_order_data_dir, "test_seqs_2sites_dna.txt"
    )
    pheno_filename = os.path.join(
        nucleotide_first_order_data_dir, "trait_test_seqs_2sites_dna.txt"
    )

    return Params(
        seqs_filename,
        pheno_filename,
        "DNA",
        "first",
        False,
        "",
        None,
        72,
    )


@pytest.fixture
def nucleotide_params_second_order(nucleotide_second_order_data_dir):
    seqs_filename = os.path.join(
        nucleotide_second_order_data_dir, "test_seqs_2sites_dna.txt"
    )
    pheno_filename = os.path.join(
        nucleotide_second_order_data_dir, "trait_test_seqs_2sites_dna.txt"
    )

    return Params(
        seqs_filename,
        pheno_filename,
        "DNA",
        "second",
        False,
        "",
        None,
        72,
    )


@pytest.fixture
def protein_params_first_order(protein_data_dir):
    seqs_filename = os.path.join(protein_data_dir, "protein_seqs_nopad.txt")
    pheno_filename = os.path.join(protein_data_dir, "protein_pheno_nopad.txt")

    return Params(
        seqs_filename,
        pheno_filename,
        "protein",
        "first",
        False,
        "",
        None,
        72,
    )


@pytest.fixture
def protein_params_first_order_padded(protein_data_dir):
    seqs_filename = os.path.join(protein_data_dir, "protein_seqs_padded.txt")
    pheno_filename = os.path.join(protein_data_dir, "protein_pheno_padded.txt")

    return Params(
        seqs_filename,
        pheno_filename,
        "protein_n",
        "first",
        False,
        "",
        None,
        72,
    )


@pytest.fixture
def protein_params_custom_aa(protein_data_dir):
    seqs_filename = os.path.join(protein_data_dir, "protein_seqs_padded_custom_aa.txt")
    pheno_filename = os.path.join(
        protein_data_dir, "protein_pheno_padded_custom_aa.txt"
    )

    return Params(
        seqs_filename,
        pheno_filename,
        "protein",
        "first",
        False,
        "",
        "YSAR",
        72,
    )


@pytest.fixture
def protein_params_custom_aa_2(protein_data_dir):
    seqs_filename = os.path.join(
        protein_data_dir, "protein_seqs_padded_custom_aa_2.txt"
    )
    pheno_filename = os.path.join(
        protein_data_dir, "protein_pheno_padded_custom_aa_2.txt"
    )

    return Params(
        seqs_filename,
        pheno_filename,
        "protein",
        "first",
        False,
        "",
        "AR,SY",
        72,
    )


@pytest.fixture
def protein_first_order_autopad(protein_data_dir):
    seqs_filename = os.path.join(protein_data_dir, "protein_seqs_padded_no_n.txt")
    pheno_filename = os.path.join(protein_data_dir, "protein_pheno_padded.txt")

    return Params(
        seqs_filename,
        pheno_filename,
        "protein_n",
        "first",
        False,
        "",
        None,
        72,
    )


@pytest.fixture
def protein_first_order_customaa(protein_data_dir):
    seqs_filename = os.path.join(protein_data_dir, "insulin_seq_test.txt")
    pheno_filename = os.path.join(protein_data_dir, "insulin_phi_test.txt")
    custom_alphbt = os.path.join(protein_data_dir, "custom_aa.txt")

    return Params(
        seqs_filename,
        pheno_filename,
        "protein_n",
        "first",
        False,
        "",
        custom_alphbt,
        72,
    )
