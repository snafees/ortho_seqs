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
        "sites",
        "dm",
        "pop_size",
        "poly_order",
        "precomputed",
        "out_dir",
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
        2,
        4,
        12,
        "first",
        False,
        "",
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
        2,
        4,
        12,
        "second",
        False,
        "",
    )


@pytest.fixture
def protein_params_first_order(protein_data_dir):
    seqs_filename = os.path.join(protein_data_dir, "protein_seqs_nopad.txt")
    pheno_filename = os.path.join(protein_data_dir, "protein_pheno_nopad.txt")

    return Params(
        seqs_filename,
        pheno_filename,
        "protein",
        6,
        20,
        6,
        "first",
        False,
        "",
    )


@pytest.fixture
def protein_params_first_order_padded(protein_data_dir):
    seqs_filename = os.path.join(protein_data_dir, "protein_seqs_padded.txt")
    pheno_filename = os.path.join(protein_data_dir, "protein_pheno_padded.txt")

    return Params(
        seqs_filename,
        pheno_filename,
        "protein_n",
        6,
        21,
        10,
        "first",
        False,
        "",
    )
