import pytest
from pathlib import Path
from collections import namedtuple
import os


Params = namedtuple(
    "Params", [
        "seqs_filename",
        "pheno_filename",
        "molecule",
        "sites",
        "dm",
        "pop_size",
        "poly_order",
        "precomputed",
        "out_dir",
    ])


@pytest.fixture
def test_data_dir():
    return os.path.join(Path(__file__).parent, "tests", "data")


@pytest.fixture
def protein_data_dir(test_data_dir):
    return os.path.join(test_data_dir, "protein")


@pytest.fixture
def nucleotide_data_dir(test_data_dir):
    return os.path.join(test_data_dir, "nucleotide")


@pytest.fixture
def nucleotide_two_sites_data_dir(nucleotide_data_dir):
    return os.path.join(nucleotide_data_dir, "second_order")


@pytest.fixture
def protein_expected_output_dir(protein_data_dir):
    return os.path.join(protein_data_dir, "expected_outputs")


@pytest.fixture
def protein_nopad_expected_output_dir(protein_expected_output_dir):
    return os.path.join(protein_expected_output_dir, "protein_seqs_nopad")


@pytest.fixture
def nucleotide_expected_output_dir(nucleotide_two_sites_data_dir):
    return os.path.join(nucleotide_two_sites_data_dir, "expected_outputs")


@pytest.fixture
def protein_seqs_no_padding(protein_data_dir):
    return os.path.join(protein_data_dir, 'protein_seqs_nopad.txt')


@pytest.fixture
def protein_pheno_no_padding(protein_data_dir):
    return os.path.join(protein_data_dir, 'protein_pheno_nopad.txt')


@pytest.fixture
def nucleotide_params_first_order(nucleotide_two_sites_data_dir):
    seqs_filename = os.path.join(
        nucleotide_two_sites_data_dir,
        'test_seqs_2sites_dna.txt'
    )
    pheno_filename = os.path.join(
        nucleotide_two_sites_data_dir,
        'trait_test_seqs_2sites_dna.txt'
    )

    return Params(
        seqs_filename,
        pheno_filename,
        'DNA',
        2,
        4,
        12,
        'first',
        False,
        '/tmp'
    )


@pytest.fixture
def protein_params_first_order(protein_data_dir):
    seqs_filename = os.path.join(
        protein_data_dir,
        'protein_seqs_nopad.txt'
    )
    pheno_filename = os.path.join(
        protein_data_dir,
        'protein_pheno_nopad.txt'
    )

    return Params(
        seqs_filename,
        pheno_filename,
        'protein',
        6,
        20,
        6,
        'first',
        False,
        '/tmp',
    )
