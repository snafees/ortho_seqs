import pytest
from pathlib import Path
import os


@pytest.fixture(scope="session")
def test_data_dir():
    return Path(__file__).parent / "tests" / "data"


@pytest.fixture(scope="session")
def protein_data_dir(test_data_dir):
    return os.path.join(test_data_dir, "protein")

@pytest.fixture(scope="session")
def nucleotide_data_dir(test_data_dir):
    return os.path.join(test_data_dir, "nucleotide")

@pytest.fixture(scope="session")
def nucleotide_two_sites_data_dir(nucleotide_data_dir):
    return os.path.join(nucleotide_data_dir, "two_sites")


@pytest.fixture(scope="session")
def protein_expected_output_dir(test_data_dir):
    return test_data_dir / "protein" / "expected_outputs"

@pytest.fixture(scope="session")
def protein_nopad_expected_output_dir(protein_expected_output_dir):
    return protein_expected_output_dir / "protein_seqs_nopad" 


@pytest.fixture(scope="session")
def nucleotide_expected_output_dir(test_data_dir):
    return test_data_dir / "nucleotide" / "two_sites" / "expected_outputs"


@pytest.fixture
def protein_seqs_no_padding(protein_data_dir):
    return os.path.join(protein_data_dir, 'protein_seqs_nopad.txt')

@pytest.fixture
def protein_pheno_no_padding(protein_data_dir):
    return os.path.join(protein_data_dir, 'protein_pheno_nopad.txt')
