import pytest
from pathlib import Path
import os



@pytest.fixture(scope="session")
def test_data_dir():
    return Path(__file__) / "tests" / "data"


@pytest.fixture(scope="session")
def protein_data_dir(test_data_dir):
    return os.path.join(test_data_dir, "protein")


@pytest.fixture(scope="session")
def protein_expected_output_dir(test_data_dir):
    return test_data_dir / "protein" / "expected_outputs"


@pytest.fixture(scope="session")
def nucleotide_expected_outputs_dir(test_data_dir):
    return test_data_dir / "nucleotide" / "expected_outputs"


@pytest.fixture
def protein_seqs_no_padding(protein_data_dir):
    return os.path.join(protein_data_dir, 'protein_seqs_nopad.txt')

@pytest.fixture
def protein_pheno_no_padding(protein_data_dir):
    return os.path.join(protein_data_dir, 'protein_pheno_nopad.txt')
