import pytest
from pathlib import Path
import os



@pytest.fixture(scope="session")
def test_data_dir():
    return Path(__file__) / "tests" / "data"

@pytest.fixture(scope="session")
def protein_expected_outputs_dir(test_data_dir):
    return test_data_dir / "protein" / "expected_outputs"

@pytest.fixture(scope="session")
def nucleotide_expected_outputs_dir(test_data_dir):
    return test_data_dir / "nucleotide" / "expected_outputs"
