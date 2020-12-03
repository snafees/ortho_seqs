import pytest
from pathlib import Path
import os



@pytest.fixture(scope="session")
def test_data_dir():
    return Path(__file__) / "tests" / "data"
