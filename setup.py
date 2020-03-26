import os
from setuptools import find_packages, setup


# -------------------------------------------------------------

NAME = 'ortho_seq_code'
PACKAGES = find_packages()
META_PATH = os.path.join('ortho_seq_code', '__init__.py')
KEYWORDS = ['regression', 'covariance', 'rna']

setup(
    entry_points="""
        [console_scripts]
        ortho_seq=ortho_seq_code:cli
    """,
)