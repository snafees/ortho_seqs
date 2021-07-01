from setuptools import find_packages, setup

with open("README.md") as f:
    long_description = f.read()

setup(
    author="Saba Nafees",
    author_email="saba.nafees314@gmail.com",
    description="A PyPI package to compute multivariate tensor-based orthogonal polynomials for sequence data and map phenotypes onto sequence space.",
    name="ortho_seqs",
    version="1.0.1",
    url="https://github.com/snafees/ortho_seqs",
    packages=find_packages(),
    keywords=[
        "regression",
        "projections" "covariance",
        "rna",
        "dna",
        "sequence-phenotype" "orthogonal polynomials",
        "higher order interactions",
    ],
    py_modules=["orthogonal_polynomial"],
    install_requires=["Click", "numpy"],
    include_package_data=True,
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={"console_scripts": ["ortho_seq = ortho_seq_code:cli"]},
)
