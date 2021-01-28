from setuptools import find_packages, setup

setup(
    name="ortho_seq_code",
    version="1.0",
    packages=find_packages(),
    keywords=[
        "regression",
        "covariance",
        "rna",
        "dna",
        "orthogonal polynomials",
        "higher order interactions",
    ],
    py_modules=["orthogonal_polynomial"],
    install_requires=["Click", "numpy"],
    include_package_data=True,
    entry_points={"console_scripts": ["ortho_seq = ortho_seq_code:cli"]},
)
