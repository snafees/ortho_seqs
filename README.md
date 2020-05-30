converting sequence data to tensor-valued orthogonal polynomials

# to install an environment with dependencies for this package
```
conda create -n ortho_seqs pip
pip install -r requirements.txt
conda activate ortho_seqs
```

or 

```
conda env create -f conda_environment.yml
conda activate ortho_seqs
```

# then to install the package
``python setup.py install``

# to run the commandline tool
To use the code, run the below command::

```
ortho_seq orthogonal-polynomial ./ortho_seq_code/testdata/pho4_r4_s1_site_3-4.csv --phenotype ./ortho_seq_code/testdata/pho4_r4_s1_ddg.txt  --pop-size 1000 --out-dir ../results_ortho_seq_testing
```