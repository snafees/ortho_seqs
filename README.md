# ortho_seqs
Converting sequence data (DNA/RNA/protein) to tensor-valued orthogonal polynomials and projecting phenotypes onto the polynomial space.

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

# input files needed
-You'll need a file with sequence data (.txt or .dat or .csv). 
-You'll need a file with corresponding phenotypes (.txt or .dat or .csv) with the same length as the number of sequences.

# to run the commandline tool
To use the code, you can run the sample command below::

```
ortho_seq orthogonal-polynomial ./ortho_seq_code/testdata/pho4_r4_s1_site_3-4.csv --phenotype ./ortho_seq_code/testdata/pho4_r4_s1_ddg.txt  --pop-size 1000 --out-dir ../results_ortho_seq_testing
```
