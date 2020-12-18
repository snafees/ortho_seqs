# ortho_seqs
Converting sequence data (DNA/RNA/protein) to tensor-valued orthogonal polynomials and projecting phenotypes onto the polynomial space.
Find out more about the approach in this paper [Analyzing genomic data using tensor-based orthogonal polynomials with application to synthetic RNAs](https://academic.oup.com/nargab/article/2/4/lqaa101/6030984) The paper gives an example of this method as applied to a case of synthetic RNA from a previously published dataset. Another manuscript detailing the use of this method to understand binding affinities of transcription factors (TFs) is currently in progress. The tool can also be used for protein sequence data and functionality to analyze amino acid sequences has been added to the CLI. 

# to install an environment with dependencies for this package
```
conda create -n ortho_seqs pip
pip install -r requirements.txt
conda activate ortho_seq
```

or 

```
conda env create -f conda_environment.yml
conda activate ortho_seq
```

# then to install the package
``python setup.py install``

# input files needed
1. You'll need a file with sequence data (.txt or .dat or .csv). 
2. You'll need a file with corresponding phenotypes (.txt or .dat or .csv) with the same length as the number of sequences.

# to run the commandline tool
To use the code, you can run the sample command below::

```
ortho_seq orthogonal-polynomial ./ortho_seq_code/testdata/pho4_r4_s1_site_3-4.csv --phenotype ./ortho_seq_code/testdata/pho4_r4_s1_ddg.txt  --pop-size 10000 --out-dir ../results_ortho_seq_testing
```
