# ortho_seqs
Ortho_seqs is a command line to convert sequence data (DNA/RNA/protein) to tensor-valued orthogonal polynomials and project phenotypes onto the polynomial space.
We do this by first converting the sequence information into 4-dimensional (for DNA/RNA) or 20-dimensional (for amino acids) vectors. The method can also be used for padded sequences to deal with unequal sequence lengths. 
Find out more about the approach in this paper [Analyzing genomic data using tensor-based orthogonal polynomials with application to synthetic RNAs](https://academic.oup.com/nargab/article/2/4/lqaa101/6030984). The paper gives an example of this method as applied to a case of synthetic RNA from a previously published dataset. Another manuscript detailing the use of this method to understand binding affinities of transcription factors (TFs) is currently in progress. The tool can also be used for protein sequence data and functionality to analyze amino acid sequences has been added to the CLI. 

# Directions for installing and running the tool are listed below.

## First, install an environment with dependencies for this package:
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

## Then, install the package:
``python setup.py install``

## Gather the input files needed.
1. You'll have one file with sequence data (.txt or .dat or .csv). See repo's data folder for examples of what these look like for DNA/RNA or proteins. 
2. You'll have one file with corresponding phenotypes (.txt or .dat or .csv) with the same length as the number of sequences. Phenotypes here are defined as real numbers.

## Then, to run the commandline tool:
To start with a test example, you can run the sample command below::

```
ortho_seq orthogonal-polynomial ./ortho_seq_code/tests/data/nucleotide/firstorder/first_order/test_seqs_2sites_dna.txt --pop_size 12 --dm 4 --sites 2 --moledule DNA --pheno_file ./ortho_seq_code/tests/data/nucleotide/firstorder/first_order/trait_test_seqs_2sites_dna.txt --poly_order second --out-dir ../results_ortho_seq_testing/DNA_2sites_test_run/
```
The above sample command line is building the tensor-valued orthogonal polynomial space based on the sequence data which consists of 12 sequences, each with two sites. Since these are DNA sequences, the vectors are 4-dimensional. Corresponding to each sequence is a phenotype value (a real number) as given in the phenotype file. For DNA, the tool can run first and second order analyses currently. We'll implement third order in a future version. For amino acids, the current version supports first order analysis and we hope to expand this in the future. 
Along with regressions on each site independent of one another and onto two sites at a time, the above command also computes *Fest* which is the phenotype estimated by the regressions. This shows that the mathematical calculations are done correctly as we now have an equation that accurately captures our initial data points. This only works here for sequences with 2 sites. If we had more sites, we'd need to do higher order calculations in order to capture all our combinations. Therefore, when runing the tool with more sites, as will probably be the case for most users, even just going up to second order gives us useful information about our system. First order tells us the importance of each site (independent of any correlations it might have with another site) and second order tells the importance of pairs of nucleotides independent of other pairs. Please take a look at the paper linked above to learn more about this method.  
