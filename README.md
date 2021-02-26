# ortho_seqs
Ortho_seqs is a command line to convert sequence data (DNA/protein) to tensor-valued orthogonal polynomials and project phenotypes onto the polynomial space.
We do this by first converting the sequence information into 4-dimensional (for DNA) or 20-dimensional (for amino acids) vectors. The method can also be used for padded sequences to deal with unequal sequence lengths. 
Find out more about the approach in this paper [Analyzing genomic data using tensor-based orthogonal polynomials with application to synthetic RNAs](https://academic.oup.com/nargab/article/2/4/lqaa101/6030984). The paper gives an example of this method as applied to a case of synthetic RNA from a previously published dataset. Another manuscript detailing the use of this method to understand binding affinities of transcription factors (TFs) is currently in progress.  

For example, the sample data inputs for this tool are shown in this image. Here, each site in a sequence is first converted to a 4-dimensional vector. The input data includes phenotype values for each sequence.
![Figure showing sequence data that gets converted to vectors. Here, each sequence has a corresponding phenotype which is represented as a real number.](https://raw.githubusercontent.com/snafees/ortho_seqs/master/vec_methods_explanation-2.png)

# Usage

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
1. You'll have one file with sequence data as seen in the first column in the image above (.txt or .dat or .csv). See repo's data folder for examples of what these look like for DNA or proteins. 
2. You'll have one file with corresponding phenotypes as seen in the last column in the image above (.txt or .dat or .csv) with the same length as the number of sequences. Phenotypes here are defined as real numbers.

## Then, to run the commandline tool:
To start with a test example, you can run the sample command below::

```
ortho_seq orthogonal-polynomial ./ortho_seq_code/tests/data/nucleotide/firstorder/first_order/test_seqs_2sites_dna.txt --pop_size 12 --dm 4 --sites 2 --moledule DNA --pheno_file ./ortho_seq_code/tests/data/nucleotide/firstorder/first_order/trait_test_seqs_2sites_dna.txt --poly_order second --out-dir ../results_ortho_seq_testing/DNA_2sites_test_run/
```
The above sample command line is building the tensor-valued orthogonal polynomial space based on the sequence data which consists of 12 sequences, each with two sites. Since these are DNA sequences, the vectors are 4-dimensional. Corresponding to each sequence is a phenotype value (a real number) as given in the phenotype file. For DNA, the tool can run first and second order analyses currently. We'll implement third order in a future version. For amino acids, the current version supports first order analysis and we hope to expand this in the future. 
Along with regressions on each site independent of one another and onto two sites at a time, the above command also computes *Fest* which is the phenotype estimated by the regressions. This shows that the mathematical calculations are done correctly as we now have an equation that accurately captures our initial data points. This only works here for sequences with 2 sites. If we had more sites, we'd need to do higher order calculations in order to capture all our combinations. Therefore, when runing the tool with more sites, as will probably be the case for most users, even just going up to second order gives us useful information about our system. First order tells us the importance of each site (independent of any correlations it might have with another site) and second order tells the importance of pairs of nucleotides independent of other pairs. Please take a look at the paper linked above to learn more about this method.  

## Flags & Functionality
```
--pheno_file
```
Input a file with phenotype values corresponding to each sequence in the sequence file.
            
 ```
 --molecule
 ```
 Currently, you can provide DNA or protein sequences.
 ```
 --sites
 ```
 The number of sites in a sequence. If you have sequences with unequal lengths, please pad them with a lowercase 'n'. See examples in the ortho_seq_code/tests/ folder. 
 ```
 --dm
 ```
 The dimension of the vector corresponding to a site along a sequence. This is 4 for DNA and 20 for protein. For protein, you can provide a dimension of 21 in the case that you have padded protein sequences. See test data in the repo for examples.
```
--pop_size
```
The number of total sequences. 
```            
--poly_order
```
The order of the polynomials that will be constructed. Currently, one can do first and second order for DNA and first order for protein.
```
--out_dir
```
Directory where results can be stored.
```
--precomputed
```
Let's say you have a case where you have the same set of sequences but two different corresponding sets of phenotypes. You can build your sequence space and then project the first set of phenotypes onto this space. Then, if you wish to see how the other set of phenotypes maps onto the same sequence space, you can use this flag so that you're not wasting time and memory to recompute the space. When doing this, be sure to add your results from the first run to the **out_dir** when rerunning the command with the **precomputed** flag. 

# Support
If you have specific or general questions, feel free to open an issue and we'll do our best to address them.
If you have any comments, suggestions or would like to chat about this method or similar ideas, feel free to reach out via email at saba.nafees314@gmail.com.

# Roadmap
We hope to implement third order analysis for DNA in the near future. For amino acids, we hope to implement second order analysis. We'll add visualization ideas soon but if you have any thoughts on this, please feel free to reach out.

# Contribution
We hope to make the tool run faster as with higher dimenions and higher order analysis of longer sequence data, we can run into memory and time issues. Any thoughts on this or visualization are welcome. 

# Authors and acknowledgements
The derivation of the method and the construction of an initial version of the program was done by Dr. Sean Rice who served as Saba's Ph.D. advisor. Thank you to Isaac Griswold-Steiner for helping write the function to compute generalized inner and outer products. Thank you to Pranathi Vemuri for helping with the very initial draft of the CLI, adding CI integration testing, and to Phoenix Logan for helping write unit-tests. Their efforts are deeply appreciated!

# License
[MIT](https://choosealicense.com/licenses/mit/)


