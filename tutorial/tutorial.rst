.. include:: includes.rst.txt
*********************************************
ortho_seqs
*********************************************
This document will walk you through the steps of how to run a dataset on ortho_seqs, and what the various outputs are.

.. _downloads:
1. Setting Up Your Computer to Run ortho_seqs
-----------------------------------------------------------

The first thing you have to do (aside from gathering data!) is set up your computer to run ortho_seqs.

  You first need to have Miniconda installed on your computer, in order to do the shell commands. To do so, follow the link `here <https://docs.conda.io/en/latest/miniconda.html>`_, and choose the appropriate version, with regards to your computer.

After you have installed Minoconda, open up Terminal, or an equivalent Command-Line Interface (CLI). Run either this:

.. code-block:: shell-session

  conda create -n ortho_seqs pip
  pip install -r requirements.txt
  conda activate ortho_seq

Or, alternatively:

.. code-block:: shell-session

  conda env create -f conda_environment.yml
  conda activate ortho_seq

To activate ortho_seqs on your device. You will also need to run:

.. code-block:: shell-session

  python setup.py install

This line must be run every time ortho_seqs is updated, so you are using the most recent version. If the above steps have worked, congrats! You now have ortho_seqs on your computer. It's time to input some data.

.. _dataset_input:
2. Your Dataset
-----------------------------------------------------------

The data that is input to ortho_seqs must include a column of sequences, and a column of their corresponding phenotype values. These two columns can either be separate .txt files, or a single .xlsx or .csv file. Take, for instance, our toy example, which is a dataset originating from a paper titled `The Intrinsic Contributions of Tyrosine, Serine, Glycine, and Arginine to the Affinity and Specificity of Antibodies <https://www.sciencedirect.com/science/article/pii/S0022283608001691?casa_token=Qs608NJVJggAAAAA:-PruJ8_0_3pBtf4NHSVo0POYtzErFcDoqJYMxJQZER51_uZNtRYvBoWIMa9j3oIZJ18uY0rS3g>`_ by Sidhu et al., and measures the Specificity ELISA Signal Optical Density, with regards to its sequences (Figure 4a) (will be referred to as the "Sidhu Dataset" for this tutorial). The dataset, when input into ortho_seqs, should look like

.. image:: sidhu_image_txt.png
  :width: 400px
  :height: 400px

Note that for .xlsx (and .csv) files, the first column must be the sequences, and the second column must be the phenotypes. In addition, there must not be any header names for any files.

.. _parameter_definitions:
3. Executing Ortho_Seqs
-----------------------------------------------------------

We now turn towards our CLI to execute ortho_seqs.
Using the Sidhu Dataset, our input would look like:

.. code-block:: shell-session

  ortho_seq orthogonal-polynomial ortho_seq_code/tests/data/nucleotide/onefile_tests/sidhu.xlsx --molecule protein --poly_order first --out_dir ../onefile_tests/sidhu --alphbt_input SYG,R --min_pct 40

Let's explore what these flags are, and how you can use them to assist you.

The file input (ortho_seq_code/.../sidhu.xlsx) is our sequence AND phenotype data.

.. code-block::

  --molecule

This flag is where you indicate what kind of molecule this is. This can be DNA, RNA, or protein, as of now. For the Sidhu Dataset, the molecules are protein molecules.

.. code-block::

  --poly_order

This flag is to indicate the highest degree of polynomial order you want. Currently, DNA and RNA can go up to 2, and protein can only be 1. For the Sidhu Dataset, we are only interested in first-order calculations.

.. code-block::

  --pheno_file

This flag is not in the example, because we don't need it. If you were to present your data as two separate .txt files, then this would be where you put the file path for the phenotype data, and the first file path is for your sequence data.

.. code-block::

  --out_dir

This flag indicates where you want the output files to go (more on what exactly is saved there later). If the folder path already exists, ortho_seqs will create a new directory with a very similar name, and it will tell you what the new path's name is.

.. code-block::

  --alphbt_input

(Note: "Characters" in the following section refer to the nucleotides for DNA, the bases for RNA, and all 21 amino acids for proteins, plus one additional character, "n", which indicates nothing is at that spot)
This flag indicates the groupings of characters you want. The default will be no groupings, or every character gets counted on its own. If you include (uppercase) letters here, then only those characters will be noted (every other character, except "n", gets converted to a "z" and treated as one group). If you comma-separate somewhere in that group, then characters will be grouped based on what comma(s) they are in between. For the Sidhu Dataset, the groupings will be:

1. SYG \
2. R \
3. Everything else (z) \
4. n \

If I were to leave out the commas, the groups would be:

1. S
2. Y
3. R
4. G
5. Everything else (z)
6. n

.. code-block::

  --min_pct

One output will be a .xlsx file containing all of the first-order covariances. However, this file can get pretty big pretty quick, and will probably have a lot of covariance values of zero. Therefore, this flag will only print out covariance values whose magnitudes are at or above the PERCENTILE value specified. The default is 75, meaning it will only save the covariances who range from the 75th to the 100th percentiles in magnitude. To keep it at the default, leave out this flag when inputting what you want. For the Sidhu Dataset, we want all magnitudes at or above the 40th percentile.

.. _outputs:
4. Obtained Outputs
-----------------------------------------------------------

Insert text here
