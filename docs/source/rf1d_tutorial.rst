.. include:: includes.rst.txt
*********************************************
Tutorial: Visualizing your rFon1D results
*********************************************
This document will walk you through the steps of how to visualize your rFon1D outputs from *ortho_seqs* using the *rf1d-viz* CLI command.

  **Note:** *rf1d-viz* assumes that you have already run *orthogonal_polynomial* on the dataset. For a tutorial on how to run *orthogonal_polynomial*, view the tutorial `here <https://github.com/snafees/ortho_seqs/blob/plot/docs/source/orthogonal_polynomial_tutorial.rst>`_.

.. _Necessities:
1. Requirements for *rf1d-viz*
-----------------------------------------------------------

+ The *{trait_file_name}_regressions.npz* file that is returned from *orthogonal-polynomial*.

+ The *rf1d* form of the alphabet input.

When you run *orthogonal-polynomial*, the CLI will output the following text towards the beginning:

.. code-block:: shell-session

  rf1d form of alphabet input:

The line **beneath** that line is the *rf1d* form of the alphabet input.

+ The molecule type of the sequence (mostly *DNA* or *protein*).

+ What the phenotype values are representing.

.. _Flags:
2. *rf1d-viz* flags:
-----------------------------------------------------------

*rf1d-viz* will require you to input the following flags, many of which have counterparts in *orthogonal-polynomial*:

.. code-block:: shell-session

  --filename

This will be the *{trait_file_name}_regressions.npz* file that is returned from *orthogonal-polynomial*.

.. code-block:: shell-session

  --alphbt_input

This will be the *rf1d* form of the alphabet input.

.. code-block:: shell-session

  --molecule

This is the molecule type.

.. code-block:: shell-session

  --phenotype

This is the phenotype type. It will be used for labelling the graphs.

.. code-block:: shell-session

  --out_dir

This is where you want the graphs stored. **Note:** the path must exist prior to running *rf1d-viz*.

.. code-block:: shell-session

  --action

This is where you specify what kind of visualization you want. The current options are:

1. *barplot* - This will create a barplot of the rFon1D values, grouped by site and alphabet input. This is called automatically when you run *orthogonal-polynomial*.
1. *histogram* - This will create a histogram of the rFon1D values.
1. *summary* - Prints out the number of sites and dimensions, the alphabet input, the molecule, and calls *sort* (another *rf1d-viz* action that is explained in further detail below). This is called in *orthogonal-polynomial* automatically, and will not be saved.
1. *heatmap* - This will create a heatmap of the rFon1D values, grouped by site and alphabet input.
1. *boxplot* - This will create a boxplot of the rFon1D values, grouped by .
1. *sort* -

.. _input:
3. Running *rf1d-viz*
-----------------------------------------------------------



.. _example:
Guided example with the Sidhu dataset
-----------------------------------------------------------

The example uses the *Sidhu* dataset, which is the same as was used for the *orthogonal-polynomial* tutorial. Recall that the input for *orthogonal-polynomial* was:

.. code-block:: shell-session

  ortho_seq orthogonal-polynomial ortho_seq_code/Sidhu/Sidhu.xlsx --molecule protein --poly_order first --out_dir docs/source/tutorial_outputs --alphbt_input SYG,R --min_pct 40 --pheno_name IC50

The regression file that will be used for *rf1d-viz* will thus be called

.. code-block:: shell-session

  Sidhu_regressions.npz

Using the CLI output, we obtain

.. code-block:: shell-session

  rf1d form of alphabet input:
  SYG,R,z,n

which reveals that the *rf1d* form of the alphabet input is **SYG,R,z,n**.

With these in mind, the CLI input for *rf1d-viz* for a **barplot** will be

.. code-block:: shell-session

  ortho_seq rf1d-viz docs/source/tutorial_outputs/Sidhu_regressions.npz --alphbt_input SYG,R,z,n --molecule protein --phenotype IC50 --out_dir docs/source/tutorial_outputs --action barplot

This line of code will reproduce the graph that is automatically run, and looks like

.. image:: tutorial_outputs/rFon1D_Regressions_of_IC50_values.png
  :height: 250px
