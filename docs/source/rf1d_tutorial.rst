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

+ The molecule type (*DNA* or *protein*).

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

This will be the item groupings you specified in *orthogonal-polynomial*.

.. code-block:: shell-session

  conda install openpyxl
  python setup.py install

This line must be run every time ortho_seqs is updated, so you are using the most recent version. If the above steps have worked, congrats! You now have ortho_seqs on your computer. It's time to input some data.

.. _input:
3. Running *rf1d-viz*
-----------------------------------------------------------



.. _example:
Worked example with the Sidhu dataset
-----------------------------------------------------------
