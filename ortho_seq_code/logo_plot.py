import click
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import PurePath

import logomaker as lm

import ortho_seq_code.utils as utils
from ortho_seq_code.constants_orthoseqs import (
    colors_for_dna_nucleotides,
    colors_for_rna_nucleotides,
    colors_for_amino_acids,
)


def logo_plot(
    filename,
    molecule,
    out_dir,
):

    """Program to generate a logo plot of sequence data"""
    # ---- Import the sequence data. ----
    (dm, sites, pop_size, seq, alphabets, _, _) = utils.get_seq_info(
        filename, None, molecule, True
    )

    # ---- Initialize terms that we will use. ----
    phi = np.zeros((sites, pop_size, dm))
    mean = np.zeros((sites, dm))
    range_sites = range(sites)
    range_popsize = range(pop_size)

    # ---- Calculate the means. ----
    for alphabet_index in range(dm):
        for i in range_popsize:
            for j in range_sites:
                if seq[i][j] == alphabets[alphabet_index]:
                    phi[j][i][alphabet_index] = 1.0

    for i, j in itertools.product(range_popsize, range_sites):
        mean[j] += phi[j][i] / pop_size

    mean_df = pd.DataFrame(mean, columns=alphabets)

    # ---- Generate frequency logo plot. ----
    seq_file_prefix = PurePath(filename).stem
    out_filename = seq_file_prefix + "_freq_logo.png"
    out_filepath = PurePath(out_dir).joinpath(out_filename)

    # generate and save logo plot as png
    if molecule == "DNA":
        color_scheme = colors_for_dna_nucleotides
    elif molecule == "RNA":
        color_scheme = colors_for_rna_nucleotides
    elif molecule == "protein":
        color_scheme = colors_for_amino_acids
    else:
        raise ValueError("Provide a valid molecule type (DNA, RNA, or protein).")

    lm.Logo(mean_df, color_scheme=color_scheme, stack_order="small_on_top")
    plt.savefig(out_filepath)


@click.command(help="program to generate logo plots from sequence data")
@click.argument("filename", type=click.Path(exists=True))
@click.option(
    "--molecule", default="DNA", help="can provide DNA, RNA, or amino acid sequence"
)
@click.option(
    "--out_dir",
    help="directory to save output/debug files to",
    type=str,
)  # noqa
def logo_plot_run(filename, molecule, out_dir):
    filename = click.format_filename(filename)
    logo_plot(filename, molecule, out_dir)
