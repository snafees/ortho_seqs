import numpy as np
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt
import ortho_seq_code.utils as utils
import itertools
import click
from ortho_seq_code.logo_colors import all_dna, all_rna, all_prot


def logo_plot(
    filename,
    molecule,
    out_dir,
):

    """Program to generate a logo plot of sequence data"""
    # ---- Import the sequence data. ----
    (dm, sites, pop_size, seq, alphabets, custom_aa, exc) = utils.get_seq_info(
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

    print(mean_df)

    # ---- Generate frequency logo plot. ----
    # generate file name
    seq_filename = filename.split("/")[-1]
    seq_file_prefix = seq_filename.split(".")[0]
    if out_dir[-1] != "/":
        out_dir = out_dir + "/"
    out_filename = out_dir + seq_file_prefix + "_freq_logo.png"

    # generate and save logo plot as png
    if molecule == "DNA":
        lm.Logo(mean_df, color_scheme=all_dna, stack_order="small_on_top")
        plt.savefig(out_filename)
    elif molecule == "RNA":
        lm.Logo(mean_df, color_scheme=all_rna, stack_order="small_on_top")
        plt.savefig(out_filename)
    elif molecule == "protein":
        lm.Logo(mean_df, color_scheme=all_prot, stack_order="small_on_top")
        plt.savefig(out_filename)
    else:
        print("Provide molecule type.")


@click.command(help="program to generate logo plots from sequence data")
@click.argument("filename", type=str)
@click.option(
    "--molecule", default="DNA", help="can provide DNA or amino acid sequence"
)
@click.option(
    "--out_dir",
    help="directory to save output/debug files to",
    type=str,
)  # noqa
def logo_plot_run(filename, molecule, out_dir):
    logo_plot(filename, molecule, out_dir)
