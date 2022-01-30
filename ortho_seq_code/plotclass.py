import os
import numpy as np
import ortho_seq_code.constants_orthoseqs as constants
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import click
from ortho_seq_code.utils import create_dir_if_not_exists


class rf1d:
    # Initialize rf1d object
    def __init__(
        self, arr, alphbt_input, molecule="protein", phenotype=None, out_dir=None
    ):
        try:
            self.x = arr
            self.x_flat = list(arr.flatten())
            self.sites = arr.shape[0]
            self.dim = arr.shape[1]
            self.site_range = np.arange(self.sites)
            self.num_dm = np.arange(self.dim)
            self.alphbt_input = alphbt_input.split(",")[0 : arr.shape[1]]
            self.molecule = molecule
            self.complist = ["<", ">", "<>", "><"]
            self.phenotype = phenotype
            self.out_dir = out_dir
        except:
            print(
                "Error: Please provide a valid ndarray object and molecule type when initializing."
            )

    def summary(self):
        print("rf1d Object:\n")
        print("Number of sites:", str(self.sites))
        print("Number of dimensions:", str(self.dim))
        print("Alphabet input:", str(self.alphbt_input))
        print("Molecule:", str(self.molecule) + "\n")
        if self.phenotype is not None:
            print("Phenotype represents", self.phenotype, "values")
        if self.out_dir is not None:
            print("Image output directory:", self.out_dir)
        print("Highest rFon1D magnitudes:")
        self.sort(by_magnitude=True)

    # rFon1D bar plot
    def barplot(
        self,
        xlab="Sequence Site",
        title=None,
        fixed_width=False,
        include_borders=True,
        out_dir=None,
    ):
        if any(i != 0 for i in self.x_flat):
            if not fixed_width:
                data_null = np.where(
                    np.array(self.x_flat) == float(0), float("nan"), self.x_flat
                )
            else:
                data_null = self.x_flat

            # Constants/constant arrays
            width = 1 / self.sites
            sd = self.sites * self.dim

            # Re-vectorization with null values (converts to dict)
            dim_num = dict()
            for i in self.site_range:
                dim_num[i] = [
                    data_null[j]
                    for j in np.arange(self.dim * i, self.dim * i + self.dim)
                ]

            # Remove all null data
            # dim_loc is the indeces of all non-null data
            dim_na = dict()
            dim_loc = dict()
            for i in self.site_range:
                dim_na[i] = np.array(dim_num[i])[
                    np.array(np.isnan(dim_num[i])) == False
                ]
                dim_loc[i] = np.arange(len(dim_num[i]))[
                    np.array(np.isnan(dim_num[i])) == False
                ]
            # Color dictionary with corresponding letters

            dim_aa = dict()

            # dim_aa

            for i in self.num_dm:
                dim_aa[i] = [data_null[j] for j in range(i, sd, self.dim)]

            col_len = len(constants.colors)
            alpb_d = dict()
            for i in self.num_dm:
                if any(i != 0 and i for i in dim_aa[i]):
                    alpb_d[i] = constants.colors[i % col_len]
                    alpb_d[self.alphbt_input[i]] = alpb_d.pop(i)

            # Creating plots
            fig, ax = plt.subplots(figsize=(8, 6))
            dim = dict()
            pi = dict()
            for i in range(self.sites + 1):
                ax.axvline(i, color="lightgray", linewidth=0.8, zorder=0)
            for i in self.site_range:
                if len(dim_na[i]) == 0:
                    ln = 1
                else:
                    ln = 1 / len(dim_na[i])
                rn = np.arange(1 / ln)
                if include_borders:
                    thickness = 1 - self.dim / 40
                else:
                    thickness = 0
                pi[i] = ax.bar(
                    x=i + np.array([j for j in rn]) * ln,
                    height=[j for j in dim_na[i]],
                    width=ln,
                    align="edge",
                    color=[constants.colors[i % col_len] for i in list(dim_loc[i])],
                    edgecolor="black",
                    zorder=3,
                    linewidth=thickness,
                )
            ax.axhline(color="black", linewidth=0.64)

            ax.set_xticks(self.site_range + width + 0.5)
            ax.set_xticklabels(np.arange(1, self.sites + 1))
            markers = [
                plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
                for color in alpb_d.values()
            ]
            if self.dim < 6:
                dim = self.dim
            else:
                dim = self.dim // 3

            ax.legend(
                markers,
                alpb_d.keys(),
                loc="best",
                ncol=dim,
                prop={"size": 20 - sum([len(i) for i in self.alphbt_input]) / 2},
            )
            ax.tick_params(width=0.8)
            ax.xaxis.label.set_size(32 - (self.sites) / 2)
            # width of the tick and the size of the tick labels
            # Regressions of off values onto each site of target RNA (orthogonalized within)
            plt.xlabel("Sequence Site")
            if self.phenotype is None:
                ylab = "Regressions of phenotype onto each site and amino acid"
            else:
                ylab = "Regressions of " + self.phenotype + " values"
            plt.ylabel(ylab)
            plt.tight_layout()
            figure = ax.get_figure()
            if out_dir is not None:
                path_sav = "rFon1D_" + str(ylab) + ".png"
                path_sav = path_sav.replace(" ", "_")
                figure.savefig(os.path.join(str(out_dir), path_sav), dpi=400)
                print(
                    "saved regression graph as",
                    str(os.path.join(str(out_dir), path_sav)),
                )
            elif self.out_dir is not None:
                path_sav = "rFon1D_" + str(ylab) + ".png"
                path_sav = path_sav.replace(" ", "_")
                figure.savefig(os.path.join(str(self.out_dir), path_sav), dpi=400)
                print(
                    "saved regression graph as",
                    str(os.path.join(str(self.out_dir), path_sav)),
                )

        else:
            print("Nothing to graph for rFon1D")

    def sort(self, n=10, by_magnitude=True, ascending=True):
        if by_magnitude:
            x_flat = abs(np.array(self.x_flat))
            x = abs(self.x)
        else:
            x_flat = self.x_flat
            x = self.x
        x_sort = sorted(x_flat, reverse=ascending)[0:n]
        for i in x_sort:
            z = np.where(x == i)
            s = z[0][0]
            k = z[1][0]
            print(
                str(round(self.x[s, k], 4))
                + "\tSite: "
                + str(s)
                + "\t\tKey: "
                + str(self.alphbt_input[k])
            )

    def trim(self, span, comp):
        if comp not in self.complist:
            print(
                "comp parameter was not valid. Valid parameters are:\n"
                + str(self.complist)
            )
            return
        if type(span) not in [list, float, int]:
            print(
                "span parameter was not valid. Must be either a 2d list of numbers, or a number."
            )
            return
        if type(list) == list:
            if len(span) > 2 or not np.all(
                [type(i) == float or type(i) == int for i in self.x_flat]
            ):
                print("Span must be a 2d list of ints or floats.")
                return
        x_flat = np.array(self.x_flat)
        x = self.x
        if comp == "<":
            if type(span) == list:
                x[x < span[0]] = 0
                x_flat[x_flat < span[0]] = 0
            else:
                x[x < span] = 0
                x_flat[x_flat < span] = 0
        elif comp == ">":
            if type(span) == list:
                x[x > span[1]] = 0
                x_flat[x_flat > span[1]] = 0
            else:
                x[x > span] = 0
                x_flat[x_flat > span] = 0
        elif comp == "<>":
            if type(span) != list:
                print("<> and >< can only be used for list ranges.")
                return
            x[abs(x - np.mean(span)) >= np.mean(span)] = 0
            x_flat[abs(x_flat - np.mean(span)) >= np.mean(span)] = 0
        else:
            if type(span) != list:
                print("<> and >< can only be used for list ranges.")
                return
            x[abs(x - np.mean(span)) < np.mean(span)] = 0
            x_flat[abs(x_flat - np.mean(span)) < np.mean(span)] = 0
        print("Successfully trimmed array.")
        self.x_flat = list(x_flat)
        self.x = x

    def histogram(
        self,
        bins=None,
        site=None,
        alphabet_item=None,
        omit_zeroes=True,
        bin_color=None,
        border=True,
        out_dir=None,
    ):
        if alphabet_item is not None and site is not None:
            print("Only one of site and alphabet_item may not be null.")
            return
        if alphabet_item is not None and alphabet_item not in self.alphbt_input:
            print("Alphabet item must be one of:")
            print(self.alphbt_input)
            return
        if site is not None:
            if site > self.sites or site < 0:
                print("Site must be between 0 and", str(self.sites))
                return
        if alphabet_item is not None:
            a = self.alphbt_input.index(alphabet_item) * self.sites
            x_red = self.x_flat[a : a + self.dim]
        elif site is not None:
            x_red = self.x[site]
        else:
            x_red = self.x_flat
        if omit_zeroes:
            x_red = np.array(x_red)[np.array(x_red) != 0].tolist()
        if border:
            if bins is not None:
                width = 1 - bins / 200
            else:
                width = 1 - len(x_red) / 1200
        else:
            width = 0
        plt.hist(x_red, bins=bins, color=bin_color, lw=width, ec="black")
        plt.ylabel(self.phenotype)
        if out_dir is not None:
            path_sav = "rFon1D_hist_" + str(self.phenotype) or "" + ".png"
            path_sav = path_sav.replace(" ", "_")
            plt.savefig(os.path.join(str(out_dir), path_sav), dpi=400)
            print(
                "saved regression graph as", str(os.path.join(str(out_dir), path_sav)),
            )
        elif self.out_dir is not None:
            path_sav = "rFon1D_hist_" + str(self.phenotype) or "" + ".png"
            path_sav = path_sav.replace(" ", "_")
            plt.savefig(os.path.join(str(self.out_dir), path_sav), dpi=400)
            print(
                "saved regression graph as",
                str(os.path.join(str(self.out_dir), path_sav)),
            )

    def heatmap(self, out_dir=None):
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(self.x)
        ax.set_xticks(np.arange(len(self.alphbt_input)))
        ax.set_xticklabels(self.alphbt_input)
        ax.set_yticks(np.arange(self.sites))
        fig.tight_layout()
        for i in range(self.sites):
            for j in range(self.dim):
                text = ax.text(
                    j, i, round(self.x[i, j], 2), ha="center", va="center", color="w"
                )
        if out_dir is not None:
            path_sav = "rFon1D_heatmap_" + str(self.phenotype) or "" + ".png"
            path_sav = path_sav.replace(" ", "_")
            plt.savefig(os.path.join(str(out_dir), path_sav), dpi=400)
            print(
                "saved regression graph as", str(os.path.join(str(out_dir), path_sav)),
            )
        elif self.out_dir is not None:
            path_sav = "rFon1D_heatmap_" + str(self.phenotype) or "" + ".png"
            path_sav = path_sav.replace(" ", "_")
            plt.savefig(os.path.join(str(self.out_dir), path_sav), dpi=400)
            print(
                "saved regression graph as",
                str(os.path.join(str(self.out_dir), path_sav)),
            )
        plt.show()

    def boxplot(self, out_dir=None):
        fig, ax = plt.subplots()
        df = pd.melt(pd.DataFrame(self.x, columns=self.alphbt_input))
        df.columns = ["Item", self.phenotype]
        sns.boxplot(x=Item, y=self.phenotype, data=df)
        sns.stripplot(x="variable", y="value", color="black", data=df, alpha=0.8)
        if out_dir is not None:
            path_sav = "rFon1D_boxplot_" + str(self.phenotype) or "" + ".png"
            path_sav = path_sav.replace(" ", "_")
            plt.savefig(os.path.join(str(out_dir), path_sav), dpi=400)
            print(
                "saved regression graph as", str(os.path.join(str(out_dir), path_sav)),
            )
        elif self.out_dir is not None:
            path_sav = "rFon1D_boxplot_" + str(self.phenotype) or "" + ".png"
            path_sav = path_sav.replace(" ", "_")
            plt.savefig(os.path.join(str(self.out_dir), path_sav), dpi=400)
            print(
                "saved regression graph as",
                str(os.path.join(str(self.out_dir), path_sav)),
            )
        plt.show()

    def set_out_dir(self, new_out_dir):
        self.out_dir = new_out_dir


@click.command()
@click.argument("filename", type=str)
@click.option("--alphbt_input", type=str, help="A list form of the alphabet input")
@click.option("--molecule", type=str, default="DNA", help="Molecule type")
@click.option("--phenotype", type=str, help="What the phenotype values represent")
@click.option("--out_dir", type=str, help="The place where you want to save files")
@click.option(
    "--action",
    type=str,
    help="What you want to be output, can be one of 'summary', 'barplot', 'histogram', and 'heatmap'",
    default="barplot",
)
def rf1d_run(filename, alphbt_input, molecule, phenotype, out_dir, action):
    arr = np.load(filename)[[i for i in np.load(filename)][1]]
    x = rf1d(arr, alphbt_input, molecule, phenotype, out_dir)
    if action == "barplot":
        x.barplot()
    elif action == "histogram":
        x.histogram()
    elif action == "summary":
        x.summary()
    elif action == "heatmap":
        x.heatmap()
    elif action == "boxplot":
        x.boxplot()
    else:
        print("Please provide a valid action.")
