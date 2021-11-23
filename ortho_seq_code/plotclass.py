import os
import numpy as np
from ortho_seq_code.constants_orthoseqs import *
import matplotlib.pyplot as plt


class rf1d:
    # Initialize rf1d object
    def __init__(self, ndarray, alphbt_input, molecule="protein", custom=False):
        try:
            self.x = ndarray
            self.x_flat = list(ndarray.flatten())
            self.s = ndarray.shape[0]
            self.d = ndarray.shape[1]
            self.ind = np.arange(self.s)
            self.num_dm = np.arange(self.d)
            self.alphbt_input = alphbt_input
            self.m = molecule
            self.is_custom = custom
            self.complist = ["<", ">", "<>", "><"]
        except:
            print(
                "Error: Please provide an ndarray object and molecule type when initializing."
            )

    def summary(self):
        print("rf1d Object:\n")
        print("Number of sites:", str(self.s))
        print("Number of dimensions:", str(self.d))
        print("Alphabet inupt:", str(self.alphbt_input))
        print("Molecule:", str(self.m) + "\n")
        print("Highest rFon1D magnitudes:")
        self.sort(by_magnitude=True)

    # rFon1D bar plot
    def plot_bar(
        self,
        xlab="Sequence Site",
        ylab=None,
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
            width = 1 / self.s
            sd = self.s * self.d

            # Re-vectorization with null values (converts to dict)
            dim_num = dict()
            for i in self.ind:
                dim_num[i] = [
                    data_null[j] for j in np.arange(self.d * i, self.d * i + self.d)
                ]
            # some_dim = [data_array_flat[i], i for i in range(0, 160, 4)]

            # Remove all null data
            # dim_loc is the indeces of all non-null data
            dim_na = dict()
            dim_loc = dict()
            for i in self.ind:
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
                dim_aa[i] = [data_null[j] for j in range(i, sd, self.d)]

            col_len = len(colors)
            alpb_d = dict()
            for i in self.num_dm:
                if any(i != 0 and i for i in dim_aa[i]):
                    alpb_d[i] = colors[i % col_len]
                    alpb_d[self.alphbt_input[i]] = alpb_d.pop(i)

            # Creating plots
            fig, ax = plt.subplots(figsize=(8, 6))
            dim = dict()
            pi = dict()
            for i in range(self.s + 1):
                ax.axvline(i, color="lightgray", linewidth=0.8, zorder=0)
            for i in self.ind:
                if len(dim_na[i]) == 0:
                    ln = 1
                else:
                    ln = 1 / len(dim_na[i])
                rn = np.arange(1 / ln)
                if include_borders:
                    thickness = 1 - self.d / 40
                else:
                    thickness = 0
                pi[i] = ax.bar(
                    x=i + np.array([j for j in rn]) * ln,
                    height=[j for j in dim_na[i]],
                    width=ln,
                    align="edge",
                    color=[colors[i % col_len] for i in list(dim_loc[i])],
                    edgecolor="black",
                    zorder=3,
                    linewidth=thickness,
                )
            ax.axhline(color="black", linewidth=0.64)

            ax.set_xticks(self.ind + width + 0.5)
            ax.set_xticklabels(np.arange(1, self.s + 1))

            color_map = [color for color in list(alpb_d.values())]
            markers = [
                plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
                for color in alpb_d.values()
            ]
            if self.d < 6:
                dim = self.d
            else:
                dim = self.d // 3

            ax.legend(
                markers, alpb_d.keys(), loc=1, ncol=dim, prop={"size": 100 / self.d},
            )
            ax.tick_params(width=0.8)
            ax.xaxis.label.set_size(32 - (self.s) / 2)
            # width of the tick and the size of the tick labels
            # Regressions of off values onto each site of target RNA (orthogonalized within)
            # plt.savefig('rFon1D_off_star.png', bbox_inches='tight')
            plt.xlabel("Sequence Site")
            # plt.title("")
            if ylab is None:
                ylab = "Regressions of phenotype onto each site and amino acid"
            plt.ylabel(ylab)
            plt.tight_layout()
            figure = ax.get_figure()
            if out_dir is not None:
                path_sav = "rFon1D_" + str(ylab) + ".png"
                figure.savefig(os.path.join(str(out_dir), path_sav), dpi=400)
                print(
                    "saved regression graph as",
                    str(os.path.join(str(out_dir), path_sav)),
                )

        else:
            print("Nothing to graph for rFon1D")

    def plot_heatmap(self):
        print("Coming soon!")

    def sort(self, n=10, by_magnitude=False, ascending=True):
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
                str(x[s, k])
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

    def plot_hist(
        self,
        bins=None,
        site=None,
        alphabet_item=None,
        omit_zeroes=True,
        bin_color=None,
        border=True,
    ):
        if alphabet_item is not None and site is not None:
            print("Only one of site and alphabet_item may not be null.")
            return
        if alphabet_item is not None and alphabet_item not in self.alphbt_input:
            print("Alphabet item must be one of:")
            print(self.alphbt_input)
            return
        if site is not None:
            if site > self.s or site < 0:
                print("Site must be between 0 and", str(self.s))
                return
        if alphabet_item is not None:
            a = self.alphbt_input.index(alphabet_item) * self.s
            x_red = y.x_flat[a : a + self.d]
        elif site is not None:
            x_red = x[site]
        else:
            x_red = self.x_flat
        if omit_zeroes:
            lst = list(x_red)
            x_red = np.array([i for i in lst if i != 0])
        if border:
            if bins is not None:
                width = 1 - bins / 200
            else:
                width = 1 - len(x_red) / 1200
        else:
            width = 0
        plt.hist(x_red, bins=bins, color=bin_color, lw=width, ec="black")
