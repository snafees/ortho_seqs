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
        except:
            print("Error: Please provide an ndarray object and molecule type when initializing.")

    # rFon1D bar plot
    def print(self):
        print("rf1d Object:\n")
        print("Number of sites:", str(self.s))
        print("Number of dimensions:", str(self.d))
        print("Alphabet inupt:", str(self.alphbt_input))
        print("Molecule:", str(self.m) + "\n")

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
                    alpb_d[self.alphabets[i]] = alpb_d.pop(i)

            print(alpb_d)
            # Creating plots
            fig, ax = plt.subplots()
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
            print("Color map:")
            print(color_map)
            markers = [
                plt.Line2D([0, 0], [0, 0], color=color_map, marker="o", linestyle="")
                for color in alpb_d.values()
            ]
            if self.d < 6:
                dim = self.d
            else:
                dim = self.d // 3

            ax.legend(
                markers, alpb_d, loc=1, ncol=self.d, prop={"size": 60 / self.d},
            )
            ax.tick_params(width=0.8, labelsize=80 / self.s)
            # width of the tick and the size of the tick labels
            # Regressions of off values onto each site of target RNA (orthogonalized within)
            # plt.savefig('rFon1D_off_star.png', bbox_inches='tight')

            plt.xlabel(
                xlab
                + "\nGroupings according to --alphbt_input:\n"
                + str(alpb_d)
                .replace("'", "")
                .replace(", ", " | ")
                .replace(": ", " is "),
                fontsize=4.8,
            )
            # plt.title("")
            if ylab is None:
                ylab = "Regressions of phenotype onto each site and amino acid"
            plt.ylabel(ylab)
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
