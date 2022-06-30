import numpy as np
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt
import ortho_seq_code.utils as utils
from ortho_seq_code.constants_orthoseqs import DNA_ALPHABETS, PROTEIN_ALPHABETS
import itertools

# TODO: write this code into a function/functions ...

# arguments to put into function
input_file = "/Users/olivia.yoo/Desktop/code/ortho_seqs/logo_plot_assets/data/nucleotide/sample_dna_equal.xlsx"
input_molecule = "dna"

# ---- Importing the sequence data. ----
(dm, sites, pop_size, seq, alphabets, custom_aa, exc) = utils.get_seq_info(
    input_file, None, input_molecule, True
)

# print("dm")
# print(dm)
# print("sites")
# print(sites)
# print("pop_size")
# print(pop_size)
# print("seq")
# print(seq)
# print("alphabets")
# print(alphabets)
# print("custom_aa")
# print(custom_aa)
# print("exc")
# print(exc)

# ---- Initializing terms that we will use. ----
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

# find means
for i, j in itertools.product(range_popsize, range_sites):
    mean[j] += phi[j][i] / pop_size

print(mean)

mean_df = pd.DataFrame(mean, columns=alphabets)
print(mean_df)

# create logo plots of the mean/frequency.
# lm.Logo(mean_df)
# plt.show()

# ---- Calculate the heights of the letters. (information content) ----
# initialize appropriate terms
if input_molecule == "dna":
    s = 4
elif input_molecule == "protein":
    s = 20
else:
    print(
        "Molecule type not supported."
    )  # possibly unnecessary? depends on the amount of error checking desired
    sys.exit(1)

e_n = 1 / np.log(2) * (s - 1) / (2 * pop_size)
H_mtx = np.zeros((sites,))
R_mtx = np.zeros((sites,))

# uncertainty matrix. H_mtx[i] = uncertainty value for position i
for i in range(sites):
    for letter in alphabets:
        f_bi = mean_df[letter][i]
        if f_bi != 0:
            new_H = f_bi * np.log2(f_bi)
        else:
            new_H = 0
        H_mtx[i] = H_mtx[i] + new_H

H_mtx = np.multiply(H_mtx, -1)

# information content matrix, R_mtx[i] = total information content of position i
for i in range(sites):
    R_mtx[i] = np.log2(s) - (H_mtx[i] + e_n)

# heights of the letters. heights_df[i][j] is the height of letter j at position i
heights_df = pd.DataFrame.copy(mean_df)
heights_df = heights_df.mul(R_mtx, axis=0)

# plot sequence logo
lm.Logo(heights_df)
plt.show()
