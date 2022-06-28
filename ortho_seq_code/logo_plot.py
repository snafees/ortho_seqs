import numpy as np
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt
import ortho_seq_code.utils as utils
from ortho_seq_code.constants_orthoseqs import DNA_ALPHABETS, PROTEIN_ALPHABETS

# TODO: write this code into a function/functions ...

# CODE THAT MAKES A LOGO PLOT!
input_file = "/Users/olivia.yoo/Desktop/code/ortho_seqs/logo_plot_assets/data/nucleotide/sample_dna_dm.xlsx"

# import the sequence data
(
    dm,
    sites,
    pop_size,
    seq,
    alphabets,
    custom_aa,
    exc
) = utils.get_seq_info(input_file, None, "dna", True)

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

# convert sequence data into position frequency matrix (pandas dataframe) and create logo plot
# assuming that the DNA sequences are aligned
# TODO: error checking for unequal/nonaligned sequences
freq_data = {'A': [0] * sites, 'C': [0] * sites, 'G': [0] * sites, 'T': [0] * sites}

# i = sequence number
for i in range(pop_size):
    # j = position number
    for j in range(sites):
        base = seq[i][j]
        freq_data[base][j] = freq_data[base][j] + 1

pfm_df = pd.DataFrame(data=freq_data)

print(pfm_df)

# FREQUENCY LOGO PLOT
# lm.Logo(pfm_df)
# plt.show()

# convert position frequency matrix into position probability matrix (i.e. normalize)
ppm_df = pd.DataFrame.copy(pfm_df)
ppm_df = ppm_df / pop_size

# RELATIVE FREQUENCY LOGO PLOT
# lm.Logo(ppm_df)
# plt.show()

# create the height matrix for the sequence logo
# four different nucleotide bases (I think this should be the same as dm but I'm not sure)
s = 4
e_n = 1 / np.log(2) * (s - 1) / (2 * pop_size)

# uncertainty matrix (one value for each column)
H_mtx = np.zeros((sites, ))

for i in range(sites):
    for base in DNA_ALPHABETS:
        current_H = H_mtx[i]
        f_bi = ppm_df[base][i]
        if f_bi != 0:
            new_H = f_bi * np.log2(f_bi)
        else:
            new_H = 0
        H_mtx[i] = current_H + new_H

H_mtx = np.multiply(H_mtx, -1)

R_mtx = np.zeros((sites, ))
for i in range(sites):
    R_mtx[i] = np.log2(s) - H_mtx[i] + e_n

print(R_mtx)

print(ppm_df)

heights_df = pd.DataFrame.copy(ppm_df)
heights_df = heights_df.mul(R_mtx, axis=0)
print(heights_df)

# ACTUAL SEQUENCE LOGO PLOT
lm.Logo(heights_df)
plt.show()