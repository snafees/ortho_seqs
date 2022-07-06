import pandas as pd

# ---. Create file containing all mutant sequences and phenotypes. ----
# df = pd.read_csv('/Users/olivia.yoo/Desktop/code/ortho_seqs_work/amyloid_beta/data/amyloid_data_mavenn.csv')
# df.drop(['set'], axis=1, inplace=True)
# df.columns = ['num_mut', 'F', 'dF', 'seq']
# df.to_csv('/Users/olivia.yoo/Desktop/code/ortho_seqs_work/amyloid_beta/data/amyloid_data.csv', index=False)

# ---- Get information about the dataset. ----
my_df = pd.read_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs_work/amyloid_beta/data/amyloid_data.csv"
)
# total_mutants = len(my_df.index)
# print(total_mutants)
#
# single_mutants = len(my_df[my_df['num_mut'] == 1])
# print(single_mutants)
#
# double_mutants = len(my_df[my_df['num_mut'] == 2])
# print(double_mutants)

# important constants for amyloid beta modules
wt_seq = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
n_sites = len(wt_seq)

# ---- Generate gatekeeper module dataset. ----


# def splice_gk(seq):
#     gk = seq[0] + seq[2] + seq[6] + seq[10] + seq[16] + seq[21] + seq[41]
#     return gk
#
#
# wt_gk = splice_gk(wt_seq)
# gk_df = my_df.copy()
# gk_df['seq'] = gk_df['seq'].apply(splice_gk)
# single_gk_df = gk_df[gk_df['num_mut'] == 1]
# single_gk_df = single_gk_df[single_gk_df['seq'] != wt_gk]
# single_gk_df.drop(['num_mut'], axis=1, inplace=True)
# single_gk_df.to_csv('/Users/olivia.yoo/Desktop/code/ortho_seqs_work/amyloid_beta/data/amyloid_single_gk.csv', index=False)
# n_gk_mutants = len(single_gk_df)


# ---- Generate N-terminal module dataset. ----


def splice_nterm(seq):
    nterm = seq[1] + seq[3:6] + seq[7:10] + seq[11:16] + seq[17:21] + seq[22:26]
    return nterm


wt_nterm = splice_nterm(wt_seq)
nterm_df = my_df.copy()
nterm_df["seq"] = nterm_df["seq"].apply(splice_nterm)
single_nterm_df = nterm_df[nterm_df["num_mut"] == 1]
single_nterm_df = single_nterm_df[single_nterm_df["seq"] != wt_nterm]
single_nterm_df.drop(["num_mut"], axis=1, inplace=True)
single_nterm_df.to_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs_work/amyloid_beta/data/amyloid_single_nterm.csv",
    index=False,
)
n_nterm_mutants = len(single_nterm_df)
#
print(n_nterm_mutants)


# ---- Generate C-terminal module dataset. ----


def splice_cterm(seq):
    cterm = seq[26:41]
    return cterm


wt_cterm = splice_cterm(wt_seq)
cterm_df = my_df.copy()
cterm_df["seq"] = cterm_df["seq"].apply(splice_cterm)
single_cterm_df = cterm_df[cterm_df["num_mut"] == 1]
single_cterm_df = single_cterm_df[single_cterm_df["seq"] != wt_cterm]
single_cterm_df.drop(["num_mut"], axis=1, inplace=True)
single_cterm_df.to_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs_work/amyloid_beta/data/amyloid_single_cterm.csv",
    index=False,
)
n_cterm_mutants = len(single_cterm_df)

print(n_cterm_mutants)
