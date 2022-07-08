import pandas as pd

# ---. Create file containing all mutant sequences and phenotypes. ----
df = pd.read_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs/amyloid_beta/data/amyloid_data_mavenn.csv"
)
df.drop(["set", "dy"], axis=1, inplace=True)
df.columns = ["num_mut", "F", "seq"]
df.to_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs/amyloid_beta/data/amyloid_data.csv",
    index=False,
)


# ---- Get information about the dataset. ----
my_df = pd.read_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs/amyloid_beta/data/amyloid_data.csv"
)

# all mutants
n_total_muts = len(my_df.index)
print(f"Total number of mutants: {n_total_muts}")

# single mutants
single_df = my_df[my_df["num_mut"] == 1].copy()
single_df.drop(["num_mut"], axis=1, inplace=True)
n_single_muts = len(single_df.index)
print(f"Total number of single AA mutants: {n_single_muts}")
n_single_stop_muts = len(single_df[single_df["seq"].str.contains("\*")])
print(f"Number of single AA mutants with stop codons: {n_single_stop_muts}")
print(
    f"Number of single AA mutants without stop codons: {n_single_muts - n_single_stop_muts}"
)

# double mutants
double_df = my_df[my_df["num_mut"] == 2].copy()
double_df.drop(["num_mut"], axis=1, inplace=True)
n_double_muts = len(double_df.index)
print(f"Total number of double AA mutants: {n_double_muts}")
n_double_stop_muts = len(double_df[double_df["seq"].str.contains("\*")])
print(f"Number of double AA mutants with stop codons: {n_double_stop_muts}")
print(
    f"Number of double AA mutants without stop codons: {n_double_muts - n_double_stop_muts}"
)


# ---- Set up for amyloid beta module separating. ----
wt_seq = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
n_sites = len(wt_seq)

no_stop_df = my_df.copy()
no_stop_df = no_stop_df[no_stop_df["seq"].str.contains("\*") == False]
single_df = no_stop_df[no_stop_df["num_mut"] == 1].copy()
double_df = no_stop_df[no_stop_df["num_mut"] == 2].copy()
single_df.drop(["num_mut"], axis=1, inplace=True)
single_df = single_df[["seq", "F"]]
double_df.drop(["num_mut"], axis=1, inplace=True)
double_df = double_df[["seq", "F"]]


def num_diff(seq1, seq2):
    count = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return count


# ---- Generate gatekeeper module dataset. ----


def splice_gk(seq):
    gk = seq[0] + seq[2] + seq[6] + seq[10] + seq[16] + seq[21] + seq[41]
    return gk


wt_gk = splice_gk(wt_seq)

# single gatekeeper mutants
single_gk_df = single_df.copy()
single_gk_df["seq"] = single_gk_df["seq"].apply(splice_gk)
single_gk_df = single_gk_df[single_gk_df["seq"] != wt_gk]

# double gatekeeper mutants
double_gk_df = double_df.copy()
double_gk_df["seq"] = double_gk_df["seq"].apply(splice_gk)
double_gk_seqs = list(double_gk_df["seq"])

wt_gk_list = [wt_gk] * (n_double_muts - n_double_stop_muts)
n_diffs = list(map(num_diff, double_gk_seqs, wt_gk_list))
double_gk_df["n_diffs"] = n_diffs
double_gk_df = double_gk_df[double_gk_df["n_diffs"] == 2]
double_gk_df.drop(["n_diffs"], axis=1, inplace=True)

gk_df = pd.concat([single_gk_df, double_gk_df], axis=0, join="outer")
gk_df.to_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs/amyloid_beta/data/amyloid_gk_nostop.csv",
    index=False,
    header=False,
)
print(f"Number of gatekeeper module-only mutants: {len(gk_df)}")


# ---- Generate N-terminal module dataset. ----


def splice_nterm(seq):
    nterm = seq[1] + seq[3:6] + seq[7:10] + seq[11:16] + seq[17:21] + seq[22:26]
    return nterm


wt_nterm = splice_nterm(wt_seq)

# single N-terminus mutants
single_nterm_df = single_df.copy()
single_nterm_df["seq"] = single_nterm_df["seq"].apply(splice_nterm)
single_nterm_df = single_nterm_df[single_nterm_df["seq"] != wt_nterm]

# double N-terminus mutants
double_nterm_df = double_df.copy()
double_nterm_df["seq"] = double_nterm_df["seq"].apply(splice_nterm)
double_nterm_seqs = list(double_nterm_df["seq"])

wt_nterm_list = [wt_nterm] * (n_double_muts - n_double_stop_muts)
n_diffs = list(map(num_diff, double_nterm_seqs, wt_nterm_list))
double_nterm_df["n_diffs"] = n_diffs
double_nterm_df = double_nterm_df[double_nterm_df["n_diffs"] == 2]
double_nterm_df.drop(["n_diffs"], axis=1, inplace=True)

nterm_df = pd.concat([single_nterm_df, double_nterm_df], axis=0, join="outer")
nterm_df.to_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs/amyloid_beta/data/amyloid_nterm_nostop.csv",
    index=False,
    header=False,
)
print(f"Number of N-terminus module-only mutants: {len(nterm_df)}")


# ---- Generate C-terminal module dataset. ----


def splice_cterm(seq):
    cterm = seq[26:41]
    return cterm


wt_cterm = splice_cterm(wt_seq)

# single C-terminus mutants
single_cterm_df = single_df.copy()
single_cterm_df["seq"] = single_cterm_df["seq"].apply(splice_cterm)
single_cterm_df = single_cterm_df[single_cterm_df["seq"] != wt_cterm]

# double C-terminus mutants
double_cterm_df = double_df.copy()
double_cterm_df["seq"] = double_cterm_df["seq"].apply(splice_cterm)
double_cterm_seqs = list(double_cterm_df["seq"])

wt_cterm_list = [wt_cterm] * (n_double_muts - n_double_stop_muts)
n_diffs = list(map(num_diff, double_cterm_seqs, wt_cterm_list))
double_cterm_df["n_diffs"] = n_diffs
double_cterm_df = double_cterm_df[double_cterm_df["n_diffs"] == 2]
double_cterm_df.drop(["n_diffs"], axis=1, inplace=True)

cterm_df = pd.concat([single_cterm_df, double_cterm_df], axis=0, join="outer")
cterm_df.to_csv(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs/amyloid_beta/data/amyloid_cterm_nostop.csv",
    index=False,
    header=False,
)
print(f"Number of C-terminus module-only mutants: {len(cterm_df)}")


# ---- Double check to make sure the right mutants are included. ----
gk_seqs = list(gk_df["seq"])
n_diffs = list(map(num_diff, gk_seqs, wt_gk_list))
print(
    f"Gatekeeper contains both single and double mutants and no others: {set(n_diffs) == set([1, 2])}"
)

nterm_seqs = list(nterm_df["seq"])
n_diffs = list(map(num_diff, nterm_seqs, wt_nterm_list))
print(
    f"N-terminus contains both single and double mutants and no others: {set(n_diffs) == set([1, 2])}"
)

cterm_seqs = list(cterm_df["seq"])
n_diffs = list(map(num_diff, cterm_seqs, wt_cterm_list))
print(
    f"C-terminus contains both single and double mutants and no others: {set(n_diffs) == set([1, 2])}"
)
