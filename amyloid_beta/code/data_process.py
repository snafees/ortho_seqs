import pandas as pd

# SET THIS DEPENDING ON WHERE YOU'RE RUNNING THE SCRIPT
location = "server"

if location == "local":
    ortho_dir = "/Users/olivia.yoo/Desktop/code/ortho_seqs/"
elif location == "server":
    ortho_dir = "/hpc/projects/data_lg/olivia.yoo/ortho_seqs/"
else:
    print("Set appropriate location.")
    exit(1)


# ---- Create file containing all mutant sequences and phenotypes. ----
df = pd.read_csv(ortho_dir + "amyloid_beta/data/amyloid_data_mavenn.csv")
df.drop(["set", "dy"], axis=1, inplace=True)
df.columns = ["num_mut", "F", "seq"]
df.to_csv(
    ortho_dir + "amyloid_beta/data/amyloid_data.csv",
    index=False,
)


# ---- Get information about the dataset. ----
my_df = pd.read_csv(ortho_dir + "amyloid_beta/data/amyloid_data.csv")

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


# ---- Pad the sequences with stop codons. ----
padded_df = my_df.copy()


def pad_seq(seq):
    try:
        stop_pos = seq.index("*")
        pad_len = len(seq) - stop_pos
        padding = "n" * pad_len
        padded_seq = seq[0:stop_pos] + padding
        return padded_seq
    except:
        return seq


padded_df["seq"] = padded_df["seq"].apply(pad_seq)
padded_df.drop(["num_mut"], axis=1, inplace=True)
padded_df = padded_df[["seq", "F"]]

padded_df.to_csv(
    ortho_dir + "amyloid_beta/data/amyloid_padded_data.csv",
    index=False,
)


# ----. Set up for separating amyloid beta modules. ----
wt_seq = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
n_sites = len(wt_seq)

separable_df = padded_df[padded_df["seq"].str.contains("n") == False].copy()
separable_df["not_seq"] = separable_df["seq"]


# ---- Generate gatekeeper dataset without stop codons. ----


def splice_gk(seq):
    gk = seq[0] + seq[2] + seq[6] + seq[10] + seq[16] + seq[21] + seq[41]
    return gk


def splice_not_gk(seq):
    not_gk = seq[1] + seq[3:6] + seq[7:10] + seq[11:16] + seq[17:21] + seq[22:41]
    return not_gk


wt_gk = splice_gk(wt_seq)
wt_not_gk = splice_not_gk(wt_seq)

gk_df = separable_df.copy()
gk_df["seq"] = gk_df["seq"].apply(splice_gk)
gk_df["not_seq"] = gk_df["not_seq"].apply(splice_not_gk)

gk_df.drop(gk_df[gk_df["not_seq"] != wt_not_gk].index, inplace=True)
gk_df.drop(gk_df[gk_df["seq"] == wt_gk].index, inplace=True)
gk_df.drop(["not_seq"], axis=1, inplace=True)

gk_df.to_csv(
    ortho_dir + "amyloid_beta/data/amyloid_gk_nostop.csv", index=False, header=False
)

print(f"Number of gatekeeper module-only mutants: {len(gk_df)}")

gk_stops = gk_df["seq"].str.contains("n").any()

print(f"Gatekeeper module contains padding: {gk_stops}")


# ---- Generate N-terminus dataset without stop codons. ----


def splice_nterm(seq):
    nterm = seq[1] + seq[3:6] + seq[7:10] + seq[11:16] + seq[17:21] + seq[22:26]
    return nterm


def splice_not_nterm(seq):
    not_nterm = seq[0] + seq[2] + seq[6] + seq[10] + seq[16] + seq[21] + seq[26:]
    return not_nterm


wt_nterm = splice_nterm(wt_seq)
wt_not_nterm = splice_not_nterm(wt_seq)

nterm_df = separable_df.copy()
nterm_df["seq"] = nterm_df["seq"].apply(splice_nterm)
nterm_df["not_seq"] = nterm_df["not_seq"].apply(splice_not_nterm)

nterm_df.drop(nterm_df[nterm_df["not_seq"] != wt_not_nterm].index, inplace=True)
nterm_df.drop(nterm_df[nterm_df["seq"] == wt_nterm].index, inplace=True)
nterm_df.drop(["not_seq"], axis=1, inplace=True)

nterm_df.to_csv(
    ortho_dir + "amyloid_beta/data/amyloid_nterm_nostop.csv", index=False, header=False
)

print(f"Number of N-terminus module-only mutants: {len(nterm_df)}")

nterm_stops = nterm_df["seq"].str.contains("n").any()

print(f"N-terminus module contains padding: {nterm_stops}")


# ---- Generate C-terminal module dataset. ----


def splice_cterm(seq):
    cterm = seq[26:41]
    return cterm


def splice_not_cterm(seq):
    not_cterm = seq[0:26] + seq[41]
    return not_cterm


wt_cterm = splice_cterm(wt_seq)
wt_not_cterm = splice_not_cterm(wt_seq)

cterm_df = separable_df.copy()
cterm_df["seq"] = cterm_df["seq"].apply(splice_cterm)
cterm_df["not_seq"] = cterm_df["not_seq"].apply(splice_not_cterm)

cterm_df.drop(cterm_df[cterm_df["not_seq"] != wt_not_cterm].index, inplace=True)
cterm_df.drop(cterm_df[cterm_df["seq"] == wt_cterm].index, inplace=True)
cterm_df.drop(["not_seq"], axis=1, inplace=True)

cterm_df.to_csv(
    ortho_dir + "amyloid_beta/data/amyloid_cterm_nostop.csv", index=False, header=False
)

print(f"Number of C-terminus module-only mutants: {len(cterm_df)}")

cterm_stops = cterm_df["seq"].str.contains("n").any()

print(f"C-terminus module contains padding: {cterm_stops}")


# ---- Generate N- and C-terminus datasets that include gatekeeper sites. ----


def splice_ngk(seq):
    return seq[0:26]


def splice_cgk(seq):
    return seq[26:]


wt_ngk = splice_ngk(wt_seq)
wt_cgk = splice_cgk(wt_seq)

split_df = separable_df.copy()
split_df["seq"] = split_df["seq"].apply(splice_ngk)
split_df["not_seq"] = split_df["not_seq"].apply(splice_cgk)

# Create N-terminus file.
ngk_df = split_df.copy()
ngk_df.drop(ngk_df[ngk_df["not_seq"] != wt_cgk].index, inplace=True)
ngk_df.drop(ngk_df[ngk_df["seq"] == wt_ngk].index, inplace=True)
ngk_df.drop(["not_seq"], axis=1, inplace=True)

ngk_df.to_csv(
    ortho_dir + "amyloid_beta/data/amyloid_ngk_nostop.csv", index=False, header=False
)

print(f"Number of N-terminus (including gatekeeper) only mutants: {len(ngk_df)}")

ngk_stops = ngk_df["seq"].str.contains("n").any()

print(f"N-terminus and gatekeeper module contains padding: {ngk_stops}")

# Create C-terminus file.
cgk_df = split_df.copy()

cgk_df.drop(cgk_df[cgk_df["seq"] != wt_ngk].index, inplace=True)
cgk_df.drop(cgk_df[cgk_df["not_seq"] == wt_cgk].index, inplace=True)
cgk_df.drop(["seq"], axis=1, inplace=True)
cgk_df = cgk_df[["not_seq", "F"]]

cgk_df.to_csv(
    ortho_dir + "amyloid_beta/data/amyloid_cgk_nostop.csv", index=False, header=False
)

print(f"Number of C-terminus (including gatekeeper) only mutants: {len(cgk_df)}")

cgk_stops = cgk_df["not_seq"].str.contains("n").any()

print(f"C-terminus and gatekeeper module contains padding: {cgk_stops}")
