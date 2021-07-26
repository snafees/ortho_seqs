import numpy as np
import pandas as pd
from ortho_seq_code.constants_orthoseqs import *


def get_seq_info(seqf, alphbt_input, molecule):
    with open(seqf) as f:
        seq = f.readlines()
    seq_series_rm = pd.Series(seq).str.replace("\n", "")
    seq_series_nospace = seq_series_rm.str.replace(" ", "")
    seq_series = seq_series_nospace[seq_series_nospace != ""]
    sites = max(seq_series.str.len())
    pop_size = len(seq_series)
    seq_list = list(np.unique(list("".join(list(seq_series)))))
    print(seq)
    for i in seq:
        if i == "\n":
            pop_size -= 1
            seq.remove(i)

    # Autopadding lowercase n's
    if len(min(seq, key=len)) != len(max(seq, key=len)):
        incomplete_seq_series = seq_series[seq_series.str.len() < sites]
        while len(incomplete_seq_series) > 0:
            incomplete_seq_series = seq_series[seq_series.str.len() < sites]
            seq_series[incomplete_seq_series.index] += "n"
        seq_series += "\n"
        seq = list(seq_series)

    # Custom alphabet
    seq_oneline = "".join(seq_series)
    seq_list = list("".join(list(seq_series)))
    if alphbt_input is not None:
        alphbt = alphbt_input.upper()
        if alphbt == "PROTEIN_PNP":
            alphbt_input = "RNDCEQHKSTY,AGILMFPWV"
        elif alphbt == "ESSENTIAL":
            alphbt_input = "ILVFWHKTM,AGPYDERSCNQ"
        elif alphbt == "ACIDIC":
            alphbt_input = "DE,RHK,AGILPVFWYSTCMNQ"
        elif alphbt == "HYDROPHOBIC":
            alphbt_input = "AVLIPFC,"
        if "," in alphbt_input:
            alphbt = alphbt_input.upper()
            if "-" in alphbt:
                alphbt = alphbt.replace("-", "")
            # Adding on remaining letters as the last group
            alphbt_excluded = np.array([i for i in alphbt if i != ","])
            if "protein" in molecule:
                alphbt_last_group = "".join(
                    np.setdiff1d(
                        np.array(PROTEIN_ALPHABETS).ravel(), np.array(alphbt_excluded)
                    )
                )
            else:
                alphbt_last_group = "".join(
                    np.setdiff1d(
                        np.array(DNA_ALPHABETS).ravel(), np.array(alphbt_excluded)
                    )
                )
            if "-" in alphbt_input:
                alphbt_last_group += "n"
            alphbt += "," + str(alphbt_last_group)
            custom_aa = alphbt.split(",")
            if "" in custom_aa:
                custom_aa.remove("")
            # Assign group names to the group
            alphbt_count = 0
            aa_dict = dict()
            for i in range(len(custom_aa)):
                aa_dict[str(alphbt_count)] = list(np.unique(list(custom_aa[i])))
                if aa_dict[str(alphbt_count)] == []:
                    del aa_dict[str(alphbt_count)]
                    alphbt_count -= 1
                alphbt_count += 1
            if "n" in seq_list and "-" not in alphbt_input:
                aa_dict[str(alphbt_count)] = ["n"]
                custom_aa.append("n")
                alphbt_count += 1
            # Replaces amino acids with groups
            for i in range(len(seq_list)):
                for j in range(alphbt_count):
                    if seq_list[i] in aa_dict[str(j)]:
                        seq_list[i] = str(list(aa_dict.keys())[j])
            seq_list_sub = seq_list
            alphabets = list(aa_dict.keys())
            custom_dict = {alphabets[i]: custom_aa[i] for i in range(len(custom_aa))}

        else:
            alphabets = sorted(list(alphbt_input))
            alphabets_other = np.setdiff1d(np.array(seq_list), np.array(alphabets))
            if len(alphabets_other) > 0 and list(alphabets_other) != ["n"]:
                if "protein" in molecule:
                    alphbt_last_group = list(
                        np.setdiff1d(
                            np.array(PROTEIN_ALPHABETS).ravel(), np.array(alphabets)
                        )
                    )
                else:
                    alphbt_last_group = list(
                        np.setdiff1d(
                            np.array(DNA_ALPHABETS).ravel(), np.array(alphabets)
                        )
                    )
                seq_list_sub = []
                for i in range(len(seq_list)):
                    if seq_list[i] in alphbt_last_group:
                        seq_list_sub.append("n")
                    else:
                        seq_list_sub.append(seq_list[i])
            if "n" in seq_list_sub and "n" not in alphabets:
                alphabets.append("n")
            custom_aa = alphabets
        seq_adj = "".join(seq_list_sub)
        seq = [seq_adj[i : i + sites] for i in range(0, len(seq_adj), sites)]
    else:
        alphabets = list(np.unique(seq_list))
    while "" in alphabets:
        alphabets.rm("")
    while " " in alphabets:
        alphabets.rm(" ")
    while "\n" in alphabets:
        alphabets.rm("\n")
    dm = len(alphabets)
    return [dm, sites, pop_size, seq, seq_series, alphabets]
