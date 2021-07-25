import numpy as np
import pandas as pd


def get_dsp(seqf):
    with open(seqf) as f:
        seq = f.readlines()
    seq_series_rm = pd.Series(seq).str.replace("\n", "")
    seq_series_nospace = seq_series_rm.str.replace(" ", "")
    seq_series = seq_series_nospace[seq_series_nospace != ""]
    sites = max(seq_series.str.len())
    pop_size = len(seq_series)
    seq_list = list(np.unique(list("".join(list(seq_series)))))
    dm = len(seq_list)
    for i in seq:
        if i == "\n":
            pop_size -= 1
            seq.remove(i)
    return [dm, sites, pop_size, seq]
