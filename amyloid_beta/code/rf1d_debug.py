import numpy as np

data = np.load(
    "/Users/olivia.yoo/Desktop/code/ortho_seqs/amyloid_beta/ortho_seq_results/gk_charge_first/amyloid_gk_nostop_regressions.npz"
)

lst = data.files

for item in lst:
    print(item)
    print(data[item])
