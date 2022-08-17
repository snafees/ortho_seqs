from ortho_seq_code.constants_orthoseqs import colors, PROTEIN_ALPHABETS

# ---- DNA color schemes. ----
all_dna_colors = ['green', 'blue', 'yellow', 'red']
all_dna = dict(zip(['A', 'C', 'G', 'T'], all_dna_colors))
all_dna['n'] = 'beige'

# ---- RNA color schemes. ----
all_rna = dict(zip(['A', 'C', 'G', 'U'], all_dna_colors))
all_rna['n'] = 'beige'

# ---- Protein color schemes. ----
all_prot = dict(zip(PROTEIN_ALPHABETS, colors))
all_prot['n'] = 'beige'