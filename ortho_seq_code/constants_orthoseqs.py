DNA_ALPHABETS = ["A", "C", "G", "T"]
PROTEIN_ALPHABETS = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]

colors = [
    "tab:blue",
    "blueviolet",
    "saddlebrown",
    "tab:orange",
    "tab:green",
    "mediumvioletred",
    "coral",
    "tab:red",
    "tab:purple",
    "dodgerblue",
    "tab:brown",
    "gold",
    "tab:pink",
    "limegreen",
    "tab:gray",
    "chocolate",
    "tab:olive",
    "mediumvioletred",
    "goldenrod",
    "tab:cyan",
    "violet",
]

nucleotide_colors = [
    "green",
    "blue",
    "yellow",
    "red",
]

colors_for_dna_nucleotides = dict(zip(["A", "C", "G", "T"], nucleotide_colors))
colors_for_dna_nucleotides["n"] = "black"

colors_for_rna_nucleotides = dict(zip(["A", "C", "G", "U"], nucleotide_colors))
colors_for_dna_nucleotides["n"] = "black"

colors_for_amino_acids = dict(zip(PROTEIN_ALPHABETS, colors))
colors_for_amino_acids["n"] = "black"
