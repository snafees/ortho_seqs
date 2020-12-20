DNA_ALPHABETS = ["A", "C", "G", "T"]
DNA_ALPHABETS_N = DNA_ALPHABETS + ["n"]
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
PROTEIN_ALPHABETS_N = PROTEIN_ALPHABETS + ["n"]
PROTEIN_ALPHABETS_POLAR = ["R", "N", "D", "C", "E", "Q", "H", "K", "S", "T", "Y"]
PROTEIN_ALPHABETS_NONPOLAR = ["A", "G", "I", "L", "M", "F", "P", "W", "V"]

DM_ALPHABETS = {
    4: DNA_ALPHABETS,
    5: DNA_ALPHABETS_N,
    20: PROTEIN_ALPHABETS,
    21: PROTEIN_ALPHABETS_N,
}
