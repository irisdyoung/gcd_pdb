std_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
rare_amino_acids = ['FME', 'MSE', 'SEC', 'PYL', 'XAA', 'UNK'] # includes unknown or other unspecified
modified_amino_acids = ['4D4', 'D2T', 'MEQ']
amino_acids = std_amino_acids + rare_amino_acids + modified_amino_acids

std_nucleotides = ['A', 'U', 'T', 'C', 'G']
modified_nucleotides = ['5MU', 'OMG', 'OMU', 'G7M', 'OMC', '2MA', '2MG', '5MC',
'6MZ', '3TD', '1MG', 'PSU', 'MA6', 'H2U', 'UR3', '4OC', 'A2M']
nucleotides = std_nucleotides + modified_nucleotides

residues = amino_acids + nucleotides

disordered_solvent = ['HOH', 'H2O', 'WAT']
ordered_solvent = ['OOO']
solvent_molecules = disordered_solvent + ordered_solvent

metal_ions = ["ZN", "MG", "FE", "FE2"]

modified_all = modified_amino_acids + modified_nucleotides

name_to_three = {
  "ALANINE": "ALA",
  "ARGININE": "ARG",
  "ASPARAGINE": "ASN",
  "ASPARTIC ACID": "ASP",
  "CYSTEINE": "CYS",
  "GLUTAMIC_ACID": "GLU",
  "GLUTAMINE": "GLN",
  "GLYCINE": "GLY",
  "HISTIDINE": "HIS",
  "ISOLEUCINE": "ILE",
  "LEUCINE": "LEU",
  "LYSINE": "LYS",
  "METHIONINE": "MET",
  "PHENYLALANINE": "PHE",
  "PROLINE": "PRO",
  "SERINE": "SER",
  "THREONINE": "THR",
  "TRYPTOPHAN": "TRP",
  "TYROSINE": "TYR",
  "VALINE": "VAL",
  "N-FORMYLMETHIONINE": "FME", # nonstandard
  "SELENOMETHIONINE": "MSE", # nonstandard
  "SELENOCYSTEINE":"SEC", # nonstandard
  "PYRROLYSINE":"PYL", # nonstandard
  "UNKNOWN":"XAA" # unknown or unspecified
}

three_to_name = {
  "ALA": "ALANINE",
  "ARG": "ARGININE",
  "ASN": "ASPARAGINE",
  "ASP": "ASPARTIC ACID",
  "CYS": "CYSTEINE",
  "GLU": "GLUTAMIC_ACID",
  "GLN": "GLUTAMINE",
  "GLY": "GLYCINE",
  "HIS": "HISTIDINE",
  "ILE": "ISOLEUCINE",
  "LEU": "LEUCINE",
  "LYS": "LYSINE",
  "MET": "METHIONINE",
  "PHE": "PHENYLALANINE",
  "PRO": "PROLINE",
  "SER": "SERINE",
  "THR": "THREONINE",
  "TRP": "TRYPTOPHAN",
  "TYR": "TYROSINE",
  "VAL": "VALINE",
  "FME": "N-FORMYLMETHIONINE", # nonstandard
  "MSE": "SELENOMETHIONINE", # nonstandard
  "SEC": "SELENOCYSTEINE", # nonstandard
  "PYL": "PYRROLYSINE", # nonstandard
  "XAA": "UNKNOWN", # unknown or unspecified
  "UNK": "UNKNOWN"  # unknown or unspecified
}

three_to_one = {
  "ALA": "A",
  "ARG": "R",
  "ASN": "N",
  "ASP": "D",
  "CYS": "C",
  "GLU": "E",
  "GLN": "Q",
  "GLY": "G",
  "HIS": "H",
  "ILE": "I",
  "LEU": "L",
  "LYS": "K",
  "MET": "M",
  "PHE": "F",
  "PRO": "P",
  "SER": "S",
  "THR": "T",
  "TRP": "W",
  "TYR": "Y",
  "VAL": "V",
  "FME": "X",  # nonstandard, no 1-letter code. Assigned same code for unknown/unspecified.
  "MSE": "X",  # nonstandard, no 1-letter code. Assigned same code for unknown/unspecified.
  "SEC": "U",  # nonstandard
  "PYL": "O",  # nonstandard
  "XAA": "X",  # unknown or unspecified
  "UNK": "X"   # unknown or unspecified
}

one_to_three = {
  "A": "ALA",
  "R": "ARG",
  "N": "ASN",
  "D": "ASP",
  "C": "CYS",
  "E": "GLU",
  "Q": "GLN",
  "G": "GLY",
  "H": "HIS",
  "I": "ILE",
  "L": "LEU",
  "K": "LYS",
  "M": "MET",
  "F": "PHE",
  "P": "PRO",
  "S": "SER",
  "T": "THR",
  "W": "TRP",
  "Y": "TYR",
  "V": "VAL",
  "X": "XAA", # unknown, unspecified, or nonstandard amino acid without unique 1-letter code.
  "U": "SEC",
  "O": "PYL"
}

formal_charges = {
  "ALA": {},
  "ARG": {"NH1":+1},
  "ASN": {},
  "ASP": {"OD1":-1},
  "CYS": {},
  "GLU": {"OE1":-1},
  "GLN": {},
  "GLY": {},
  "HIS": {"ND1":+1},
  "ILE": {},
  "LEU": {},
  "LYS": {"NZ":+1},
  "MET": {},
  "PHE": {},
  "PRO": {},
  "SER": {},
  "THR": {},
  "TRP": {},
  "TYR": {},
  "VAL": {},
  "FME": {},
  "MSE": {},
  "SEC": {},
  "PYL": {},
  "BCR": {},
  "BCT": {"O1":-1},
  "CL":  {"CL":-1},
  "CLA": {"NB":-1, "ND":-1, "MG":+2},
  "DGD": {},
  "FE":  {"FE":+2},
  "FE2": {"FE":+2},
  "HEM": {"NA":-1, "NC":-1},
  "LBA": {"NB":-1, "ND":-1, "MG":+2},
  "LHG": {"O4":-1},
  "LMG": {},
  "OEX": {"MN1":"?", "MN2":"?", "MN3":"?", "MN4":"?", "CA1":"?", "O1":"?", "O2":"?", "O3":"?", "O4":"?", "O5":"?"},
  "PHO": {},
  "PL9": {},
  "SQD": {"O8":-1},
  "UNL": {}
}
