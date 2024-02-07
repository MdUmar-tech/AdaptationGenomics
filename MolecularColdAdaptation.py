from Bio.SeqUtils import ProtParam
from scipy.stats import ttest_1samp
import numpy as np

# Calculate Amino Acid Indices
def calculate_indices(protein_sequence):
    param_calc = ProtParam.ProteinAnalysis(str(protein_sequence))
    aa_percent = param_calc.get_amino_acids_percent()
    
    arg_to_lys_ratio = aa_percent['R'] / aa_percent['K'] if 'K' in aa_percent and aa_percent['K'] != 0 else 1  # Ensure non-zero denominator
    acidic_residue_freq = aa_percent['D'] + aa_percent['E']
    proline_residue_freq = aa_percent['P']
    aromaticity = param_calc.aromaticity()
    gravy = param_calc.gravy()
    
    return arg_to_lys_ratio, acidic_residue_freq, proline_residue_freq, aromaticity, gravy

# Example protein sequences for mesophilic proteins
mesophilic_protein_sequences = [
    "MKLVGFNGCSTRKTAGHESIYAINKLETVDTSNGKTFAWDVFLHSTQDLRFIPFQQMVGTSFLGHNTEVAVRDPALKLKDYTPEYISKRSSSIQEE",
    "MAEGNSGDSGQQQQQQQQTTSNNGASSHDLQQL",
    "MSPAMEAPMAAAAAGPLAGKPTPEQSDEEREQGDDEAEGDDSGGDRSPADGGEDRDNQ",
    "MKSLLLLFLIAGVGIPSAEQAAAAPFEQSFPDFSSWLS"
]
print(mesophilic_protein_sequences)
# Calculate indices for each mesophilic protein sequence
mesophilic_indices = []
for seq in mesophilic_protein_sequences:
    indices = calculate_indices(seq)
    mesophilic_indices.append(indices)

mesophilic_indices = np.array(mesophilic_indices)

# A hypothetical example of a protein sequence from Rhodococcus sp. JG3
rhodococcus_protein_sequence = "MKLVGFNGCSTRKTAGHESIYAINKLETVDTSNGKTFAWDVFLHSTQDLRFIPFQQMVGTSFLGHNTEVAVRDPALKLKDYTPEYISKRSSSIQEE"

# Calculate indices for the Rhodococcus sp. JG3 protein sequence
rhodococcus_indices = calculate_indices(rhodococcus_protein_sequence)

# Calculate average mesophilic indices
average_mesophilic_indices = np.mean(mesophilic_indices, axis=0)

# Perform one-sample t-test comparing the averaged indices of mesophilic proteins with the corresponding index of the protein from Rhodococcus sp. JG3
t_statistic, p_value = ttest_1samp(mesophilic_indices, rhodococcus_indices)

print("Averaged Mesophilic Indices:", average_mesophilic_indices)
print("T-Statistic:", t_statistic)
print("P-Value:", p_value)
