"""Shared structural constants for the hydrogel system."""

# One CHOH-CH2 repeat unit (C, H, O, H, C, H, H): matches fully hydrolyzed
# PVA's one-hydroxyl-per-two-backbone-carbons spacing.
PVA_ATOMS_PER_MONOMER = 7
GLU_ATOMS_PER_MOLECULE = 31
# Chosen to keep total atom count comfortably bounded at n=25 with
# n_pva = 2 * n_glu.
# Total atoms ~= n_glu * (2 * PVA_ATOMS_PER_MONOMER * 25 + GLU_ATOMS_PER_MOLECULE)
# => 150 * (2*7*25 + 31) = 150 * 381 = 57,150 atoms.
DEFAULT_N_GLU = 150
