"""Shared structural constants for the hydrogel system."""

# One CHOH-CH2 repeat unit (C, H, O, H, C, H, H): matches fully hydrolyzed
# PVA's one-hydroxyl-per-two-backbone-carbons spacing.
PVA_ATOMS_PER_MONOMER = 7
GLU_ATOMS_PER_MOLECULE = 31

# An uncapped PVA strand (cap=False in pva_builder.build_pva) is
# CH2-(CHOH-CH2)*n: one extra leading CH2 group (C, H, H -- 3 atoms) in
# addition to the n repeat units. Total atoms per uncapped chain is
# therefore PVA_ATOMS_PER_MONOMER * n + PVA_LEADING_CH2_ATOMS, not a bare
# PVA_ATOMS_PER_MONOMER * n.
PVA_LEADING_CH2_ATOMS = 3


def pva_atoms_per_chain(n: int) -> int:
    """Total atom count of one uncapped PVA strand with ``n`` repeat units."""
    return PVA_ATOMS_PER_MONOMER * n + PVA_LEADING_CH2_ATOMS


# Chosen to keep total atom count comfortably bounded at n=25 with
# n_pva = 2 * n_glu.
# Total atoms ~= n_glu * (2 * pva_atoms_per_chain(25) + GLU_ATOMS_PER_MOLECULE)
# => 150 * (2*178 + 31) = 150 * 387 = 58,050 atoms.
DEFAULT_N_GLU = 150
