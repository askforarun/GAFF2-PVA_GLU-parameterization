# Example PDB files

This directory contains generated example PDB files for the corrected PVA--GLU workflow.

Single-chain PVA examples:
- `PVA7_trim.pdb`: uncapped PVA strand with `n = 7`.
- `PVA17_trim.pdb`: uncapped PVA strand with `n = 17`.
- `PVA21_trim.pdb`: uncapped PVA strand with `n = 21`.
- `PVA25_trim.pdb`: uncapped PVA strand with `n = 25`.

Representative Packmol combined-system examples:
- `packed_system_pva17_glu150.pdb`: 300 uncapped PVA strands with `n = 17` and 150 GLU molecules.
- `packed_system_pva21_glu150.pdb`: 300 uncapped PVA strands with `n = 21` and 150 GLU molecules.
- `packed_system_pva25_glu150.pdb`: 300 uncapped PVA strands with `n = 25` and 150 GLU molecules.

The molecule ordering in the combined PDB files follows the workflow convention used by `amber_to_lammps.py`: all PVA molecules first, followed by all GLU molecules. These files are examples of the combined coordinate input used for LAMMPS data-file generation; regenerated Packmol configurations may differ in coordinates because Packmol placement is stochastic unless a seed is fixed.
