# Example combined PDB files

This directory contains representative Packmol-generated combined PDB files for the corrected PVA--GLU workflow.

Files:
- `packed_system_pva17_glu150.pdb`: 300 uncapped PVA strands with `n = 17` and 150 GLU molecules.
- `packed_system_pva21_glu150.pdb`: 300 uncapped PVA strands with `n = 21` and 150 GLU molecules.
- `packed_system_pva25_glu150.pdb`: 300 uncapped PVA strands with `n = 25` and 150 GLU molecules.

The molecule ordering follows the workflow convention used by `amber_to_lammps.py`: all PVA molecules first, followed by all GLU molecules.
These files are examples of the combined coordinate input used for LAMMPS data-file generation; regenerated Packmol configurations may differ in coordinates because Packmol placement is stochastic unless a seed is fixed.
