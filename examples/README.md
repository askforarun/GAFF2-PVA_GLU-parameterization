# Example PDB files

This directory contains generated example PDB files for the corrected PVA--GLU workflow.

Single-chain PVA examples:
- `PVA7_trim.pdb`: uncapped PVA strand with `n = 7`.
- `PVA17_trim.{pdb,mol2,frcmod,top,crd}`: uncapped PVA strand with `n = 17`; the corrected no-cap strand contains 122 atoms.
- `PVA21_trim.{pdb,mol2,frcmod,top,crd}`: uncapped PVA strand with `n = 21`; the corrected no-cap strand contains 150 atoms.
- `PVA25_trim.{pdb,mol2,frcmod,top,crd}`: uncapped PVA strand with `n = 25`; the corrected no-cap strand contains 178 atoms.

The `*_trim` files are no-cap/reactive-strand examples used before
crosslinking. Their atom counts follow the corrected fully hydrolyzed PVA
repeat, `CH2-CH(OH)`, with terminal methyl caps removed for reaction with GLU.

Representative Packmol combined-system examples:
- `packed_system_pva17_glu150.pdb`: 300 uncapped PVA strands with `n = 17` and 150 GLU molecules.
- `packed_system_pva21_glu150.pdb`: 300 uncapped PVA strands with `n = 21` and 150 GLU molecules.
- `packed_system_pva25_glu150.pdb`: 300 uncapped PVA strands with `n = 25` and 150 GLU molecules.

The molecule ordering in the combined PDB files follows the workflow convention used by `amber_to_lammps.py`: all PVA molecules first, followed by all GLU molecules. These files are examples of the combined coordinate input used for LAMMPS data-file generation; regenerated Packmol configurations may differ in coordinates because Packmol placement is stochastic unless a seed is fixed.
