# Example PDB files

This directory contains generated example files for the PVA--GLU workflow.
The packed systems below are small demonstration inputs, not the current
large-system production boxes. See the main README for production counts.

Single-chain PVA examples:
- `PVA7_trim.pdb`: uncapped PVA strand with `n = 7`.
- `PVA17_trim.{pdb,mol2,frcmod,top,crd}`: uncapped PVA strand with `n = 17`; the no-cap strand contains 122 atoms.
- `PVA21_trim.{pdb,mol2,frcmod,top,crd}`: uncapped PVA strand with `n = 21`; the no-cap strand contains 150 atoms.
- `PVA25_trim.{pdb,mol2,frcmod,top,crd}`: uncapped PVA strand with `n = 25`; the no-cap strand contains 178 atoms.

The `*_trim` files are no-cap/reactive-strand examples used before
crosslinking. Their atom counts follow the fully hydrolyzed PVA
repeat, `CH2-CH(OH)`, with terminal methyl caps removed for reaction with GLU.

Small Packmol demonstration systems:
- `packed_system_pva17_glu150.pdb`: 300 uncapped PVA strands with `n = 17` and 150 GLU molecules.
- `packed_system_pva21_glu150.pdb`: 300 uncapped PVA strands with `n = 21` and 150 GLU molecules.
- `packed_system_pva25_glu150.pdb`: 300 uncapped PVA strands with `n = 25` and 150 GLU molecules.

The molecule ordering in the combined PDB files follows the workflow convention used by `amber_to_lammps.py`: all PVA molecules first, followed by all GLU molecules. These files are examples of the combined coordinate input used for LAMMPS data-file generation; regenerated Packmol configurations may differ in coordinates because Packmol placement is stochastic unless a seed is fixed.

## Current production-size packed inputs

These compressed PDB files each contain one packed, uncrosslinked starting
configuration for the current production counts. They are complete coordinate
files, not links to a separate simulation repository.

| File | PVA strands | GLU junctions | Atoms | Source configuration ID |
|---|---:|---:|---:|---|
| `packed_system_pva17_glu1350.pdb.gz` | 2700 | 1350 | 371,250 | `1180c166f2e323e3e0efb9efd5931fbe` |
| `packed_system_pva21_glu1050.pdb.gz` | 2100 | 1050 | 347,550 | `90a60f1bf1d090ebe85d2fc4030c00c7` |
| `packed_system_pva25_glu900.pdb.gz` | 1800 | 900 | 348,300 | `4e4d5e4ec34deb6ef6af84bbbca1f5e4` |

All PVA strands precede all GLU junctions, as in the small examples.
These are initial Packmol structures, not equilibrated or crosslinked networks.
Decompress before passing a file to the parameterization example:

```bash
gzip -dk examples/packed_system_pva17_glu1350.pdb.gz
python example_parametrization.py --chain-length 17 --n-pva 2700 --n-glu 1350 \
  --combined-pdb examples/packed_system_pva17_glu1350.pdb --output-prefix pva17_glu1350
```

Reference charge loading and AMBER-to-LAMMPS conversion remain distinct steps;
see the main README's charge-insertion note before using converted inputs.
