# Regenerating parametrization outputs

`pva_builder.py` and `system_constants.py` define the fully hydrolyzed PVA
repeat as an alternating `-CH2-CHOH-` backbone with one hydroxyl per two
backbone carbons. See `pva_builder.py`'s module docstring for the structure
and atom ordering used by the current workflow.

Generated/parametrized files can be regenerated from the current builder
before production use. This includes:

- `combined_pva17.pdb`, `combined_pva21.pdb`, `combined_pva25.pdb`
- `PVA17_trim*.{pdb,mol2,frcmod,top,crd}`
- `pva17_glu.lammps`, `pva17_glu_parm.lammps`
- `charge_data/PVA7.pdb`, `charge_data/PVA7_min.pdb`, `charge_data/PVA7_min.mol2`
- `charge_data/glutaraldehyde*.{pdb,mol2,top,crd,frcmod}` and
  `charge_data/crosslinked_struct*.{pdb,mol2}`
- Antechamber scratch files (`ANTECHAMBER*`, `ATOMTYPE.INF`)

**Reference charge files**:
- `charge_data/PVA_monomercharges.txt` — charges for the 7-atom
  `CHOH-CH2` monomer definition (`PVA_ATOMS_PER_MONOMER = 7`)
- `charge_data/glutaraldehyde_charges.txt` — charges for the GLU reference

## To regenerate

```bash
conda activate AmberTools25
cd /users/ass2009/sharedscratch/GAFF2-PVA_GLU-parameterization

# Reference structures used by charge_data/extract_charges.py
cd charge_data && python extract_charges.py && cd ..

# Per chain-length example outputs (repeat for n=17,21,25 as needed)
python example_parametrization.py --chain-length 17 --n-pva 300 --n-glu 150 \
  --combined-pdb combined_pva17.pdb --output-prefix pva17_glu
```

These outputs are gitignored going forward (see `.gitignore`) — regenerate
locally rather than committing them.
