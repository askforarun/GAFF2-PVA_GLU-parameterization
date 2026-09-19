# Regenerating parametrization outputs

`pva_builder.py` and `system_constants.py` were corrected on 2026-09-19 to
build the true fully-hydrolyzed PVA repeat unit (`-CH2-CHOH-` alternating,
one hydroxyl per two backbone carbons) instead of the previous
`-CH2-CHOH-CH2-` pattern (one hydroxyl per three backbone carbons). See
`pva_builder.py`'s module docstring for the corrected structure.

Because of this fix, every generated/parametrized file that was built from
the old code is now chemically stale and has been removed from version
control (see the removal commit). This includes:

- `combined_pva17.pdb`, `combined_pva21.pdb`, `combined_pva25.pdb`
- `PVA17_trim*.{pdb,mol2,frcmod,top,crd}`
- `pva17_glu.lammps`, `pva17_glu_parm.lammps`
- `charge_data/PVA7.pdb`, `charge_data/PVA7_min.pdb`, `charge_data/PVA7_min.mol2`
- `charge_data/glutaraldehyde*.{pdb,mol2,top,crd,frcmod}` and
  `charge_data/crosslinked_struct*.{pdb,mol2}` (the GLU/junction chemistry
  itself was already correct and unchanged, but these files were
  parametrized alongside the old PVA structure and should be regenerated
  for consistency)
- Antechamber scratch files (`ANTECHAMBER*`, `ATOMTYPE.INF`)

**Not removed** (already corrected in place this session, still valid):
- `charge_data/PVA_monomercharges.txt` — regenerated from the new 7-atom
  `CHOH-CH2` monomer definition (`PVA_ATOMS_PER_MONOMER = 7`)
- `charge_data/glutaraldehyde_charges.txt` — junction chemistry was already
  correct; charges unchanged

## To regenerate

```bash
conda activate AmberTools25
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization

# Reference structures used by charge_data/extract_charges.py
cd charge_data && python extract_charges.py && cd ..

# Per chain-length example outputs (repeat for n=17,21,25 as needed)
python example_parametrization.py --chain-length 17 --n-pva 300 --n-glu 150 \
  --combined-pdb combined_pva17.pdb --output-prefix pva17_glu
```

These outputs are gitignored going forward (see `.gitignore`) — regenerate
locally rather than committing them.
