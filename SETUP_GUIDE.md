# Setup & Usage Guide

## Quick Start

### 1. Activate Environment
```bash
conda activate AmberTools25
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization
```

### 2. Run Example (All 3 Steps)
```bash
python example_parametrization.py
```

### 3. Output Files
After running the example, you'll have:
```
PVA7_trim.pdb           ← Step 1: Raw structure
PVA7_trim.mol2          ← Step 2: MOL2 format
PVA7_trim.frcmod        ← Step 2: Force field corrections
PVA7_trim.top           ← Step 2: AMBER topology
PVA7_trim.crd           ← Step 2: AMBER coordinates
```

Plus: Loaded charges for 2 PVA chains (7 monomers each) + 1 GLU molecule

---

## Files Included

### Core Parametrization Modules
| File | Purpose |
|------|---------|
| `pva_builder.py` | Step 1: Build PVA with hard-coded geometry |
| `molecular_utils.py` | Step 2 & 3: Parametrize with GAFF2, load charges |
| `system_constants.py` | Configuration constants |
| `genhydrogel.py` | Full integration (optional) |
| `amber_to_lammps.py` | AMBER→LAMMPS conversion (optional) |

### Reference Data (charge_data/)
| File | Purpose |
|------|---------|
| `PVA_monomercharges.txt` | Step 3: Pre-extracted PVA charges |
| `glutaraldehyde_charges.txt` | Step 3: Pre-extracted GLU charges |
| `PVA7_min.pdb` | Reference minimized PVA |
| `PVA7_min.mol2` | Reference parametrized PVA |
| `glutaraldehyde.pdb` | Reference minimized GLU |
| `extract_charges.py` | Utility to extract charges (if needed) |

### Documentation
| File | Purpose |
|------|---------|
| `README.md` | Complete overview |
| `PARAMETRIZATION_STEPS.md` | Detailed step-by-step breakdown |
| `SETUP_GUIDE.md` | This file |
| `example_parametrization.py` | Example script |

---

## Installation Requirements

### Required
```bash
# Install AmberTools (for Antechamber, Parmchk2, tLeap)
conda install -c conda-forge ambertools
```

### Optional (for full workflow)
```bash
# If running full hydrogel simulation workflow
conda install -c conda-forge signac signac-flow lammps gromacs
```

---

## Step-by-Step Workflow

### Step 1: Build PVA Structure

```python
from pva_builder import build_pva

# Build 7-monomer PVA chain
atoms, bonds = build_pva(
    n=7,                              # 7 monomers
    output_file="PVA7_trim.pdb",     # Output file
    cap=False                         # Uncapped (no terminal CH3)
)
```

**Output:** `PVA7_trim.pdb` (uncapped PVA structure)

**Characteristics:**
- Hard-coded geometry (no minimization)
- Ready for Antechamber parametrization
- Will be relaxed in LAMMPS later

---

### Step 2: Parametrize with GAFF2

```python
from molecular_utils import getfiles

# Run Antechamber + Parmchk2 + tLeap
getfiles("PVA7_trim.pdb")
```

**Requires:** AmberTools installed

**Process:**
1. `antechamber` → assigns GAFF2 atom types, initial charges
2. `parmchk2` → fills in missing force field parameters
3. `tleap` → generates AMBER topology

**Outputs:**
- `PVA7_trim.mol2` — Parametrized structure
- `PVA7_trim.frcmod` — Force field corrections
- `PVA7_trim.top` — AMBER topology
- `PVA7_trim.crd` — AMBER coordinates

---

### Step 3: Load Pre-Extracted Charges

```python
from molecular_utils import load_system_charges
import numpy as np

# Load charges for 2 PVA chains (7 monomers each) + 1 GLU
charges = load_system_charges(
    chain_length=7,   # 7 monomers per chain
    n_pva=2,          # 2 PVA chains
    n_glu=1           # 1 GLU molecule
)

# Returns: numpy array of shape (228,)
# Layout: [pva1_charges, pva2_charges, glu_charges]
```

**Reference Files Used:**
- `charge_data/PVA_monomercharges.txt` ← charges for 1 PVA monomer (10 atoms)
- `charge_data/glutaraldehyde_charges.txt` ← charges for 1 GLU (18 atoms)

**Why Pre-Extracted?**
- Consistent across all jobs
- From minimized reference (PVA7_min.pdb)
- Avoids re-running Antechamber per job

---

## Common Tasks

### Task 1: Parametrize Different Chain Length

```python
from pva_builder import build_pva
from molecular_utils import getfiles

# Build and parametrize 21-monomer PVA
build_pva(n=21, output_file="PVA21_trim.pdb", cap=False)
getfiles("PVA21_trim.pdb")

# Output: PVA21_trim.mol2, .frcmod, .top, .crd
```

### Task 2: Load Charges for Different System Size

```python
from molecular_utils import load_system_charges

# Load for 3 chains (25 monomers) + 2 GLU molecules
charges = load_system_charges(chain_length=25, n_pva=3, n_glu=2)
print(f"Total charges: {len(charges)}")  # (25*10*3) + (2*18) = 786
```

### Task 3: Check Charge Neutrality

```python
from molecular_utils import load_system_charges

charges = load_system_charges(chain_length=7, n_pva=2, n_glu=1)
net_charge = sum(charges)
print(f"Net charge: {net_charge:.6f}")  # Should be ≈ 0.0
```

### Task 4: Validate Charge Files

```python
from molecular_utils import _read_corrected_charges
from pathlib import Path

charge_data = Path("charge_data")

# Load PVA charges
pva_charges = _read_corrected_charges(str(charge_data / "PVA_monomercharges.txt"))
print(f"PVA monomer: {len(pva_charges)} atoms, total charge = {sum(pva_charges):.5f}")

# Load GLU charges
glu_charges = _read_corrected_charges(str(charge_data / "glutaraldehyde_charges.txt"))
print(f"GLU: {len(glu_charges)} atoms, total charge = {sum(glu_charges):.5f}")
```

---

## Troubleshooting

### Error: "antechamber not found"
**Solution:** Install AmberTools
```bash
conda install -c conda-forge ambertools
conda activate AmberTools25
```

### Error: "PVA_monomercharges.txt not found"
**Solution:** Ensure you're in the correct directory
```bash
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization
```

### Error: "charge file is empty"
**Solution:** Check that charge_data files are present
```bash
ls -la charge_data/
# Should show: PVA_monomercharges.txt, glutaraldehyde_charges.txt, etc.
```

### Output files not generated
**Solution:** Check for Antechamber errors in leap.log
```bash
cat leap.log
# Look for error messages from tLeap
```

---

## Next Steps

After parametrization, you can:

1. **Use in full hydrogel workflow** (from main repo)
   ```bash
   cd /mnt/scratch/users/ass2009/hydrogel_simulation
   python hydrogel_signac.py init
   python hydrogel_signac.py run -o generate_hydrogel
   ```

2. **Convert to LAMMPS** (using `amber_to_lammps.py`)
   ```python
   from src.amber_to_lammps import amber2lammps
   amber2lammps(
       data_file="data.lammps",
       param_file="parm.lammps",
       topologies=["PVA7_trim_mod.top"],
       pdb_file="packed_system.pdb",
       charges_target=[0],
   )
   ```

3. **Run LAMMPS simulations** (compression, crosslinking)
   ```bash
   mpirun -np 64 lmp -in in.compression_lammps
   ```

---

## References

- **GAFF2 Force Field**: Wang et al., J. Comp. Chem. 2004
- **Antechamber**: Part of AmberTools
- **Hard-Coded Geometry**: Standard organic chemistry bond lengths & angles

---

## Contact

For questions or issues, refer to:
- Main repo: `/mnt/scratch/users/ass2009/hydrogel_simulation/`
- CLAUDE.md: Full workflow documentation
