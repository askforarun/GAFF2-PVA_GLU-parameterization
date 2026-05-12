# GAFF2 PVA Parameterization Pipeline

This repository contains the complete parametrization workflow for PVA (Polyvinyl Alcohol) using GAFF2 force field. The pipeline consists of three main stages:

## Overview

The parametrization workflow converts a hard-coded PVA polymer structure into a fully parametrized LAMMPS-compatible system with GAFF2 force field parameters and partial charges.

### Three Main Steps

| Step | Module | Input | Output | Purpose |
|------|--------|-------|--------|---------|
| **1. Build PVA with hard-coded geometry** | `pva_builder.py` | Chain length `n` | `PVA{n}_trim.pdb` | Generate uncapped PVA chain using fixed bond lengths and tetrahedral angles |
| **2. Run Antechamber + Parmchk2** | `molecular_utils.py::getfiles()` | `PVA{n}_trim.pdb` | `.mol2`, `.frcmod`, `.top`, `.crd` | Assign GAFF2 atom types, generate force field parameters (bonds, angles, dihedrals) |
| **3. Assign pre-extracted charges** | `molecular_utils.py::load_system_charges()` | Reference charges file | Per-atom partial charges | Reuse charges from minimized reference (PVA7_min.pdb) |

---

## Step 1: Build PVA with Hard-Coded Geometry

**File:** `pva_builder.py`

**Function:** `build_pva(n, output_file, cap=False)`

**Description:**
Builds a PVA polymer chain with the structure: `CH3-(CH2-CHOH-CH2)n-CH3`

The geometry is **hard-coded** (no external minimization), using:
- Fixed bond lengths (C-C: ~1.54 Å, C-O: ~1.42 Å, C-H: ~1.09 Å)
- Tetrahedral angles (109.47°)
- Van der Waals-based dihedral placement

**Key Parameters:**
- `n`: Number of CH2-CHOH-CH2 repeat units (typical: 3–25 monomers)
- `cap`: If `False`, generates uncapped structure (removes final 8 atoms/terminal CH3 groups)
- `output_file`: Output PDB file name (default: `PVA{n}.pdb` or `PVA{n}_trim.pdb`)

**Example Usage:**
```bash
python3 << 'EOF'
from pva_builder import build_pva
atoms, bonds = build_pva(n=7, output_file="PVA7_trim.pdb", cap=False)
print(f"Built PVA with {len(atoms)} atoms")
EOF
```

**Output:**
- `PVA{n}_trim.pdb` — Uncapped PVA structure (ready for Antechamber)
- Contains CONECT records for bond information

**Note:** This structure is NOT minimized; it's a scaffold for subsequent parametrization.

---

## Step 2: Run Antechamber + Parmchk2 (GAFF2 Parametrization)

**Files:** `molecular_utils.py`, depends on AmberTools

**Function:** `getfiles(pdb_file)`

**Description:**
Converts raw PVA structure into fully parametrized AMBER/LAMMPS topology by:

1. **Antechamber** → Assigns GAFF2 atom types and generates initial partial charges
   ```bash
   antechamber -j 4 -at gaff2 -dr no -fi pdb -fo mol2 \
     -i PVA7_trim.pdb -o PVA7_trim.mol2
   ```

2. **Parmchk2** → Fills in missing force field parameters
   ```bash
   parmchk2 -i PVA7_trim.mol2 -o PVA7_trim.frcmod -f mol2 -a Y
   ```

3. **tLeap** → Generates AMBER topology file
   ```bash
   source leaprc.gaff2
   SUS = loadmol2 PVA7_trim.mol2
   check SUS
   loadamberparams PVA7_trim.frcmod
   saveamberparm SUS PVA7_trim.top PVA7_trim.crd
   ```

**Pipeline (from `genhydrogel.py`):**
```python
pva_mod, glu_mod = parameterise_molecules(polymer, crosslinker)
```

Inside `parameterise_molecules()`:
1. Calls `getfiles(polymer)` to run Antechamber/Parmchk2/tLeap
2. Hand-corrects atom types for PVA:
   ```python
   content = content.replace("ha", "hc").replace("c2", "c3")
   ```
3. Returns corrected MOL2 base names

**Outputs:**
- `PVA{n}_trim.mol2` — Parametrized PVA in MOL2 format
- `PVA{n}_trim.frcmod` — AMBER force field corrections
- `PVA{n}_trim.top` — AMBER topology file
- `PVA{n}_trim.crd` — AMBER coordinate file

**Requirements:**
- AmberTools (Antechamber, Parmchk2, tLeap)
- Install via: `conda install -c conda-forge ambertools`

---

## Step 3: Assign Pre-Extracted Partial Charges

**File:** `molecular_utils.py`

**Function:** `load_system_charges(chain_length, n_pva, n_glu, pva_charge_file=None, glu_charge_file=None)`

**Description:**
Reuses partial charges from pre-minimized reference structures rather than using Antechamber-generated charges. This provides consistency across all jobs and avoids redundant calculation.

**Reference Charge Files:**
- `charge_data/PVA_monomercharges.txt` — Partial charges for one PVA monomer (extracted from PVA7_min.pdb)
- `charge_data/glutaraldehyde_charges.txt` — Partial charges for GLU molecule (from glutaraldehyde.pdb)

**Format of Charge Files:**
```
Atom  Type  OrigCharge  CorrectedCharge
C1    c3    0.123       0.125
C2    c3    0.089       0.090
...
```

**Example Usage:**
```python
from molecular_utils import load_system_charges
charges = load_system_charges(chain_length=7, n_pva=2, n_glu=1)
# Returns 1-D numpy array: [PVA_chain_1 charges | PVA_chain_2 charges | GLU_1 charges]
```

**Why Pre-Extract Charges?**
1. **Consistency**: All jobs use the same reference charges (from minimized PVA7)
2. **Speed**: Avoids re-running Antechamber for charge assignment
3. **Accuracy**: Hard-coded geometries differ slightly; using pre-minimized reference ensures chemical consistency

**The Charge Array Layout:**
```
[ c1, h1, c2, h2, o2, h3, c3, h4, h5, ...  (n*10 atoms for monomer),
  ..., (repeated n_pva times),
  ..., (GLU charges, n_glu times) ]
```

---

## Complete Workflow: From PDB to LAMMPS Data

### Step-by-Step Example

```bash
#!/bin/bash

# Set up environment
conda activate AmberTools25

# Step 1: Build PVA polymer with hard-coded geometry
python3 << 'EOF'
from pva_builder import build_pva
build_pva(n=7, output_file="PVA7_trim.pdb", cap=False)
EOF

# Step 2: Parametrize with GAFF2 (Antechamber + Parmchk2)
python3 << 'EOF'
from molecular_utils import getfiles
getfiles("PVA7_trim.pdb")
# Generates: PVA7_trim.mol2, PVA7_trim.frcmod, PVA7_trim.top, PVA7_trim.crd
EOF

# Step 3: Load pre-extracted partial charges
python3 << 'EOF'
import numpy as np
from molecular_utils import load_system_charges

# Load charges for 1 PVA chain (7 monomers) and 1 GLU
charges = load_system_charges(chain_length=7, n_pva=1, n_glu=1)
print(f"Total charges array: {len(charges)} atoms")
print(f"Charges: {charges}")
EOF
```

---

## File Dependencies

### Core Modules

- **`pva_builder.py`** (Step 1)
  - Standalone; no external dependencies beyond numpy
  
- **`molecular_utils.py`** (Steps 2 & 3)
  - Depends on: `system_constants.py`
  - Requires: AmberTools (Antechamber, Parmchk2, tLeap)
  
- **`system_constants.py`** (Configuration)
  - Hardcoded constants: `DEFAULT_N_GLU`, `PVA_ATOMS_PER_MONOMER`, `GLU_ATOMS_PER_MOLECULE`
  
- **`genhydrogel.py`** (Full orchestration, optional)
  - Integrates all three steps plus additional steps (LAMMPS conversion, crosslink parametrization)
  - Depends on: all above modules + `amber_to_lammps.py`

- **`amber_to_lammps.py`** (Optional, for LAMMPS conversion)
  - Converts AMBER topology → LAMMPS format

### Reference Data

- **`charge_data/PVA_monomercharges.txt`** — Pre-extracted PVA charges
- **`charge_data/glutaraldehyde_charges.txt`** — Pre-extracted GLU charges
- **`charge_data/PVA7_min.pdb`** — Minimized PVA reference (for charge extraction)
- **`charge_data/glutaraldehyde.pdb`** — Minimized GLU reference
- **`charge_data/crosslinked_struct_min.pdb`** — Minimized crosslink motif reference

---

## Directory Structure

```
GAFF2-PVA-parameterization/
├── pva_builder.py                    # Step 1: Build PVA with hard-coded geometry
├── molecular_utils.py                # Step 2 & 3: Parametrization + charge loading
├── genhydrogel.py                    # Full pipeline orchestrator (optional)
├── system_constants.py               # Configuration constants
├── amber_to_lammps.py                # AMBER→LAMMPS conversion (optional)
├── charge_data/
│   ├── PVA_monomercharges.txt        # Step 3: Pre-extracted PVA charges
│   ├── glutaraldehyde_charges.txt    # Step 3: Pre-extracted GLU charges
│   ├── PVA7_min.pdb                  # Reference minimized PVA
│   ├── PVA7_min.mol2                 # Reference parametrized PVA
│   ├── glutaraldehyde.pdb            # Reference minimized GLU
│   ├── crosslinked_struct_min.pdb    # Crosslink reference
│   └── extract_charges.py            # Utility to extract charges from PDB
└── README.md                         # This file
```

---

## Quick Start

### Minimal Example (Just the 3 Parametrization Steps)

```bash
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization

# Activate AmberTools environment
conda activate AmberTools25

# Run complete workflow
python3 << 'EOF'
import numpy as np
from pva_builder import build_pva
from molecular_utils import getfiles, load_system_charges

# Step 1: Build PVA
build_pva(n=7, output_file="PVA7_trim.pdb", cap=False)

# Step 2: Parametrize
getfiles("PVA7_trim.pdb")

# Step 3: Load charges
charges = load_system_charges(chain_length=7, n_pva=1, n_glu=1)
print(f"✓ Complete: {len(charges)} charges loaded")
EOF
```

---

## References

- **GAFF2 Force Field**: General AMBER Force Field v2
  - J. Wang et al., J. Comp. Chem. 2004 (GAFF)
  - Additional GAFF2 info: http://ambermd.org/
  
- **Antechamber**: Automated topology builder
  - Part of AmberTools
  
- **Hard-Coded PVA Geometry**: Fixed bond lengths and tetrahedral angles
  - Standard organic chemistry: C-C ~1.54 Å, C-O ~1.42 Å, C-H ~1.09 Å
  - Angles: sp³ tetrahedral ~109.47°

---

## Notes

1. **Why hard-coded geometry in Step 1?**
   - Fast and reproducible
   - Provides initial scaffold for Antechamber
   - Will be relaxed later during LAMMPS compression

2. **Why pre-extracted charges in Step 3?**
   - Consistent across all jobs
   - Avoids computational redundancy
   - Uses charges from minimized reference (PVA7_min.pdb)

3. **Antechamber-generated vs. pre-extracted charges?**
   - Antechamber charges from hard-coded geometry may differ slightly
   - Pre-extracted charges from minimized reference provide chemical consistency
   - Workflow uses pre-extracted for reproducibility

---

## Contact & Attribution

Pipeline developed for PVA-GLU hydrogel MD simulations using Signac workflow management.
