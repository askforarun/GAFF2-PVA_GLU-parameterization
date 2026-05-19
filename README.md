# GAFF2 PVA–GLU Parameterization Pipeline

This repository contains the complete parametrization workflow for PVA (Polyvinyl Alcohol) and GLU (Glutaraldehyde) using GAFF2 force field. The pipeline consists of four main stages, converting hard-coded molecular structures to fully parametrized LAMMPS-ready systems.

## Overview

The parametrization workflow converts hard-coded PVA polymer and GLU crosslinker structures into fully parametrized LAMMPS-compatible systems with GAFF2 force field parameters and pre-extracted partial charges.

### Four Main Steps (+ Optional GLU Preparation)

| Step | Module | Input | Output | Purpose |
|------|--------|-------|--------|---------|
| **1. Build PVA with hard-coded geometry** | `pva_builder.py` | Chain length `n` | `PVA{n}_trim.pdb` | Generate uncapped PVA chain using fixed bond lengths and tetrahedral angles |
| **1b. (Optional) Prepare GLU crosslinker** | Pre-minimized files | `glutaraldehyde.pdb` | `.mol2`, `.frcmod`, `.prmtop` (if needed) | Prepare tetrafunctional GLU crosslinker for parametrization |
| **2. Run Antechamber + Parmchk2** | `molecular_utils.py::getfiles()` | PVA and/or GLU PDB | `.mol2`, `.frcmod`, `.prmtop` | Assign GAFF2 atom types, generate force field parameters (bonds, angles, dihedrals) |
| **3. Assign pre-extracted charges** | `molecular_utils.py::load_system_charges()` | Reference charges files | Per-atom partial charges | Reuse charges from minimized references (PVA7_min.pdb, glutaraldehyde.pdb) |
| **4. Convert to LAMMPS format** | `amber_to_lammps.py` | `.prmtop` files + combined PDB | `system.lammps`, `parm.lammps` | Convert AMBER topology to LAMMPS data and parameter files for MD simulation |

---

## Step 1: Build PVA with Hard-Coded Geometry

**File:** `pva_builder.py`

**Function:** `build_pva(n, output_file, cap=False)`

**Description:**
Builds a PVA polymer chain with fixed geometry.

**Structure:**
- **With `cap=True`:** `CH3-(CH2-CHOH-CH2)n-CH3` (capped with terminal methyl groups)
- **With `cap=False`:** `(CH2-CHOH-CH2)n` (uncapped/trimmed backbone only)

The geometry is **hard-coded** (no external minimization), using:
- Fixed bond lengths (C-C: ~1.54 Å, C-O: ~1.42 Å, C-H: ~1.09 Å)
- Tetrahedral angles (109.47°)
- Van der Waals-based dihedral placement

**Key Parameters:**
- `n`: Number of CH2-CHOH-CH2 repeat units (typical: 3–25 monomers)
- `cap`: If `True`, generates capped structure with CH3 terminals; if `False`, generates uncapped backbone (removes terminal CH3 groups)
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

## Step 1b: Prepare Pre-Polymerized Glutaraldehyde (GLU) Crosslinker (Optional)

**File:** Pre-minimized structure available in `charge_data/glutaraldehyde.pdb`

**Description:**
Unlike PVA (which is generated algorithmically), the **pre-polymerized** glutaraldehyde crosslinker is a tetrafunctional molecule provided as a minimized reference structure. This is not monomeric glutaraldehyde, but a condensed oligomeric form designed to function as an ideal four-branch junction in the network.

**Structure:**
- **Pre-polymerized tetramer:** 4 reactive junction carbons that connect to PVA strand ends
- **Tetrafunctional connectivity:** Each GLU molecule serves as one junction node
- **Pre-minimized geometry:** Ensures consistent charges and conformations across all jobs

**Variants:**

1. **Saturated (fully bonded, no reactive sites):** `charge_data/glutaraldehyde.pdb`
   - Used for charge extraction and analysis
   - Has hydrogen atoms on all junction carbons
   - For reference only in charge/structural analysis
   
2. **Activated (reactive for crosslinking):** Derived by removing H atoms from junction carbons
   - Used in actual network construction (genhydrogel.py)
   - Defines the four reactive sites for C–C bonding with PVA strand ends
   - One H removed per junction carbon = 4 reactive sites total

**Usage:**
If parametrizing GLU separately, use the pre-minimized structure directly with Antechamber:
```bash
python3 << 'EOF'
from molecular_utils import getfiles
getfiles("charge_data/glutaraldehyde.pdb")
# Generates: glutaraldehyde.mol2, glutaraldehyde.frcmod, glutaraldehyde.prmtop
EOF
```

**Important Notes:**
1. **Pre-extracted charges are for unsaturated GLU:** The charges in `glutaraldehyde_charges.txt` are extracted for the **unsaturated form** (reactive junction carbons with H removed), but retain 31 total atoms.
2. **Saturated vs. Unsaturated:** 
   - **Saturated** (in `glutaraldehyde.pdb`): Reference structure with explicit H on all junction carbons — used for parametrization
   - **Unsaturated** (in `glutaraldehyde_charges.txt`): Charges for reactive form with H removed from junction carbons — used in network construction
3. **For most workflows, skip this step:** Use pre-extracted GLU charges from `charge_data/glutaraldehyde_charges.txt` directly (Step 3).

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
See the **Notes** section below for detailed explanation. Briefly: this ensures consistency across all jobs, avoids redundant Antechamber runs, and guarantees chemical accuracy since hard-coded geometries differ from minimized structures.

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
import sys
sys.path.insert(0, "/users/ass2009/sharedscratch/GAFF2-PVA-parameterization")
from pva_builder import build_pva

build_pva(n=7, output_file="PVA7_trim.pdb", cap=False)
print("✓ Step 1: PVA7_trim.pdb generated")
EOF

# Step 2: Parametrize PVA with GAFF2 (Antechamber + Parmchk2)
python3 << 'EOF'
import sys
sys.path.insert(0, "/users/ass2009/sharedscratch/GAFF2-PVA-parameterization")
from molecular_utils import getfiles

getfiles("PVA7_trim.pdb")
print("✓ Step 2a: PVA7_trim.top generated (70 atoms, 69 bonds)")

# Step 2b: Parametrize GLU crosslinker (required for combined systems)
getfiles("charge_data/glutaraldehyde.pdb")
print("✓ Step 2b: glutaraldehyde.top generated (31 atoms, 32 bonds)")
EOF

# Step 3: Load pre-extracted partial charges for both molecules
python3 << 'EOF'
import sys
sys.path.insert(0, "/users/ass2009/sharedscratch/GAFF2-PVA-parameterization")
from molecular_utils import load_system_charges

# Load charges for 1 PVA chain (7 monomers) and 1 GLU
charges = load_system_charges(chain_length=7, n_pva=1, n_glu=1)
print(f"✓ Step 3: Charges loaded for {len(charges)} atoms")
print(f"  PVA (70 atoms) + GLU (31 atoms) = {len(charges)} total")
EOF

# Create combined PDB with both molecules
python3 << 'EOF'
# Read PVA and GLU structures, combine them into one PDB file
with open("PVA7_trim.pdb") as f:
    pva_lines = f.readlines()

with open("charge_data/glutaraldehyde.pdb") as f:
    glu_lines = f.readlines()

# Count PVA atoms
pva_atoms = sum(1 for line in pva_lines if line.startswith(("ATOM", "HETATM")))

# Write combined PDB
combined = []
for line in pva_lines:
    if line.startswith(("ATOM", "HETATM", "CONECT")):
        break
    combined.append(line)

# Extract PVA ATOM records
for line in pva_lines:
    if line.startswith(("ATOM", "HETATM")):
        combined.append(line)

# Extract GLU ATOM records with renumbered indices
glu_atom_count = 0
for line in glu_lines:
    if line.startswith(("ATOM", "HETATM")):
        atom_num = pva_atoms + glu_atom_count + 1
        new_line = line[:6] + f"{atom_num:5d}" + line[11:]
        combined.append(new_line)
        glu_atom_count += 1
    elif line.startswith("CONECT"):
        break

combined.append("END\n")

with open("combined.pdb", "w") as f:
    f.writelines(combined)
print(f"✓ Created combined.pdb: {pva_atoms} PVA + {glu_atom_count} GLU = {pva_atoms + glu_atom_count} atoms")
EOF

# Step 4: Convert AMBER topologies to LAMMPS format
python3 amber_to_lammps.py system.lammps parm.lammps combined.pdb \
  -t PVA7_trim.top charge_data/glutaraldehyde.top \
  -c 1 1 \
  --charges 0 0 \
  --verbose

echo "✓ Step 4: LAMMPS files generated (system.lammps, parm.lammps)"
```

**Output Files:**
- `system.lammps` — LAMMPS data file (101 atoms: 70 PVA + 31 GLU)
- `parm.lammps` — LAMMPS parameter file with force field coefficients

See **Step 4** below for details on the LAMMPS conversion.

---

## Step 4: Convert AMBER Topology to LAMMPS Format

**File:** `amber_to_lammps.py`

**Function:** Command-line tool using ParmEd to convert AMBER `.prmtop` files to LAMMPS data/parameter files

**Description:**
Converts AMBER topology files (generated in Step 2) into LAMMPS-compatible data and parameter files. This is the final step needed to run simulations in LAMMPS.

**Inputs:**
- AMBER topology file(s) (`.prmtop`) from Step 2
- Combined PDB file with all molecules positioned (typically from Packmol)

**Outputs:**
- `system.lammps` (or custom name) — LAMMPS data file containing:
  - Atom count, bond count, angle count, dihedral count
  - Box dimensions and boundary conditions
  - Atomic coordinates and types
  - Bond, angle, dihedral lists
  
- `parm.lammps` (or custom name) — LAMMPS parameter file containing:
  - Force field coefficients (bond constants, angle constants, dihedral coefficients)
  - Nonbonded Lennard-Jones parameters (sigma, epsilon)
  - Atom mass definitions

**Command-Line Usage:**
```bash
python amber_to_lammps.py <output_data> <output_parm> <combined.pdb> \
  -t <topology.prmtop> \
  -c <molecule_count> \
  --charges <net_charge>
```

**Example (single molecule type):**
```bash
python amber_to_lammps.py system.lammps parm.lammps combined.pdb \
  -t PVA7_trim.prmtop \
  -c 1 \
  --charges 0
```

**Example (multiple molecule types, e.g., PVA + GLU):**
```bash
python amber_to_lammps.py pva_glu_system.lammps pva_glu_parm.lammps combined.pdb \
  -t PVA7_trim.prmtop glutaraldehyde.prmtop \
  -c 10 5 \
  --charges 0 0 \
  --verbose
```

**Key Options:**
- `-t / --topologies` — AMBER topology files (one or more)
- `-c / --counts` — Number of each molecule type (same order as topologies)
- `--charges` — Target net charge per molecule type (use 0 for neutral)
- `-b / --buffer` — Padding around molecules (default: 3.8 Å)
- `--verbose` — Print detailed progress messages
- `--keep-temp` — Retain temporary intermediate files (bonds.txt, angles.txt, etc.)

**Requirements:**
- ParmEd: `pip install parmed` or `conda install -c conda-forge parmed`
- NumPy

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
  
- **`amber_to_lammps.py`** (Step 4: AMBER → LAMMPS conversion)
  - Converts AMBER topology files (`.prmtop`) to LAMMPS data/parameter files
  - Requires: ParmEd, NumPy
  - Command-line tool; can be called from Python or shell
  
- **`genhydrogel.py`** (Full orchestration, optional)
  - Integrates all four steps plus additional steps (network generation, crosslink formation)
  - Depends on: all above modules

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

### Complete Workflow (PVA + GLU → LAMMPS)

```bash
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization
conda activate AmberTools25

# Step 1: Build PVA
python3 << 'EOF'
from pva_builder import build_pva
build_pva(n=7, output_file="PVA7_trim.pdb", cap=False)
print("✓ PVA7_trim.pdb generated")
EOF

# Step 2: Parametrize both PVA and GLU
python3 << 'EOF'
from molecular_utils import getfiles
getfiles("PVA7_trim.pdb")
print("✓ PVA7_trim.top generated")
getfiles("charge_data/glutaraldehyde.pdb")
print("✓ glutaraldehyde.top generated")
EOF

# Step 3: Load charges
python3 << 'EOF'
from molecular_utils import load_system_charges
charges = load_system_charges(chain_length=7, n_pva=1, n_glu=1)
print(f"✓ Charges loaded: {len(charges)} atoms")
EOF

# Create combined PDB (see full workflow above for details)
python3 << 'EOF'
# Combine PVA7_trim.pdb and charge_data/glutaraldehyde.pdb into combined.pdb
# (See "Complete Workflow" section for full code)
EOF

# Step 4: Convert to LAMMPS
python3 amber_to_lammps.py system.lammps parm.lammps combined.pdb \
  -t PVA7_trim.top charge_data/glutaraldehyde.top \
  -c 1 1 --charges 0 0

echo "✓ Complete! LAMMPS files ready: system.lammps, parm.lammps"
```

**Result:**
- `system.lammps`: 101 atoms (70 PVA + 31 GLU), 101 bonds, 185 angles, 281 dihedrals
- `parm.lammps`: Force field parameters for all interactions

---

## Complete PVA–GLU Combined Workflow (Large System)

For production simulations with many PVA chains and GLU crosslinkers (e.g., 300 PVA chains + 150 GLU molecules), the workflow is identical but with scaled molecule counts:

```bash
#!/bin/bash
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization
conda activate AmberTools25

# Step 1: Build reference PVA structure (n=7)
python3 << 'EOF'
from pva_builder import build_pva
build_pva(n=7, output_file="PVA7_trim.pdb", cap=False)
EOF

# Step 2: Parametrize reference molecules (one-time step)
python3 << 'EOF'
from molecular_utils import getfiles
getfiles("PVA7_trim.pdb")
getfiles("charge_data/glutaraldehyde.pdb")
EOF

# Step 3: Load charges for large system
python3 << 'EOF'
from molecular_utils import load_system_charges
# Load charges for 300 PVA chains + 150 GLU molecules (300*70 + 150*31 atoms)
charges = load_system_charges(chain_length=7, n_pva=300, n_glu=150)
total_atoms = len(charges)
print(f"✓ Loaded charges: {total_atoms} atoms ({300*70} PVA + {150*31} GLU)")
EOF

# Step 4: Create combined PDB with Packmol and convert to LAMMPS
# (Use Packmol or similar to combine 300 copies of PVA + 150 copies of GLU)
# Then convert:
python3 amber_to_lammps.py pva_glu_large.lammps pva_glu_large.parm combined_large.pdb \
  -t PVA7_trim.top charge_data/glutaraldehyde.top \
  -c 300 150 \
  --charges 0 0 \
  --verbose

echo "✓ Large system ready for MD simulation"
echo "  LAMMPS data: pva_glu_large.lammps"
echo "  LAMMPS parm: pva_glu_large.parm"
```

**Typical System Sizes:**
| n_pva | n_glu | Total Atoms | Approx. Bonds | Use Case |
|-------|-------|------------|---------------|----------|
| 1 | 1 | 101 | 101 | Testing / Quick validation |
| 10 | 5 | 950 | 950 | Small test system |
| 300 | 150 | 24,150 | 24,150 | Production simulations |

**Key Points for PVA–GLU Systems:**
1. **One-time parametrization**: Generate `PVA7_trim.top` and `charge_data/glutaraldehyde.top` once
2. **Reuse for all system sizes**: Scale only the molecule counts (`-c 300 150`), not the parametrization
3. **Charge loading scales automatically**: `load_system_charges(n_pva=300, n_glu=150)` returns all charges needed
4. **Combined PDB**: Use Packmol or similar tool to generate `combined_large.pdb` with all molecules positioned

**Note on Saturated vs. Unsaturated GLU:**
- The parametrized topology (`glutaraldehyde.top`) is derived from the **saturated reference structure** (`charge_data/glutaraldehyde.pdb`)
- The pre-extracted charges (`glutaraldehyde_charges.txt`) are for the **unsaturated form** used in actual network construction
- Both maintain 31 atoms; the difference is in which hydrogens are explicitly present (junction carbons have H in saturated, removed in unsaturated for reactivity)

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
   - Antechamber charges from hard-coded geometry may differ slightly from those from minimized structures
   - Pre-extracted charges from minimized reference (PVA7_min.pdb) provide chemical consistency
   - Consistent across all jobs and avoids computational redundancy
   - Ensures reproducibility of the workflow

---

## Contact & Attribution

Pipeline developed for PVA-GLU hydrogel MD simulations using Signac workflow management.
