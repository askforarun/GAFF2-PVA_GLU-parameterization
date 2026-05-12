# Three-Step PVA Parametrization Pipeline

## Visual Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                    STEP 1: BUILD PVA STRUCTURE                      │
│                                                                      │
│  Module: pva_builder.py                                             │
│  Function: build_pva(n, output_file, cap=False)                    │
│                                                                      │
│  Input:  Chain length n (e.g., n=7 for 7 monomers)                 │
│  Process: Hard-coded geometry using fixed bond lengths & angles     │
│  Output: PVA{n}_trim.pdb                                            │
│                                                                      │
│  Example:                                                            │
│    from pva_builder import build_pva                               │
│    build_pva(n=7, output_file="PVA7_trim.pdb", cap=False)          │
│                                                                      │
│  ✗ Not minimized (will be done in LAMMPS later)                    │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│          STEP 2: PARAMETRIZE WITH GAFF2 (Antechamber)              │
│                                                                      │
│  Module: molecular_utils.py                                         │
│  Function: getfiles(pdb_file)                                       │
│                                                                      │
│  Input:  PVA7_trim.pdb                                              │
│  Process:                                                            │
│    1. Antechamber assigns GAFF2 atom types & initial charges        │
│       antechamber -at gaff2 -fi pdb -fo mol2 ...                    │
│    2. Parmchk2 fills in missing force field parameters              │
│       parmchk2 -i PVA7_trim.mol2 -o PVA7_trim.frcmod ...            │
│    3. tLeap builds AMBER topology                                   │
│       saveamberparm SUS PVA7_trim.top PVA7_trim.crd                 │
│                                                                      │
│  Outputs:                                                            │
│    - PVA7_trim.mol2        (MOL2 parametrized structure)            │
│    - PVA7_trim.frcmod      (Force field corrections)                │
│    - PVA7_trim.top         (AMBER topology)                         │
│    - PVA7_trim.crd         (AMBER coordinates)                      │
│                                                                      │
│  ✓ Generates complete GAFF2 parameters (bonds, angles, dihedrals)  │
│  ⚠ Uses Antechamber-derived charges (not minimized geometry)        │
│  ℹ Hand-corrects atom types: ha→hc, c2→c3                          │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│        STEP 3: ASSIGN PRE-EXTRACTED PARTIAL CHARGES                 │
│                                                                      │
│  Module: molecular_utils.py                                         │
│  Function: load_system_charges(chain_length, n_pva, n_glu)          │
│                                                                      │
│  Input:  Chain length, number of chains, number of GLU molecules    │
│  Reference Files:                                                    │
│    - charge_data/PVA_monomercharges.txt      (from PVA7_min.pdb)   │
│    - charge_data/glutaraldehyde_charges.txt  (from GLU ref)         │
│                                                                      │
│  Process:                                                            │
│    1. Load charges for 1 PVA monomer (10 atoms)                     │
│    2. Tile charges across all n_pva chains                          │
│    3. Load charges for all n_glu molecules                          │
│    4. Return concatenated 1-D numpy array                           │
│                                                                      │
│  Output: charges = [pva1_charges, pva2_charges, ..., glu_charges]  │
│                                                                      │
│  ✓ Uses charges from minimized reference (PVA7_min.pdb)            │
│  ✓ Consistent across all jobs in workflow                          │
│  ⚠ Ignores Antechamber charges from Step 2                         │
│  ℹ Reason: Minimized reference → more chemically consistent        │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
        ┌─────────────────────────────────────┐
        │ Ready for AMBER→LAMMPS Conversion   │
        │ (genhydrogel.py continues from here)│
        └─────────────────────────────────────┘
```

---

## Detailed Breakdown

### STEP 1: Hard-Coded PVA Geometry

**File:** `pva_builder.py`

```python
def build_pva(n: int, output_file: str = None, cap: bool = True):
    """
    Build PVA polymer: CH3-(CH2-CHOH-CH2)n-CH3
    
    Args:
        n: Number of repeat units (e.g., 7)
        output_file: Output PDB filename (default: PVA{n}.pdb or PVA{n}_trim.pdb)
        cap: If False, removes terminal CH3 groups (produces uncapped structure)
    
    Returns:
        tuple: (atoms list, bonds dict)
    """
```

**Geometry Details:**
- Fixed bond lengths: C-C ≈ 1.54 Å, C-O ≈ 1.42 Å, C-H ≈ 1.09 Å
- Fixed angles: sp³ tetrahedral ≈ 109.47°
- No external minimization (OpenBabel/LAMMPS, etc.)

**Structure:**
```
CH3-(CH2-CHOH-CH2)-(CH2-CHOH-CH2)-...-CH3
 ↓
Repeat unit: C-H-H-C-H-O-H-C-H-H (per 3n+2 carbon backbone)
Terminal:    CH3 at each end
```

**Example Output (PVA7_trim.pdb):**
```
ATOM      1  C1  PVAA    1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  H1  PVAA    1       0.890   0.000   0.000  1.00  0.00           H
ATOM      3  H2  PVAA    1      -0.445   0.816   0.000  1.00  0.00           H
...
CONECT    1    2    3    4
CONECT    2    1
CONECT    3    1
CONECT    4    1    5    6
...
```

---

### STEP 2: GAFF2 Parametrization via Antechamber

**File:** `molecular_utils.py` → `getfiles()` function

**Process Flow:**
```
PVA7_trim.pdb
    ↓
[Antechamber]  ← Assigns GAFF2 atom types, initial charges
    ↓
PVA7_trim.mol2
    ↓
[Parmchk2]  ← Fills in missing FF parameters (bonds, angles, dihedrals)
    ↓
PVA7_trim.frcmod  (force field corrections)
    ↓
[tLeap]  ← Generates AMBER topology
    ↓
PVA7_trim.top, PVA7_trim.crd
```

**Antechamber Command:**
```bash
antechamber -j 4 -at gaff2 -dr no \
  -fi pdb -fo mol2 \
  -i PVA7_trim.pdb -o PVA7_trim.mol2
```

**Parmchk2 Command:**
```bash
parmchk2 -i PVA7_trim.mol2 -o PVA7_trim.frcmod -f mol2 -a Y
```

**tLeap Script:**
```
source leaprc.gaff2
SUS = loadmol2 PVA7_trim.mol2
check SUS
loadamberparams PVA7_trim.frcmod
saveamberparm SUS PVA7_trim.top PVA7_trim.crd
quit
```

**Hand-Corrections Applied:**
```python
# From parameterise_molecules() in genhydrogel.py
content = content.replace("ha", "hc").replace("c2", "c3")
# ha (aromatic H) → hc (alkane H)
# c2 (sp2 carbon) → c3 (sp3 carbon)
# Ensures correct GAFF2 atom types for aliphatic PVA
```

**Example Output (PVA7_trim.frcmod):**
```
bond parameters
c3 c3   222.0   1.5380
c3 oh   320.0   1.4100
oh ho   553.0   0.9600
...

angle parameters
c3 c3 c3    63.0   109.47
c3 c3 oh    67.0   109.47
...

dihedral parameters
X  c3 c3  X      1    0.1556 180.0 -2.0
...
```

---

### STEP 3: Pre-Extracted Partial Charges

**File:** `molecular_utils.py` → `load_system_charges()` function

**Reference Data:**
```
charge_data/
├── PVA_monomercharges.txt     ← Charges for 1 PVA monomer (10 atoms)
├── glutaraldehyde_charges.txt ← Charges for GLU (1 molecule)
├── PVA7_min.pdb               ← Source: minimized PVA reference
└── glutaraldehyde.pdb         ← Source: minimized GLU reference
```

**Format of Charge File (PVA_monomercharges.txt):**
```
Atom  Type  OrigCharge  CorrectedCharge
C1    c3    0.098       0.100
H1    hc   -0.023      -0.025
H2    hc   -0.023      -0.025
C2    c3   -0.052      -0.050
H3    hc    0.045       0.045
O1    oh   -0.532      -0.535
HO    ho    0.410       0.410
C3    c3   -0.052      -0.050
H4    hc    0.045       0.045
H5    hc    0.045       0.045
```

**Charge Array Construction:**
```python
charges = load_system_charges(chain_length=7, n_pva=3, n_glu=1)

# Returns:
# [ pva_monomer_charges (10 atoms) × 7 monomers,  ← for PVA chain 1
#   pva_monomer_charges (10 atoms) × 7 monomers,  ← for PVA chain 2
#   pva_monomer_charges (10 atoms) × 7 monomers,  ← for PVA chain 3
#   glutaraldehyde_charges (18 atoms) ]            ← for GLU molecule
# Total: (10×7×3) + 18 = 228 atoms
```

**Why Pre-Extracted?**
1. **Consistency**: All jobs use same reference (PVA7_min.pdb)
2. **Chemistry**: Minimized structure → more realistic charges
3. **Speed**: No re-running Antechamber per job
4. **Validation**: Pre-checked for charge neutrality and correctness

---

## Integration with Full Workflow

The three steps fit into the complete `generatehydrogel()` pipeline:

```python
def generatehydrogel(n, n_pva, n_glu, init_box_size):
    cleanup_beginning()
    
    # Step 3: Load charges
    charges = load_system_charges(n, n_pva, n_glu)
    
    # Step 1: Build structure
    polymer, crosslinker = build_pva_system(n, n_pva, n_glu, init_box_size)
    
    # Step 2: Parametrize
    pva_mod, glu_mod = parameterise_molecules(polymer, crosslinker)
    
    # Additional steps...
    convert_to_lammps(pva_mod, glu_mod, n_pva, n_glu)
    cross_struct = build_crosslink_reference()
    prepare_lammps_params(cross_struct, charges)  # ← Uses charges from Step 3
    
    cleanup_end()
```

---

## Quick Reference

| Step | Input | Command/Function | Output | Requires |
|------|-------|------------------|--------|----------|
| **1** | `n` | `build_pva(n, cap=False)` | `PVA{n}_trim.pdb` | Python, numpy |
| **2** | `PVA{n}_trim.pdb` | `getfiles()` | `.mol2`, `.frcmod`, `.top`, `.crd` | AmberTools |
| **3** | `n, n_pva, n_glu` | `load_system_charges()` | Numpy array | charge_data/*.txt |

---

## Running the Three Steps Independently

```bash
#!/bin/bash
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization
conda activate AmberTools25

# Step 1: Build PVA
python3 -c "
from pva_builder import build_pva
build_pva(n=7, output_file='PVA7_trim.pdb', cap=False)
print('✓ Step 1 complete: PVA7_trim.pdb')
"

# Step 2: Parametrize
python3 -c "
from molecular_utils import getfiles
getfiles('PVA7_trim.pdb')
print('✓ Step 2 complete: PVA7_trim.mol2, .frcmod, .top, .crd')
"

# Step 3: Load charges
python3 -c "
from molecular_utils import load_system_charges
charges = load_system_charges(chain_length=7, n_pva=1, n_glu=1)
print(f'✓ Step 3 complete: {len(charges)} charges loaded')
"
```

---

## Notes

- **Step 1 geometry is NOT minimized**: Hard-coded scaffold; relaxation happens later in LAMMPS
- **Step 2 generates parameters but uses Antechamber charges**: These may differ from pre-minimized reference
- **Step 3 replaces charges with pre-extracted values**: Ensures consistency across workflow
- **All three steps are required** for a complete, chemically consistent parametrized PVA system
