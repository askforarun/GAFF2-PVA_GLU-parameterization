# Parameterise Molecules: Complete Guide

## What is `parameterise_molecules()`?

**Function signature:**
```python
def parameterise_molecules(polymer: str, crosslinker: str) -> tuple[str, str]:
```

**Purpose:**
Parametrize both PVA polymer and GLU crosslinker molecules independently using GAFF2 force field and return base names of parametrized files.

**Returns:**
- `pva_mod`: Base name of corrected PVA files (without extension)
- `glu_mod`: Base name of corrected GLU files (without extension)

---

## Step-by-Step Process

### Input
```python
polymer = "PVA7_trim.pdb"      # Uncapped PVA structure from build_pva()
crosslinker = "/path/to/charge_data/glutaraldehyde.pdb"  # GLU reference
```

### Process

#### For PVA Polymer:

**Step 1: Get base name**
```python
pva_base = os.path.splitext(polymer)[0]  # "PVA7_trim"
```

**Step 2: Run Antechamber to generate MOL2 and initial parametrization**
```python
getfiles(polymer)  # Runs:
# antechamber -j 4 -at gaff2 -dr no -fi pdb -fo mol2 \
#   -i PVA7_trim.pdb -o PVA7_trim.mol2
# parmchk2 -i PVA7_trim.mol2 -o PVA7_trim.frcmod -f mol2 -a Y
# tleap script to generate .top and .crd
```

**Outputs from getfiles():**
- `PVA7_trim.mol2` — MOL2 format with GAFF2 atom types
- `PVA7_trim.frcmod` — Force field parameters
- `PVA7_trim.top` — AMBER topology
- `PVA7_trim.crd` — AMBER coordinates

**Step 3: Hand-correct atom types in MOL2 file**
```python
with open(f"{pva_base}.mol2") as f:
    content = f.read().replace("ha", "hc").replace("c2", "c3")
    # ha (aromatic H) → hc (alkane H)
    # c2 (sp2 carbon) → c3 (sp3 carbon)
```

**Why correct?**
- Antechamber may assign aromatic types to aliphatic molecules
- PVA is purely aliphatic (no aromatic rings)
- Correct types improve force field accuracy

**Step 4: Save corrected MOL2 with "_mod" suffix**
```python
pva_mod = f"{pva_base}_mod"  # "PVA7_trim_mod"
with open(f"{pva_mod}.mol2", "w") as f:
    f.write(content)
```

**Step 5: Regenerate FRCMOD with corrected atom types**
```python
subprocess.run(
    f"parmchk2 -i {pva_mod}.mol2 -o {pva_mod}.frcmod -f mol2 -a Y",
    shell=True,
)
# Generates: PVA7_trim_mod.frcmod with corrected parameters
```

**Step 6: Print result**
```python
print(f"PVA FF parameters → {pva_mod}.frcmod")
# Output: "PVA FF parameters → PVA7_trim_mod.frcmod"
```

**Return value for PVA:** `"PVA7_trim_mod"`

---

#### For GLU Crosslinker (Similar Process):

**Step 1: Copy GLU from reference path**
```python
glu_src = Path(crosslinker)  # "/path/to/glutaraldehyde.pdb"
glu_pdb = glu_src.name  # "glutaraldehyde.pdb"
shutil.copy2(glu_src, glu_pdb)  # Copy to current directory
```

**Step 2: Get base name**
```python
glu_base = os.path.splitext(glu_pdb)[0]  # "glutaraldehyde"
```

**Step 3: Run Antechamber**
```python
getfiles(glu_pdb)  # Same process as PVA
# Outputs: glutaraldehyde.mol2, glutaraldehyde.frcmod, etc.
```

**Step 4: Hand-correct atom types (different from PVA)**
```python
with open(f"{glu_base}.mol2") as f:
    content = f.read().replace("c2", "c6").replace("h4", "h1")
    # c2 (sp2 carbon) → c6 (aromatic sp2 carbon)
    # h4 (aromatic H) → h1 (aromatic H on H-bonded atom)
```

**Why different corrections?**
- GLU has carbonyl groups (C=O) with sp2 carbons
- Aromatic types are appropriate here
- Different atom type substitutions from PVA

**Step 5: Save corrected MOL2**
```python
glu_mod = f"{glu_base}_mod"  # "glutaraldehyde_mod"
with open(f"{glu_mod}.mol2", "w") as f:
    f.write(content)
```

**Step 6: Regenerate FRCMOD**
```python
subprocess.run(
    f"parmchk2 -i {glu_mod}.mol2 -o {glu_mod}.frcmod -f mol2 -a Y",
    shell=True,
)
# Generates: glutaraldehyde_mod.frcmod
```

**Step 7: Print result**
```python
print(f"GLU FF parameters → {glu_mod}.frcmod")
```

**Return value for GLU:** `"glutaraldehyde_mod"`

---

### Return Statement
```python
return pva_mod, glu_mod
# Returns: ("PVA7_trim_mod", "glutaraldehyde_mod")
```

---

## How to Use the Return Values

### After `parameterise_molecules()` completes:

```python
pva_mod, glu_mod = parameterise_molecules(polymer, crosslinker)
# pva_mod = "PVA7_trim_mod"
# glu_mod = "glutaraldehyde_mod"

print(pva_mod)  # "PVA7_trim_mod"
print(glu_mod)  # "glutaraldehyde_mod"
```

### Files Created (base names):

After parametrization, these files exist in the working directory:

**For PVA:**
```
PVA7_trim_mod.mol2         ← Corrected MOL2 with aliphatic types
PVA7_trim_mod.frcmod       ← Force field parameters for PVA
```

**For GLU:**
```
glutaraldehyde_mod.mol2    ← Corrected MOL2 with sp2 types
glutaraldehyde_mod.frcmod  ← Force field parameters for GLU
```

---

## Next Step: `convert_to_lammps()`

The return values are immediately used in the next function:

```python
pva_mod, glu_mod = parameterise_molecules(polymer, crosslinker)
convert_to_lammps(pva_mod, glu_mod, n_pva, n_glu)  # ← Pass the base names
```

### Inside `convert_to_lammps()`:

**For each molecule (PVA and GLU):**

```python
for mod, label in ((pva_mod, "PVA"), (glu_mod, "GLU")):
    # mod = "PVA7_trim_mod" or "glutaraldehyde_mod"
    # label = "PVA" or "GLU"
    
    # Step 1: Create tLeap input file
    tleap_in = f"tleap_{label.lower()}.in"  # "tleap_pva.in" or "tleap_glu.in"
    with open(tleap_in, "w") as f:
        f.write(f"source leaprc.gaff2\n")
        f.write(f"MOL = loadmol2 {mod}.mol2\n")           # Load PVA7_trim_mod.mol2
        f.write("check MOL\n")
        f.write(f"loadamberparams {mod}.frcmod\n")       # Load PVA7_trim_mod.frcmod
        f.write(f"saveamberparm MOL {mod}.top {mod}.crd\n")  # Save as PVA7_trim_mod.top
        f.write("quit")
    
    # Step 2: Run tLeap
    run_tleap_with_error_check(tleap_in)
    
    # Step 3: Print result
    print(f"{label} topology → {mod}.top")
    # Output: "PVA topology → PVA7_trim_mod.top"
    #         "GLU topology → glutaraldehyde_mod.top"
```

### Files Generated by `convert_to_lammps()`:

**From PVA:**
```
PVA7_trim_mod.top          ← AMBER topology for PVA
PVA7_trim_mod.crd          ← AMBER coordinates for PVA
```

**From GLU:**
```
glutaraldehyde_mod.top     ← AMBER topology for GLU
glutaraldehyde_mod.crd     ← AMBER coordinates for GLU
```

### Then `amber2lammps()` converts to LAMMPS:

```python
amber2lammps(
    data_file="data.lammps",
    param_file="parm.lammps",
    topologies=[f"{pva_mod}.top", f"{glu_mod}.top"],  # Pass the .top files
    molecule_counts=[n_pva, n_glu],
    pdb_file="packed_system.pdb",
    charges_target=[0, 0],
    verbose=True,
)
```

---

## Complete Workflow Chain

```
generatehydrogel()
    ↓
build_pva_system()
    ↓
polymer = "PVA7_trim.pdb"
crosslinker = "glutaraldehyde.pdb"
    ↓
parameterise_molecules(polymer, crosslinker)
    │
    ├─ For PVA:
    │   ├─ Run Antechamber on PVA7_trim.pdb
    │   ├─ Correct atom types (ha→hc, c2→c3)
    │   ├─ Run Parmchk2 on corrected MOL2
    │   └─ Output: pva_mod = "PVA7_trim_mod"
    │
    └─ For GLU:
        ├─ Run Antechamber on glutaraldehyde.pdb
        ├─ Correct atom types (c2→c6, h4→h1)
        ├─ Run Parmchk2 on corrected MOL2
        └─ Output: glu_mod = "glutaraldehyde_mod"
    ↓
return (pva_mod, glu_mod)
    ↓
convert_to_lammps(pva_mod, glu_mod, n_pva, n_glu)
    │
    ├─ For each module (PVA and GLU):
    │   ├─ Create tleap_{pva,glu}.in
    │   ├─ Load {pva_mod,glu_mod}.mol2
    │   ├─ Load {pva_mod,glu_mod}.frcmod
    │   └─ Save {pva_mod,glu_mod}.top and .crd
    │
    └─ Run amber2lammps():
        ├─ Read .top files
        ├─ Combine PVA and GLU topologies
        ├─ Create LAMMPS data.lammps
        └─ Create LAMMPS parm.lammps
```

---

## Key Points

### What `parameterise_molecules()` Returns:

| Item | Value | Type |
|------|-------|------|
| `pva_mod` | `"PVA7_trim_mod"` | str (base name) |
| `glu_mod` | `"glutaraldehyde_mod"` | str (base name) |

### What These Base Names Reference:

**Files that exist after the function:**
```
PVA7_trim_mod.mol2         (input to tLeap)
PVA7_trim_mod.frcmod       (input to tLeap)

glutaraldehyde_mod.mol2    (input to tLeap)
glutaraldehyde_mod.frcmod  (input to tLeap)
```

### Why Base Names?

Return only the base name (without extension) because:
- `.mol2` and `.frcmod` files are used as inputs to tLeap
- tLeap generates `.top` and `.crd` files with the same base name
- Same variable works for both input and output files
- Flexibility: can add any extension needed

---

## Atom Type Corrections Explained

### PVA Corrections: `ha → hc`, `c2 → c3`

**Why?**
- PVA is purely aliphatic (no aromatic rings)
- Antechamber may misidentify types as aromatic
- Force field accuracy improves with correct types

**Atom types:**
- `ha` = aromatic hydrogen (not needed for aliphatic)
- `hc` = alkane hydrogen ✓ (correct for PVA)
- `c2` = sp2 carbon (carbonyl/aromatic)
- `c3` = sp3 carbon ✓ (correct for PVA backbone)

### GLU Corrections: `c2 → c6`, `h4 → h1`

**Why?**
- GLU has carbonyl groups (C=O) with sp2 carbons
- Aromatic types are more appropriate for conjugated C=O
- Different chemistry from simple aliphatic PVA

**Atom types:**
- `c2` = generic sp2 carbon
- `c6` = aromatic sp2 carbon ✓ (better for C=O)
- `h4` = aromatic hydrogen on generic atom
- `h1` = aromatic hydrogen on H-bonded atom ✓ (better for carbonyl)

---

## Example Usage in Code

```python
# Step 1: Build structures
polymer, crosslinker = build_pva_system(n=7, n_pva=2, n_glu=1, 
                                        init_box_size=500, 
                                        packmol_seed=None)
# polymer = "PVA7_trim.pdb"
# crosslinker = "/path/to/glutaraldehyde.pdb"

# Step 2: Parametrize
pva_mod, glu_mod = parameterise_molecules(polymer, crosslinker)
# pva_mod = "PVA7_trim_mod"
# glu_mod = "glutaraldehyde_mod"

# Step 3: Convert to LAMMPS using the returned base names
convert_to_lammps(pva_mod, glu_mod, n_pva=2, n_glu=1)
# Uses: PVA7_trim_mod.mol2, PVA7_trim_mod.frcmod
#       glutaraldehyde_mod.mol2, glutaraldehyde_mod.frcmod
# Produces: PVA7_trim_mod.top, glutaraldehyde_mod.top
#           data.lammps, parm.lammps
```

---

## Summary

**What:** `parameterise_molecules()` parametrizes PVA and GLU independently

**Inputs:** PDB file paths (polymer, crosslinker)

**Process:**
1. Run Antechamber (atom typing, charge assignment)
2. Correct atom types for chemical accuracy
3. Run Parmchk2 (fill in missing parameters)
4. Return base names of corrected files

**Outputs:** Base names (`pva_mod`, `glu_mod`)

**Usage:** Pass to `convert_to_lammps()` which:
- Creates tLeap scripts
- Generates AMBER topologies (.top files)
- Converts to LAMMPS format

**Files created:** `.mol2` and `.frcmod` files for each molecule

**Key insight:** Function returns base names because same names are used for both input files (`.mol2`, `.frcmod`) and output files generated by tLeap (`.top`, `.crd`)

