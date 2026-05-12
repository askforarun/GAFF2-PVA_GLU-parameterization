# GAFF2-PVA-Parameterization Repository Index

## Quick Links

**Start Here:**
1. [README.md](README.md) — Complete overview of the three-step parametrization pipeline
2. [SUMMARY.txt](SUMMARY.txt) — High-level summary and quick reference
3. [SETUP_GUIDE.md](SETUP_GUIDE.md) — Installation, usage, and troubleshooting

**Detailed Documentation:**
- [PARAMETRIZATION_STEPS.md](PARAMETRIZATION_STEPS.md) — Visual flowcharts, detailed breakdown, quick reference

**Example & Testing:**
- [example_parametrization.py](example_parametrization.py) — Runnable example demonstrating all 3 steps

---

## Repository Contents

### Core Source Files (Ready to Use)

| File | Purpose | Size |
|------|---------|------|
| `pva_builder.py` | Step 1: Build PVA with hard-coded geometry | 19 KB |
| `molecular_utils.py` | Step 2 & 3: Parametrize with GAFF2, load charges | 21 KB |
| `system_constants.py` | Configuration constants | 323 B |
| `genhydrogel.py` | Full integration pipeline (optional) | 13 KB |
| `amber_to_lammps.py` | AMBER→LAMMPS conversion (optional) | 39 KB |

### Reference Data (charge_data/)

| File | Purpose | Type |
|------|---------|------|
| `PVA_monomercharges.txt` | Pre-extracted PVA partial charges | Text |
| `glutaraldehyde_charges.txt` | Pre-extracted GLU partial charges | Text |
| `PVA7_min.pdb` | Minimized PVA reference structure | PDB |
| `PVA7_min.mol2` | Parametrized PVA reference | MOL2 |
| `PVA7.pdb` | Uncapped PVA reference | PDB |
| `glutaraldehyde.pdb` | Minimized GLU reference | PDB |
| `crosslinked_struct_min.pdb` | Crosslink reference structure | PDB |
| `crosslinked_struct_min.mol2` | Parametrized crosslink reference | MOL2 |
| `extract_charges.py` | Utility to extract charges | Python |

### Documentation (4 Files)

| File | Audience | Length | Content |
|------|----------|--------|---------|
| `README.md` | Everyone | ~300 lines | Complete reference, file dependencies, architecture |
| `PARAMETRIZATION_STEPS.md` | Technical | ~400 lines | Visual flowcharts, detailed breakdown, examples |
| `SETUP_GUIDE.md` | Users | ~300 lines | Installation, usage, troubleshooting, next steps |
| `SUMMARY.txt` | Quick reference | ~250 lines | High-level overview, dependency map, key info |

### Example & Testing

| File | Purpose | Status |
|------|---------|--------|
| `example_parametrization.py` | Runnable example of all 3 steps | Ready to use |

---

## The Three Parametrization Steps

### Step 1: Build PVA with Hard-Coded Geometry
**Module:** `pva_builder.py`  
**Function:** `build_pva(n, output_file, cap=False)`  
**Input:** Chain length `n` (e.g., 7)  
**Output:** `PVA{n}_trim.pdb`

Build a PVA polymer using fixed bond lengths and tetrahedral angles.
Structure: `CH3-(CH2-CHOH-CH2)n-CH3` (cap=True) or `(CH2-CHOH-CH2)n` (cap=False)

See: [PARAMETRIZATION_STEPS.md](PARAMETRIZATION_STEPS.md#step-1-hard-coded-pva-geometry)

### Step 2: Parametrize with GAFF2
**Module:** `molecular_utils.py`  
**Function:** `getfiles(pdb_file)`  
**Input:** `PVA{n}_trim.pdb`  
**Outputs:** `.mol2`, `.frcmod`, `.top`, `.crd`

Run Antechamber (assign atom types & charges) + Parmchk2 (fill FF parameters) + tLeap (AMBER topology).

Requires: AmberTools

See: [PARAMETRIZATION_STEPS.md](PARAMETRIZATION_STEPS.md#step-2-gaff2-parametrization-via-antechamber)

### Step 3: Assign Pre-Extracted Charges
**Module:** `molecular_utils.py`  
**Function:** `load_system_charges(chain_length, n_pva, n_glu)`  
**Input:** Chain length, # of chains, # of GLU  
**Output:** numpy array of partial charges

Load and tile pre-extracted charges from minimized reference structures.

See: [PARAMETRIZATION_STEPS.md](PARAMETRIZATION_STEPS.md#step-3-pre-extracted-partial-charges)

---

## Getting Started

### 1. Installation
```bash
conda activate AmberTools25
# or: conda install -c conda-forge ambertools
```

### 2. Run Example
```bash
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization
python example_parametrization.py
```

### 3. Expected Output
```
✓ Step 1: Built PVA7_trim.pdb (71 atoms)
✓ Step 2: Generated PVA7_trim.mol2, .frcmod, .top, .crd
✓ Step 3: Loaded 228 pre-extracted charges
```

---

## Key Information

### Force Field
- **GAFF2** (General AMBER Force Field v2)
- Broad applicability to organic molecules
- Automated via Antechamber

### Files You'll Generate
After running the 3 steps:
- `PVA{n}_trim.pdb` — Uncapped PVA structure
- `PVA{n}_trim.mol2` — MOL2 format with GAFF2 types
- `PVA{n}_trim.frcmod` — Force field corrections
- `PVA{n}_trim.top` — AMBER topology
- Numpy array of partial charges

### Integration Points
- **Full hydrogel workflow:** Use output in `/mnt/scratch/users/ass2009/hydrogel_simulation/`
- **LAMMPS simulations:** Convert to LAMMPS format with `amber_to_lammps.py`
- **GROMACS:** Further conversion to GROMACS topology

---

## Directory Tree

```
GAFF2-PVA-parameterization/
├── pva_builder.py                    # Step 1: Build structure
├── molecular_utils.py                # Step 2 & 3: Parametrize & charges
├── system_constants.py               # Configuration
├── genhydrogel.py                    # Full pipeline (optional)
├── amber_to_lammps.py                # AMBER→LAMMPS (optional)
├── example_parametrization.py        # Runnable example
│
├── charge_data/
│   ├── PVA_monomercharges.txt        # Pre-extracted PVA charges
│   ├── glutaraldehyde_charges.txt    # Pre-extracted GLU charges
│   ├── PVA7_min.pdb                  # Reference minimized PVA
│   ├── PVA7_min.mol2                 # Reference parametrized PVA
│   ├── PVA7.pdb                      # Uncapped reference PVA
│   ├── glutaraldehyde.pdb            # Reference GLU
│   ├── crosslinked_struct_min.*      # Crosslink reference
│   └── extract_charges.py            # Charge extraction utility
│
├── README.md                         # Complete reference
├── PARAMETRIZATION_STEPS.md          # Detailed breakdown
├── SETUP_GUIDE.md                    # Installation & usage
├── SUMMARY.txt                       # High-level overview
├── INDEX.md                          # This file
│
└── .git/                             # Git repository
```

---

## Quick Reference

| Task | File | Function | Command |
|------|------|----------|---------|
| Build PVA | pva_builder.py | `build_pva(7)` | `python -c "from pva_builder import build_pva; build_pva(7)"` |
| Parametrize | molecular_utils.py | `getfiles()` | `python -c "from molecular_utils import getfiles; getfiles('PVA7_trim.pdb')"` |
| Load charges | molecular_utils.py | `load_system_charges()` | `python -c "from molecular_utils import load_system_charges; load_system_charges(7, 2, 1)"` |
| Run all steps | example_parametrization.py | `main()` | `python example_parametrization.py` |

---

## Troubleshooting

**Problem:** `antechamber not found`  
**Solution:** `conda install -c conda-forge ambertools`

**Problem:** `PVA_monomercharges.txt not found`  
**Solution:** Ensure current directory is `/users/ass2009/sharedscratch/GAFF2-PVA-parameterization`

**Problem:** Files not generated in Step 2  
**Solution:** Check `leap.log` for Antechamber/tLeap errors

See [SETUP_GUIDE.md](SETUP_GUIDE.md#troubleshooting) for more.

---

## Documentation Roadmap

**For a quick overview:**
→ Read: `SUMMARY.txt`

**For complete understanding:**
→ Read: `README.md`

**For detailed technical breakdown:**
→ Read: `PARAMETRIZATION_STEPS.md`

**For hands-on setup & usage:**
→ Read: `SETUP_GUIDE.md`

**For a working example:**
→ Run: `python example_parametrization.py`

---

## Next Steps

After parametrization:

1. **Use in hydrogel workflow:**
   ```bash
   cd /mnt/scratch/users/ass2009/hydrogel_simulation
   python hydrogel_signac.py init
   python hydrogel_signac.py run -o generate_hydrogel
   ```

2. **Convert to LAMMPS:**
   ```python
   from amber_to_lammps import amber2lammps
   ```

3. **Run MD simulations:**
   ```bash
   mpirun -np 64 lmp -in in.compression_lammps
   ```

---

## Repository Metadata

- **Location:** `/users/ass2009/sharedscratch/GAFF2-PVA-parameterization`
- **Size:** ~404 KB
- **Files:** 21 (5 source + 11 reference data + 4 docs + 1 example)
- **Status:** ✓ Complete and ready to use
- **Source:** Copied from `/mnt/scratch/users/ass2009/hydrogel_simulation`

---

## License & Attribution

Developed for PVA-GLU hydrogel MD simulations using Signac workflow management.

For full context:
- Main repo: `/mnt/scratch/users/ass2009/hydrogel_simulation/`
- CLAUDE.md: Project documentation
- README_signac.md: Full workflow guide

---

## Quick Command Reference

```bash
# Activate environment
conda activate AmberTools25

# Navigate to repo
cd /users/ass2009/sharedscratch/GAFF2-PVA-parameterization

# Run example (all 3 steps)
python example_parametrization.py

# Or run steps manually:

# Step 1: Build
python -c "from pva_builder import build_pva; build_pva(7, 'PVA7_trim.pdb', cap=False)"

# Step 2: Parametrize
python -c "from molecular_utils import getfiles; getfiles('PVA7_trim.pdb')"

# Step 3: Load charges
python -c "from molecular_utils import load_system_charges; charges = load_system_charges(7, 1, 1); print(f'Loaded {len(charges)} charges')"
```

---

**Last updated:** 2026-05-12

