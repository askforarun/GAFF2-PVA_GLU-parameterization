# Installation Instructions

## Quick Install

```bash
# Clone the repository
git clone https://github.com/askforarun/GAFF2-PVA-parameterization.git
cd GAFF2-PVA-parameterization

# Activate AmberTools environment (or create one)
conda activate AmberTools25
# OR: conda install -c conda-forge ambertools

# Run the example
python example_parametrization.py
```

---

## Detailed Setup

### 1. Clone Repository

```bash
git clone https://github.com/askforarun/GAFF2-PVA-parameterization.git
cd GAFF2-PVA-parameterization
```

### 2. Install AmberTools (Required)

**Option A: Using existing conda environment**
```bash
conda activate AmberTools25
# Verify installation
which antechamber
which parmchk2
which tleap
```

**Option B: Create new conda environment**
```bash
conda create -n gaff2-pva -c conda-forge ambertools
conda activate gaff2-pva
```

**Option C: System-wide installation (advanced)**
```bash
# If AmberTools is installed system-wide, ensure it's on PATH
export AMBERHOME=/path/to/amber
export PATH=$AMBERHOME/bin:$PATH
```

### 3. Verify Installation

```bash
python3 -c "
import numpy as np
print('✓ NumPy:', np.__version__)

# Test Antechamber availability
import subprocess
result = subprocess.run('antechamber -h', shell=True, capture_output=True)
if result.returncode == 0:
    print('✓ Antechamber: available')
else:
    print('✗ Antechamber: NOT available')
    
# Test parmchk2
result = subprocess.run('parmchk2 -h', shell=True, capture_output=True)
if result.returncode == 0:
    print('✓ Parmchk2: available')
else:
    print('✗ Parmchk2: NOT available')

# Test tleap
result = subprocess.run('tleap -s', shell=True, capture_output=True, input=b'quit\n')
if result.returncode == 0:
    print('✓ tLeap: available')
else:
    print('✗ tLeap: NOT available')
"
```

Expected output:
```
✓ NumPy: 1.xx.x
✓ Antechamber: available
✓ Parmchk2: available
✓ tLeap: available
```

---

## System-Specific Installation

### Linux (Ubuntu/Debian)

```bash
# Install Miniconda (if not already installed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create and activate environment
conda create -n gaff2-pva python=3.10 -c conda-forge
conda activate gaff2-pva

# Install AmberTools
conda install -c conda-forge ambertools

# Clone and test
git clone https://github.com/askforarun/GAFF2-PVA-parameterization.git
cd GAFF2-PVA-parameterization
python example_parametrization.py
```

### macOS

```bash
# Install Miniconda (if not already installed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
# or for Intel: Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-*.sh

# Create and activate environment
conda create -n gaff2-pva python=3.10 -c conda-forge
conda activate gaff2-pva

# Install AmberTools (may take longer on macOS)
conda install -c conda-forge ambertools

# Clone and test
git clone https://github.com/askforarun/GAFF2-PVA-parameterization.git
cd GAFF2-PVA-parameterization
python example_parametrization.py
```

### HPC Cluster (e.g., Alces Flight)

```bash
# Load necessary modules
flight env activate gridware
module load apps/gromacs_cuda/2026.1  # or appropriate module

# Activate AmberTools environment
conda activate AmberTools25

# Clone and run
git clone https://github.com/askforarun/GAFF2-PVA-parameterization.git
cd GAFF2-PVA-parameterization
python example_parametrization.py
```

---

## Troubleshooting Installation

### Problem: "antechamber: command not found"

**Cause:** AmberTools not installed or not on PATH

**Solutions:**

1. Install AmberTools:
   ```bash
   conda install -c conda-forge ambertools
   ```

2. Verify installation:
   ```bash
   which antechamber
   # Should return: /path/to/conda/envs/ENV_NAME/bin/antechamber
   ```

3. If installed but not found, add to PATH:
   ```bash
   export AMBERHOME=$CONDA_PREFIX
   export PATH=$AMBERHOME/bin:$PATH
   ```

### Problem: "ModuleNotFoundError: No module named 'numpy'"

**Cause:** NumPy not installed

**Solution:**
```bash
conda install numpy
# or: pip install numpy
```

### Problem: "leap.log" shows error from tLeap

**Cause:** Antechamber or tLeap issues with MOL2 file

**Solution:**

1. Check the leap.log file:
   ```bash
   cat leap.log
   ```

2. Verify MOL2 file is valid:
   ```bash
   head -20 PVA7_trim.mol2
   ```

3. Try manual parametrization:
   ```bash
   antechamber -j 4 -at gaff2 -dr no -fi pdb -fo mol2 \
     -i PVA7_trim.pdb -o PVA7_trim.mol2 -v
   ```

### Problem: Permission denied when cloning

**Cause:** SSH key not configured or repository access issue

**Solutions:**

1. Use HTTPS instead of SSH:
   ```bash
   git clone https://github.com/askforarun/GAFF2-PVA-parameterization.git
   ```

2. If SSH is required, set up SSH key:
   ```bash
   ssh-keygen -t ed25519 -C "your-email@example.com"
   # Add public key to GitHub: https://github.com/settings/keys
   ```

### Problem: "File not found: charge_data/PVA_monomercharges.txt"

**Cause:** Repository not fully cloned or files missing

**Solution:**

Verify all files are present:
```bash
ls -la charge_data/
# Should show: PVA_monomercharges.txt, glutaraldehyde_charges.txt, etc.

# If missing, re-clone:
cd ..
rm -rf GAFF2-PVA-parameterization
git clone https://github.com/askforarun/GAFF2-PVA-parameterization.git
cd GAFF2-PVA-parameterization
```

---

## Verify Complete Installation

Run the complete example to verify everything works:

```bash
python example_parametrization.py
```

Expected output:
```
======================================================================
PVA GAFF2 Parametrization Pipeline
======================================================================

Configuration:
  Chain length:    7 monomers per PVA chain
  PVA chains:      2
  GLU molecules:   1

──────────────────────────────────────────────────────────────────────
STEP 1: Build PVA with Hard-Coded Geometry
──────────────────────────────────────────────────────────────────────
✓ Successfully built PVA structure
  Output file: PVA7_trim.pdb
  Total atoms: 71
  Total bonds: 70

──────────────────────────────────────────────────────────────────────
STEP 2: Parametrize with GAFF2 (Antechamber + Parmchk2)
──────────────────────────────────────────────────────────────────────
Running Antechamber, Parmchk2, and tLeap...

✓ Successfully parametrized PVA
  Output files:
    - PVA7_trim.mol2       (Parametrized structure in MOL2 format)
    - PVA7_trim.frcmod     (Force field corrections)
    - PVA7_trim.top        (AMBER topology file)
    - PVA7_trim.crd        (AMBER coordinate file)

──────────────────────────────────────────────────────────────────────
STEP 3: Assign Pre-Extracted Partial Charges
──────────────────────────────────────────────────────────────────────
✓ Successfully loaded partial charges
  Total atoms: 228
  Breakdown:
    - PVA atoms: 140 (2 chains × 7 monomers × 10 atoms/monomer)
    - GLU atoms: 18 (1 molecules × 18 atoms/molecule)
  
  Net system charge: 0.00000 (should be ≈0.0)

──────────────────────────────────────────────────────────────────────
PIPELINE COMPLETE
──────────────────────────────────────────────────────────────────────

✓ Step 1: Built PVA structure (PVA7_trim.pdb)
✓ Step 2: Parametrized with GAFF2 (.mol2, .frcmod, .top, .crd)
✓ Step 3: Loaded 228 pre-extracted partial charges

The system is now ready for:
  - AMBER→LAMMPS conversion (amber_to_lammps.py)
  - LAMMPS compression and crosslinking simulations
  - Full hydrogel MD workflow (hydrogel_signac.py)

======================================================================
```

If you see this output with all ✓ marks, installation is **successful**!

---

## Optional: Full Hydrogel Workflow Integration

If you want to use this in the full hydrogel simulation workflow:

```bash
# Install additional dependencies for full workflow
conda install -c conda-forge signac signac-flow lammps gromacs

# Clone main hydrogel simulation repo
git clone /mnt/scratch/users/ass2009/hydrogel_simulation
cd hydrogel_simulation

# Initialize and run workflow
python hydrogel_signac.py init
python hydrogel_signac.py run -o generate_hydrogel
```

---

## Optional: Development Setup (for contributing)

If you want to modify or extend the code:

```bash
# Clone with development mode
git clone https://github.com/askforarun/GAFF2-PVA-parameterization.git
cd GAFF2-PVA-parameterization

# Install with development dependencies
conda install -c conda-forge ambertools pytest pytest-cov

# Run tests (if available)
pytest tests/

# Make your changes and test locally
python example_parametrization.py
```

---

## Getting Help

If installation fails:

1. **Check the logs:**
   ```bash
   cat leap.log  # For tLeap errors
   conda list    # Check installed packages
   which antechamber  # Verify AmberTools on PATH
   ```

2. **Review documentation:**
   - [README.md](README.md) — Overview
   - [SETUP_GUIDE.md](SETUP_GUIDE.md) — Troubleshooting
   - [SUMMARY.txt](SUMMARY.txt) — Quick reference

3. **Common issues:**
   - See "Troubleshooting Installation" section above
   - Check GitHub issues: https://github.com/askforarun/GAFF2-PVA-parameterization/issues

4. **Verify Python version:**
   ```bash
   python --version
   # Should be 3.7 or later
   ```

---

## System Requirements

**Minimum:**
- Python 3.7+
- NumPy
- AmberTools (with antechamber, parmchk2, tleap)
- 200 MB disk space
- 1 GB RAM

**Recommended:**
- Python 3.10+
- Conda/Miniconda
- 500 MB disk space (for reference data)
- 2 GB RAM

**Optional (for full workflow):**
- Signac + Flow (workflow management)
- LAMMPS (MD simulations)
- GROMACS (equilibration and analysis)

---

## Next Steps After Installation

1. **Run example:**
   ```bash
   python example_parametrization.py
   ```

2. **Read documentation:**
   ```bash
   cat SUMMARY.txt        # Quick overview
   cat README.md          # Complete reference
   ```

3. **Use in your own work:**
   ```python
   from pva_builder import build_pva
   from molecular_utils import getfiles, load_system_charges
   
   # Build PVA chain
   build_pva(n=10, output_file="PVA10_trim.pdb", cap=False)
   
   # Parametrize
   getfiles("PVA10_trim.pdb")
   
   # Load charges
   charges = load_system_charges(chain_length=10, n_pva=2, n_glu=1)
   ```

---

## Installation Verification Checklist

- [ ] Repository cloned successfully
- [ ] AmberTools installed (antechamber, parmchk2, tleap available)
- [ ] NumPy installed
- [ ] `example_parametrization.py` runs without errors
- [ ] Output files generated (.pdb, .mol2, .frcmod, .top, .crd)
- [ ] Charges loaded successfully (net charge ≈ 0.0)
- [ ] Documentation accessible (README.md, SETUP_GUIDE.md, etc.)

All checks passing? You're ready to use the toolkit! ✓

---

## Support & Contact

For issues or questions:
- GitHub: https://github.com/askforarun/GAFF2-PVA-parameterization
- Main project: https://github.com/askforarun/hydrogel_simulation
- Email: askforarun@gmail.com

---

**Installation Guide v1.0**  
**Last Updated:** May 12, 2026
