#!/usr/bin/env python3
"""
Example script demonstrating the three-step PVA parametrization pipeline:
  Step 1: Build PVA with hard-coded geometry
  Step 2: Parametrize with GAFF2 (Antechamber + Parmchk2)
  Step 3: Assign pre-extracted partial charges
"""

import os
import sys
from pathlib import Path

# Add current directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from pva_builder import build_pva
from molecular_utils import getfiles, load_system_charges


def main():
    """Run the three-step parametrization pipeline."""

    print("=" * 70)
    print("PVA GAFF2 Parametrization Pipeline")
    print("=" * 70)

    # Configuration
    chain_length = 7          # Number of PVA monomers per chain
    n_pva = 2                 # Number of PVA chains
    n_glu = 1                 # Number of GLU molecules

    print(f"\nConfiguration:")
    print(f"  Chain length:    {chain_length} monomers per PVA chain")
    print(f"  PVA chains:      {n_pva}")
    print(f"  GLU molecules:   {n_glu}")
    print()

    # =========================================================================
    # STEP 1: Build PVA structure with hard-coded geometry
    # =========================================================================
    print("─" * 70)
    print("STEP 1: Build PVA with Hard-Coded Geometry")
    print("─" * 70)
    print(f"Module: pva_builder.py")
    print(f"Function: build_pva(n={chain_length}, cap=False)")
    print()

    try:
        output_file = f"PVA{chain_length}_trim.pdb"
        atoms, bonds = build_pva(
            n=chain_length,
            output_file=output_file,
            cap=False
        )
        print(f"✓ Successfully built PVA structure")
        print(f"  Output file: {output_file}")
        print(f"  Total atoms: {len(atoms)}")
        print(f"  Total bonds: {sum(len(b) for b in bonds.values()) // 2}")
        print()
        print(f"  Structure: CH3-(CH2-CHOH-CH2)×{chain_length}")
        print(f"  Geometry: Hard-coded (fixed bond lengths & angles)")
        print(f"  Status: NOT minimized (scaffold for parametrization)")
        print()
    except Exception as e:
        print(f"✗ Error in Step 1: {e}")
        return False

    # =========================================================================
    # STEP 2: Parametrize with GAFF2
    # =========================================================================
    print("─" * 70)
    print("STEP 2: Parametrize with GAFF2 (Antechamber + Parmchk2)")
    print("─" * 70)
    print(f"Module: molecular_utils.py")
    print(f"Function: getfiles('{output_file}')")
    print()
    print("Process:")
    print(f"  1. Antechamber → Assigns GAFF2 atom types and initial charges")
    print(f"  2. Parmchk2 → Fills in missing force field parameters")
    print(f"  3. tLeap → Generates AMBER topology")
    print()

    try:
        print(f"Running Antechamber, Parmchk2, and tLeap...")
        getfiles(output_file)

        base_name = output_file.replace(".pdb", "")
        print()
        print(f"✓ Successfully parametrized PVA")
        print(f"  Output files:")
        print(f"    - {base_name}.mol2       (Parametrized structure in MOL2 format)")
        print(f"    - {base_name}.frcmod     (Force field corrections)")
        print(f"    - {base_name}.top        (AMBER topology file)")
        print(f"    - {base_name}.crd        (AMBER coordinate file)")
        print()
        print(f"  Force field: GAFF2 (General AMBER Force Field v2)")
        print(f"  Status: Complete parameters (bonds, angles, dihedrals)")
        print()
    except Exception as e:
        print(f"✗ Error in Step 2: {e}")
        print(f"  Note: Requires AmberTools (Antechamber, Parmchk2, tLeap)")
        print(f"  Install: conda install -c conda-forge ambertools")
        return False

    # =========================================================================
    # STEP 3: Load pre-extracted partial charges
    # =========================================================================
    print("─" * 70)
    print("STEP 3: Assign Pre-Extracted Partial Charges")
    print("─" * 70)
    print(f"Module: molecular_utils.py")
    print(f"Function: load_system_charges({chain_length}, {n_pva}, {n_glu})")
    print()
    print("Reference Data:")
    print(f"  - charge_data/PVA_monomercharges.txt      (from PVA7_min.pdb)")
    print(f"  - charge_data/glutaraldehyde_charges.txt  (from GLU reference)")
    print()

    try:
        charges = load_system_charges(
            chain_length=chain_length,
            n_pva=n_pva,
            n_glu=n_glu
        )

        print(f"✓ Successfully loaded partial charges")
        print(f"  Total atoms: {len(charges)}")
        print(f"  Breakdown:")

        # Calculate component sizes
        pva_atoms_per_monomer = 10  # From system_constants
        glu_atoms_per_molecule = 18  # From charge_data

        pva_total_atoms = chain_length * pva_atoms_per_monomer * n_pva
        glu_total_atoms = n_glu * glu_atoms_per_molecule

        print(f"    - PVA atoms: {pva_total_atoms} ({n_pva} chains × {chain_length} monomers × {pva_atoms_per_monomer} atoms/monomer)")
        print(f"    - GLU atoms: {glu_total_atoms} ({n_glu} molecules × {glu_atoms_per_molecule} atoms/molecule)")
        print()

        # Show sample charges
        print(f"  Sample charges (first 5 atoms, 3 PVA atoms, 2 GLU atoms):")
        for i, charge in enumerate(charges[:5]):
            if i < 3:
                print(f"    Atom {i+1:2d} (PVA): {charge:8.5f}")
            else:
                print(f"    Atom {i+1:2d} (GLU): {charge:8.5f}")

        # Calculate net charge
        net_charge = sum(charges)
        print(f"  Net system charge: {net_charge:.5f} (should be ≈0.0)")
        print()

        print(f"  Charge source: Pre-extracted from minimized reference (PVA7_min.pdb)")
        print(f"  Status: Consistent across all jobs in workflow")
        print()

    except Exception as e:
        print(f"✗ Error in Step 3: {e}")
        return False

    # =========================================================================
    # Summary
    # =========================================================================
    print("─" * 70)
    print("PIPELINE COMPLETE")
    print("─" * 70)
    print()
    print(f"✓ Step 1: Built PVA structure ({output_file})")
    print(f"✓ Step 2: Parametrized with GAFF2 (.mol2, .frcmod, .top, .crd)")
    print(f"✓ Step 3: Loaded {len(charges)} pre-extracted partial charges")
    print()
    print("The system is now ready for:")
    print(f"  - AMBER→LAMMPS conversion (amber_to_lammps.py)")
    print(f"  - LAMMPS compression and crosslinking simulations")
    print(f"  - Full hydrogel MD workflow (hydrogel_signac.py)")
    print()
    print("=" * 70)

    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
