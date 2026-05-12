#!/usr/bin/env python3
"""
Extract charges from crosslinked structure and PVA monomer
Steps:
1) Use SMILES for crosslinked structure: CCC(O)CC1CC(CC(O)CC)OC(CCCC2OC(CC(O)CC)CC(CC(O)CC)O2)O1
2) Convert crosslinked structure to PDB using obabel
3) Minimize crosslinked structure using obabel
4) Run antechamber to produce crosslinked_struct_min.mol2 file
5) Extract PVA terminal group charges from crosslinked structure, apply correction, print to file
6) Extract glutaraldehyde charges from crosslinked structure, apply correction, print to file
7) Build PVA chain of length 7 using pva_builder.py
8) Minimize PVA chain and run antechamber
9) Extract PVA monomer charges for specified atoms, print to file
"""

import MDAnalysis as mda
import numpy as np
import subprocess
import os
import sys
import glob

# GUARD: Ensure this script can only run from the charge_data folder
current_dir = os.path.basename(os.getcwd())
if current_dir != 'charge_data':
    print("ERROR: This script must be run from the 'charge_data' folder.")
    print(f"Current directory: '{current_dir}'")
    print("Please navigate to the charge_data folder and run this script again.")
    sys.exit(1)

print("✓ Running from charge_data folder - OK")

# Import cleanup function from main hydrogel package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.file_operations import cleanup_beginning, cleanup_end
from src.molecular_utils import generate_monomer_from_smiles
from src.pva_builder import build_pva
from src.system_constants import GLU_ATOMS_PER_MOLECULE, PVA_ATOMS_PER_MONOMER


def _validate_atom_count(label: str, found: int, expected: int) -> None:
    """Ensure extracted atom sets stay aligned with shared system constants."""
    if found != expected:
        raise RuntimeError(
            f"{label} atom count mismatch: expected {expected}, found {found}. "
            "Update the selection list or the shared system constants."
        )

def step1_convert_smiles_to_pdb(smiles: str) -> str:
    """Step 1-2: Convert SMILES to PDB using obabel"""
    print("=== STEP 1-2: Converting SMILES to PDB ===")
    print(f"SMILES: {smiles}")
    
    pdb_file = 'crosslinked_struct.pdb'
    generated_file = generate_monomer_from_smiles(smiles, 'crosslinked_struct')
    
    print(f"Generated {pdb_file} from SMILES")
    return pdb_file

def step3_minimize_structure(pdb_file: str) -> str:
    """Step 3: Minimize using obabel"""
    print("\n=== STEP 3: Minimizing structure ===")
    
    base_name = os.path.splitext(pdb_file)[0]
    minimized_pdb = f'{base_name}_min.pdb'
    
    cmd = f'obabel {pdb_file} -opdb -O {minimized_pdb} --minimize --ff mmff94'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Obabel minimization failed: {result.stderr}")
    
    print(f"Minimized structure saved as {minimized_pdb}")
    return minimized_pdb

def step4_run_antechamber(pdb_file: str) -> str:
    """Step 4: Run antechamber to produce crosslinked_struct.mol2 file"""
    print("\n=== STEP 4: Running antechamber ===")
    
    mol2_file = 'crosslinked_struct_min.mol2'
    
    # Run antechamber
    cmd = f"antechamber -j 4 -at gaff2 -dr no -fi pdb -fo mol2 -i {pdb_file} -o {mol2_file} -c bcc"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Antechamber failed: {result.stderr}")
    
    print(f"Generated {mol2_file} with antechamber charges")
    
    # Check charges
    mol2 = mda.Universe(mol2_file)
    print(f"Sum of charges in {mol2_file}: {np.sum(mol2.atoms.charges)}")
    
    return mol2_file

def step5_extract_pva_terminal_group_charges(mol2_file: str):
    """Step 5: Extract PVA terminal group charges, apply correction, print to file"""
    print("\n=== STEP 5: Extracting PVA terminal group charges ===")
    
    # Load the mol2 file
    mol2 = mda.Universe(mol2_file)
    mol2.guess_TopologyAttrs(context='default', to_guess=['elements'])
    
    # Define PVA terminal group atoms in CH2-CHOH-CH2 order
    # These atoms come from crosslinked_struct_min.mol2 (minimized crosslinked PVA-glutaraldehyde structure)
    # This represents the terminal portion: CH3-(CH2-CHOH-CH2)n-CH3
    # CH2 GROUP 1 (first in order): C3, H7, H8
    # CHOH GROUP (middle): C2, H5, O, H6  
    # CH2 GROUP 2 (last in order): C1, H3, H4
    ch2_1_atoms = ['C3', 'H7', 'H8']
    choh_atoms = ['C2', 'H5', 'O', 'H6']
    ch2_2_atoms = ['C1', 'H3', 'H4']
    
    print("PVA_TERMINAL_GROUP CHARGES")
    print("="*50)
    print(f"{'Atom':<8} {'Charge':>12} {'Corrected':>12}")
    print("-"*50)
    
    pva_charges = []
    pva_corrected_charges = []
    
    # Select all PVA atoms
    all_pva_atoms = ch2_1_atoms + choh_atoms + ch2_2_atoms
    _validate_atom_count("PVA terminal group", len(all_pva_atoms), PVA_ATOMS_PER_MONOMER)
    pva_selection = " or ".join([f"name {name}" for name in all_pva_atoms])
    pva_mol_atoms = mol2.select_atoms(pva_selection)
    
    # Validate that all expected atoms are present to avoid IndexError
    missing_atoms = [name for name in all_pva_atoms if len(pva_mol_atoms[pva_mol_atoms.names == name]) == 0]
    if missing_atoms:
        raise RuntimeError(
            f"Expected PVA terminal atoms not found in {mol2_file}: {', '.join(missing_atoms)}. "
            "Ensure antechamber retained the atom names or update the selection list."
        )
    
    # Calculate residual charge and apply to C3
    total_pva_sum = np.sum([atom.charge for atom in pva_mol_atoms])
    residual_charge = total_pva_sum
    c3_correction = -residual_charge
    
    # Display in CH2-CHOH-CH2 order
    for group_name, group_atoms in [("CH2 GROUP 1", ch2_1_atoms), ("CHOH GROUP", choh_atoms), ("CH2 GROUP 2", ch2_2_atoms)]:
        if group_name == "CH2 GROUP 1":
            print(f"\n{group_name}:")
        elif group_name == "CHOH GROUP":
            print(f"\n{group_name}:")
        else:
            print(f"\n{group_name}:")
        
        for atom_name in group_atoms:
            atom = pva_mol_atoms[pva_mol_atoms.names == atom_name][0]
            original_charge = atom.charge
            corrected_charge = original_charge
            
            if atom.name == 'C3':
                corrected_charge += c3_correction
            
            print(f"{atom.name:<8} {original_charge:12.6f} {corrected_charge:12.6f}{'*' if atom.name == 'C3' else ' '}")
            pva_charges.append(original_charge)
            pva_corrected_charges.append(corrected_charge)
    
    print("\n" + "-"*50)
    print(f"{'SUM':<8} {np.sum(pva_charges):12.6f} {np.sum(pva_corrected_charges):12.6f}")
    print("* = corrected charge")
    print("="*50 + "\n")
    
    # Write to file in CH2-CHOH-CH2 order
    with open('PVA_terminal_group_charges.txt', 'w') as f:
        f.write("PVA_TERMINAL_GROUP CHARGES\n")
        f.write("="*50 + "\n")
        f.write(f"{'Atom':<8} {'Type':<8} {'Charge':>12} {'Corrected':>12}\n")
        f.write("-"*60 + "\n")
        
        for group_name, group_atoms in [("CH2 GROUP 1", ch2_1_atoms), ("CHOH GROUP", choh_atoms), ("CH2 GROUP 2", ch2_2_atoms)]:
            f.write(f"\n{group_name}:\n")
            for atom_name in group_atoms:
                atom = pva_mol_atoms[pva_mol_atoms.names == atom_name][0]
                original_charge = pva_charges[all_pva_atoms.index(atom_name)]
                corrected_charge = pva_corrected_charges[all_pva_atoms.index(atom_name)]
                star = '*' if atom.name == 'C3' else ' '
                f.write(f"{atom.name:<8} {atom.type:<8} {original_charge:12.6f} {corrected_charge:12.6f}{star}\n")
        
        f.write("\n" + "-"*60 + "\n")
        f.write(f"{'SUM':<8} {'':<8} {np.sum(pva_charges):12.6f} {np.sum(pva_corrected_charges):12.6f}\n")
        f.write("* = corrected charge\n")
        f.write("="*50 + "\n")
    
    print("PVA_terminal_group_charges.txt file generated")

def step6_extract_glutaraldehyde_charges(mol2_file: str):
    """Step 6: Extract glutaraldehyde charges, apply correction, print to file"""
    print("\n=== STEP 6: Extracting glutaraldehyde charges ===")
    
    # Load the mol2 file
    mol2 = mda.Universe(mol2_file)
    
    # Select glutaraldehyde atoms
    glutaraldehyde_atoms = mol2.select_atoms("name H10 or name H11 or name H12 or name C5 or name C6 or name O2 or name H9 or name C4 or name O7 or name C11 or name H22 or name H23 or name C12 or name C13 or name H24 or name H25 or name H26 or name H27 or name C14 or name H28 or name C15 or name H29 or name O6 or name O3 or name C16 or name H30 or name C22 or name C21 or name H40 or name H42 or name H41")
    _validate_atom_count("GLU molecule", len(glutaraldehyde_atoms), GLU_ATOMS_PER_MOLECULE)
    glutaraldehyde_atoms.write("glutaraldehyde.pdb")
    print(f"Selected {len(glutaraldehyde_atoms)} glutaraldehyde atoms")
    
    # Get charges
    raw_glutaraldehyde_charges = glutaraldehyde_atoms.charges
    
    # Adjust carbon charges to make net charge zero
    adjusted_charges = raw_glutaraldehyde_charges.copy()
    charge_adjustments = {}
    
    # Calculate charge adjustment needed: net_charge/4
    net_charge = np.sum(raw_glutaraldehyde_charges)
    charge_adjustment = -net_charge / 4
    
    print(f"Net charge: {net_charge:.6f}")
    
    if net_charge > 0:
        print(f"Net charge is positive: subtracting {abs(charge_adjustment):.6f} from each carbon")
    elif net_charge < 0:
        print(f"Net charge is negative: adding {abs(charge_adjustment):.6f} to each carbon")
    else:
        print("Net charge is already zero: no adjustment needed")
    
    # Adjust specific carbon charges
    for i, atom in enumerate(glutaraldehyde_atoms):
        if atom.name in ['C4', 'C6', 'C16', 'C22']:
            adjusted_charges[i] += charge_adjustment
            charge_adjustments[atom.name] = f"{raw_glutaraldehyde_charges[i]:.6f} → {adjusted_charges[i]:.6f}"
    
    # Print charge table
    print("GLUTARALDEHYDE CHARGES (CORRECTED)")
    print("="*50)
    print(f"Number of atoms: {len(glutaraldehyde_atoms)}")
    print("="*50)
    print(f"{'Atom':<8} {'Type':<8} {'Charge':>12} {'Corrected':>12}")
    print("-"*60)
    
    for i, atom in enumerate(glutaraldehyde_atoms):
        atom_name = atom.name
        original_charge = raw_glutaraldehyde_charges[i]
        adjusted_charge = adjusted_charges[i]
        
        if atom_name in charge_adjustments:
            print(f"{atom_name:<8} {atom.type:<8} {original_charge:12.6f} {adjusted_charge:12.6f}*")
        else:
            print(f"{atom_name:<8} {atom.type:<8} {original_charge:12.6f} {adjusted_charge:12.6f}")
    
    print("-"*60)
    print(f"{'SUM':<8} {'':<8} {np.sum(raw_glutaraldehyde_charges):12.6f} {np.sum(adjusted_charges):12.6f}")
    print("* = corrected charge")
    print("="*50 + "\n")
    
    # Write to file
    with open('glutaraldehyde_charges.txt', 'w') as f:
        f.write("GLUTARALDEHYDE CHARGES (CORRECTED)\n")
        f.write("="*50 + "\n")
        f.write(f"Number of atoms: {len(glutaraldehyde_atoms)}\n")
        f.write("="*50 + "\n")
        f.write(f"{'Atom':<8} {'Type':<8} {'Charge':>12} {'Corrected':>12}\n")
        f.write("-"*60 + "\n")
        
        for i, atom in enumerate(glutaraldehyde_atoms):
            atom_name = atom.name
            original_charge = raw_glutaraldehyde_charges[i]
            adjusted_charge = adjusted_charges[i]
            
            if atom_name in charge_adjustments:
                f.write(f"{atom_name:<8} {atom.type:<8} {original_charge:12.6f} {adjusted_charge:12.6f}*\n")
            else:
                f.write(f"{atom_name:<8} {atom.type:<8} {original_charge:12.6f} {adjusted_charge:12.6f}\n")
        
        f.write("-"*60 + "\n")
        f.write(f"{'SUM':<8} {'':<8} {np.sum(raw_glutaraldehyde_charges):12.6f} {np.sum(adjusted_charges):12.6f}\n")
        f.write("* = corrected charge\n")
        f.write("="*50 + "\n")
    
    print("glutaraldehyde_charges.txt file generated")
    return adjusted_charges

def step7_build_pva_chain() -> str:
    """Step 7: Build PVA chain of length 7 using pva_builder.py"""
    print("\n=== STEP 7: Building PVA chain of length 7 ===")
    
    # Build PVA chain
    atoms, bonds = build_pva(n=7, print_info=True)
    
    pdb_file = 'PVA7.pdb'
    print(f"Built PVA chain with {len(atoms)} atoms")
    return pdb_file

def step8_minimize_pva_structure(pdb_file: str) -> str:
    """Step 8: Minimize PVA structure using obabel"""
    print("\n=== STEP 8: Minimizing PVA structure ===")
    
    base_name = os.path.splitext(pdb_file)[0]
    minimized_pdb = f'{base_name}_min.pdb'
    
    cmd = f'obabel {pdb_file} -opdb -O {minimized_pdb} --minimize --ff mmff94'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Obabel minimization failed: {result.stderr}")
    
    print(f"Minimized PVA structure saved as {minimized_pdb}")
    return minimized_pdb

def step8_run_antechamber_pva(pdb_file: str) -> str:
    """Step 8: Run antechamber to produce PVA7_min.mol2 file"""
    print("\n=== STEP 8: Running antechamber for PVA ===")
    
    mol2_file = 'PVA7_min.mol2'
    
    # Run antechamber
    cmd = f"antechamber -j 4 -at gaff2 -dr no -fi pdb -fo mol2 -i {pdb_file} -o {mol2_file} -c bcc"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Antechamber failed: {result.stderr}")
    
    print(f"Generated {mol2_file} with antechamber charges")
    
    # Check charges
    mol2 = mda.Universe(mol2_file)
    print(f"Sum of charges in {mol2_file}: {np.sum(mol2.atoms.charges)}")
    
    return mol2_file

def step9_extract_pva_monomer_charges(mol2_file: str):
    """Step 9: Extract PVA monomer charges for specified atoms in order"""
    print("\n=== STEP 9: Extracting PVA monomer charges ===")
    
    # Load the mol2 file
    mol2 = mda.Universe(mol2_file)
    mol2.guess_TopologyAttrs(context='default', to_guess=['elements'])
    
    # Define PVA monomer atoms from PVA chain in CH2-CHOH-CH2 order
    # These atoms come from the PVA7_min.mol2 (minimized PVA chain structure)
    # CH2 GROUP 1: C10, H16, H17
    # CHOH GROUP: C11, H18, O4, H39  
    # CH2 GROUP 2: C12, H19, H20

    ch2_1_atoms = ['C10', 'H19', 'H20']
    choh_atoms = ['C11', 'H21', 'O4', 'H22']
    ch2_2_atoms = ['C12', 'H23', 'H24']

    pva_atoms = ch2_1_atoms + choh_atoms + ch2_2_atoms
    _validate_atom_count("PVA monomer", len(pva_atoms), PVA_ATOMS_PER_MONOMER)
    
    pva_charges = []
    original_charges = []
    
    # Select PVA atoms
    pva_selection = " or ".join([f"name {name}" for name in pva_atoms])
    pva_mol_atoms = mol2.select_atoms(pva_selection)
    
    # Check if we found all atoms
    if len(pva_mol_atoms) != len(pva_atoms):
        raise RuntimeError(
            f"PVA monomer extraction mismatch: expected {len(pva_atoms)} named atoms, "
            f"found {len(pva_mol_atoms)}. Available atoms: {list(pva_mol_atoms.names)}"
        )
    
    # Get charges for all atoms
    
    for atom_name in pva_atoms:
        matching_atoms = pva_mol_atoms[pva_mol_atoms.names == atom_name]
        if len(matching_atoms) > 0:
            atom = matching_atoms[0]
            charge = atom.charge
            pva_charges.append(charge)
            original_charges.append(charge)
        else:
            pva_charges.append(0.0)
            original_charges.append(0.0)
    
    # Calculate charge correction to make net charge zero
    net_charge = np.sum(pva_charges)
    if abs(net_charge) > 1e-6:  # Only adjust if net charge is significant
        correction = -net_charge / 2  # Divide by 2 for C10 and C12
        # Apply correction to C10 and C12 (indices 0 and 7 in pva_charges)
        pva_charges[0] += correction  # C10
        pva_charges[7] += correction  # C12
    
    # Display charges
    if abs(net_charge) > 1e-6:
        print("PVA_MONOMER CHARGES (CORRECTED)")
    else:
        print("PVA_MONOMER CHARGES")
    print("="*50)
    print(f"{'Atom':<8} {'Charge':>12} {'Corrected':>12}")
    print("-"*50)
    
    print("\nCH2 GROUP 1:")
    for i, atom_name in enumerate(ch2_1_atoms):
        original = original_charges[i]
        corrected = pva_charges[i]
        star = "*" if corrected != original else " "
        print(f"{atom_name:<8} {original:12.6f} {corrected:12.6f}{star}")
    
    print("\nCHOH GROUP:")
    for i, atom_name in enumerate(choh_atoms):
        original = original_charges[i+3]
        corrected = pva_charges[i+3]
        star = "*" if corrected != original else " "
        print(f"{atom_name:<8} {original:12.6f} {corrected:12.6f}{star}")
    
    print("\nCH2 GROUP 2:")
    for i, atom_name in enumerate(ch2_2_atoms):
        original = original_charges[i+7]
        corrected = pva_charges[i+7]
        star = "*" if corrected != original else " "
        print(f"{atom_name:<8} {original:12.6f} {corrected:12.6f}{star}")
    
    print("-"*50)
    print(f"{'SUM':<8} {np.sum(original_charges):12.6f} {np.sum(pva_charges):12.6f}")
    if abs(net_charge) > 1e-6:
        print("* = corrected charge")
    print("="*50 + "\n")
    
    # Write to file with grouped structure format
    with open('PVA_monomercharges.txt', 'w') as f:
        if abs(net_charge) > 1e-6:
            f.write("PVA_MONOMER CHARGES (CORRECTED)\n")
        else:
            f.write("PVA_MONOMER CHARGES\n")
        f.write("="*50 + "\n")
        f.write(f"{'Atom':<8} {'Type':<8} {'Charge':>12} {'Corrected':>12}\n")
        f.write("-"*60 + "\n")
        
        f.write("\nCH2 GROUP 1:\n")
        for i, atom_name in enumerate(ch2_1_atoms):
            original = original_charges[i]
            corrected = pva_charges[i]
            star = "*" if corrected != original else " "
            atom = pva_mol_atoms[pva_mol_atoms.names == atom_name][0]
            f.write(f"{atom_name:<8} {atom.type:<8} {original:12.6f} {corrected:12.6f}{star}\n")
        
        f.write("\nCHOH GROUP:\n")
        for i, atom_name in enumerate(choh_atoms):
            original = original_charges[i+3]
            corrected = pva_charges[i+3]
            star = "*" if corrected != original else " "
            atom = pva_mol_atoms[pva_mol_atoms.names == atom_name][0]
            f.write(f"{atom_name:<8} {atom.type:<8} {original:12.6f} {corrected:12.6f}{star}\n")
        
        f.write("\nCH2 GROUP 2:\n")
        for i, atom_name in enumerate(ch2_2_atoms):
            original = original_charges[i+7]
            corrected = pva_charges[i+7]
            star = "*" if corrected != original else " "
            atom = pva_mol_atoms[pva_mol_atoms.names == atom_name][0]
            f.write(f"{atom_name:<8} {atom.type:<8} {original:12.6f} {corrected:12.6f}{star}\n")
        
        f.write("-"*60 + "\n")
        f.write(f"{'SUM':<8} {'':<8} {np.sum(original_charges):12.6f} {np.sum(pva_charges):12.6f}\n")
        if abs(net_charge) > 0.0001:
            f.write("* = corrected charge\n")
        f.write("="*50 + "\n")
    
    print("PVA_monomercharges.txt file generated")
    return pva_charges

def main():
    """Main function to run all steps"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cleanup_beginning(script_dir)

    smiles = "CCC(O)CC1CC(CC(O)CC)OC(CCCC2OC(CC(O)CC)CC(CC(O)CC)O2)O1"
    
    try:
        # Step 1-2: Convert SMILES to PDB
        pdb_file = step1_convert_smiles_to_pdb(smiles)
        
        # Step 3: Minimize structure
        minimized_pdb = step3_minimize_structure(pdb_file)
        
        # Step 4: Run antechamber
        mol2_file = step4_run_antechamber(minimized_pdb)
        
        # Step 5: Extract PVA terminal group charges
        step5_extract_pva_terminal_group_charges(mol2_file)
        
        # Step 6: Extract glutaraldehyde charges
        charges = step6_extract_glutaraldehyde_charges(mol2_file)
        
        print(f"\n=== CROSSLINKED STRUCTURE PROCESSING COMPLETE ===")
        print(f"Processed {len(charges)} glutaraldehyde atoms")
        print("Crosslinked structure charge files generated successfully!")
        
        # Step 7: Build PVA chain of length 7
        pva_pdb_file = step7_build_pva_chain()
        
        # Step 8: Minimize PVA structure and run antechamber
        pva_minimized_pdb = step8_minimize_pva_structure(pva_pdb_file)
        pva_mol2_file = step8_run_antechamber_pva(pva_minimized_pdb)
        
        # Step 9: Extract PVA monomer charges
        pva_charges = step9_extract_pva_monomer_charges(pva_mol2_file)
        
        print(f"\n=== PVA MONOMER PROCESSING COMPLETE ===")
        print(f"Processed {len(pva_charges)} PVA monomer atoms")
        print("PVA monomer charge file generated successfully!")
        
        print(f"\n=== ALL PROCESSING COMPLETE ===")
        print("Generated charge files:")
        print("- PVA_terminal_group_charges.txt")
        print("- glutaraldehyde_charges.txt") 
        print("- PVA_monomercharges.txt")
        
        # Clean up temporary files at the end (use script dir; cwd may differ if run via wrapper)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        cleanup_end(script_dir)
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
