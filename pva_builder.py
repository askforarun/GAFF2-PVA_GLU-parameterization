#!/usr/bin/env python3
"""
PVA Builder Module

This module provides tools for building Polyvinyl Alcohol (PVA) polymer structures
with hard-coded geometry. The PVA polymer consists of CH3 terminal groups and 
CH2-CHOH-CH2 repeat units, representing polyvinyl alcohol with alternating hydroxyl-bearing carbons and proper methyl termination.

PVA Polymer Structure:
    Full Polymer: CH3-(CH2-CHOH-CH2)n-CH3
    Repeat Unit: CH2-CHOH-CH2
    Represents: Polyvinyl alcohol with alternating hydroxyl-bearing carbons and CH3 termination
    
Atom Ordering:
    - Repeat units follow CH2-CHOH-CH2 pattern in the middle section
    - Within each repeat, atoms are written as CHHCHOHCHH (C, H, H, C, H, O, H, C, H, H)
    - Last 8 atoms contain CH3 terminal groups in CHHHCHHH order
    
Usage Examples:
    # Build PVA with 3 repeat units (11 carbon backbone)
    atoms, bonds = build_pva(n=3)
    
    # Build PVA with custom output file
    atoms, bonds = build_pva(n=5, output_file="custom_pva.pdb")
    
    # Command line usage
    python pva_builder.py --n 7 --out my_pva.pdb
    
Output Files:
    - PVA{n}.pdb: Generated PVA polymer structure with CH3 termination (cap=True)
    - PVA{n}_trim.pdb: Generated uncapped structure when cap=False
    - Contains CONECT records for bond information
    
Functions:
    - build_pva(): Build PVA polymer chain with CH3-(CH2-CHOH-CH2)n-CH3 structure
    - write_pdb_with_conect(): Write PDB file with bond information
    
Note:
    - Geometry is hard-coded, no SMILES processing or external dependencies required
    - Uses standard bond lengths and tetrahedral geometry
    - Suitable as a coarse scaffold for molecular dynamics simulations
"""

import numpy as np
from collections import Counter
from typing import List, Dict, Tuple
import argparse






def build_pva(
    n: int,
    print_info: bool = False,
    output_file: str = None,
    cap: bool = True,
) -> Tuple[List[dict], Dict[int, List[int]]]:
    """
    Build PVA polymer with CH3-(CH2-CHOH-CH2)n-CH3 structure.
    
    This function creates a Polyvinyl Alcohol polymer chain with CH3 terminal groups
    and CH2-CHOH-CH2 repeat units in the middle. The structure follows the pattern:
        CH3-(CH2-CHOH-CH2)n-CH3
    
    The polymer structure:
        - CH3-(CH2-CHOH-CH2)n-CH3
        - Total carbon atoms: 3n + 2 (n repeat units + 2 CH3 terminals)
        - Terminal groups: CH3 at both ends
        - Middle units: CH2-CHOH-CH2 pattern
        - Atom ordering: Repeat units in per-unit CHHCHOHCHH order, last 8 atoms contain CH3 groups in CHHHCHHH order
    
    Args:
        n (int): Number of CH2-CHOH-CH2 repeat units in the polymer
                Must be positive integer. Typical values: 3-20
        print_info (bool): Whether to print group information and statistics
                          Default: False
        output_file (str): Custom output file path. If None, uses PVA{n}.pdb
                         (or PVA{n}_trim.pdb when cap=False). Default: None
        cap (bool): Whether to keep terminal CH3 groups. When False, the final
                    8 atoms (CH3-CH3) are removed to produce uncapped ends.
                    Default: True
    
    Returns:
        Tuple[List[dict], Dict[int, List[int]]]: 
            - List of atom dictionaries with keys: name, element, x, y, z
            - Dictionary of bond connections (1-based atom indices)
    
    Raises:
        ValueError: If n is not a positive integer
    
    Example:
        >>> atoms, bonds = build_pva(n=3)
        >>> print(f"Built PVA with {len(atoms)} atoms")
        CH3 groups: 2 at positions [1, 11]
        CH2 groups: 6 at positions [2, 4, 5, 7, 8, 10]
        CHOH groups: 3 at positions [3, 6, 9]
        Total atoms: 38
        
    Note:
        - Uses standard bond lengths: C-C=1.54Å, C-H=1.09Å, C-O=1.43Å, O-H=0.96Å
        - Atoms arranged in straight line along x-axis
        - Hydroxyl groups placed in tetrahedral geometry
        - Automatically writes PVA{n}.pdb file with CONECT records (or PVA{n}_trim.pdb if cap=False)
        - Terminal groups are CH3 (not CH2) for proper polymer termination
        - Last 8 atoms contain CH3 terminal groups in CHHHCHHH order when cap=True
        - OH hydrogens placed before CH3 groups to maintain atom ordering requirements
    """
    
    # Validate input
    if not isinstance(n, int) or n <= 0:
        raise ValueError(f"n must be a positive integer, got {n}")
    
    # Note: Geometry is hard-coded, no SMILES processing needed

    # Standard bond lengths
    cc = 1.54  # C-C bond length
    ch = 1.09  # C-H bond length  
    co = 1.43  # C-O bond length
    oh = 0.96  # O-H bond length

    # Generate backbone types for CH3-(CH2-CHOH-CH2)n-CH3 structure
    # Pattern: CH3 + (CH2-CHOH-CH2)*n + CH3
    # Build repeat units first, then terminal CH3 groups appended last
    types = ["CH3"] + ["CH2", "CHOH", "CH2"] * n + ["CH3"]
    nC = len(types)

    # Build carbon backbone in a straight line
    bpos = np.zeros((nC, 3))
    for i in range(nC):
        bpos[i] = np.array([i * cc, 0.0, 0.0])

    # Build atoms in specific order:
    # 1. All repeat units (CH2-CHOH-CH2 pattern) in CHHCHOHCHH per-repeat order
    # 2. Both CH3 groups at the end (ensures last 8 atoms are CHHHCHHH)
    atoms = []
    bonds = {}
    h_counter = 1
    o_counter = 1

    def add_atom(name: str, element: str, pos: np.ndarray) -> int:
        atoms.append({
            "name": name,
            "element": element, 
            "x": pos[0],
            "y": pos[1],
            "z": pos[2]
        })
        return len(atoms)

    def add_bond(i: int, j: int):
        if i not in bonds:
            bonds[i] = []
        if j not in bonds:
            bonds[j] = []
        bonds[i].append(j)
        bonds[j].append(i)

    c_ids = []
    previous_cid = None
    
    # Build repeat units (CH2-CHOH-CH2 pattern) first
    for i in range(1, nC-1):  # Skip first and last CH3
        t = types[i]
        cpos = bpos[i]
        
        # Add carbon atom
        cid = add_atom(f"C{i}", "C", cpos)  # Renumber to skip CH3 positions
        c_ids.append(cid)
        
        # Bond to previous carbon
        if previous_cid is not None:
            add_bond(previous_cid, cid)
        previous_cid = cid

        # Alternate substituent orientation along the backbone to reduce steric clashes
        orient = 1.0 if (i % 2 == 0) else -1.0

        if t == "CH2":
            # CH2 GROUP: Carbon, Hydrogen, Hydrogen
            h1_pos = cpos + ch * np.array([0.0, 0.945 * orient, 0.326])
            h2_pos = cpos + ch * np.array([0.0, -0.945 * orient, -0.326])
            
            add_atom(f"H{h_counter}", "H", h1_pos)
            add_bond(cid, len(atoms))
            h_counter += 1
            
            add_atom(f"H{h_counter}", "H", h2_pos)
            add_bond(cid, len(atoms))
            h_counter += 1

        elif t == "CHOH":
            # CHOH GROUP: Carbon, Hydrogen, Oxygen, Hydrogen
            
            # Hydrogen on carbon
            h_c_pos = cpos + ch * np.array([0.0, 0.945 * orient, 0.326])
            add_atom(f"H{h_counter}", "H", h_c_pos)
            add_bond(cid, len(atoms))
            h_counter += 1

            # Oxygen
            o_pos = cpos + co * np.array([0.0, -0.945 * orient, -0.326])
            oid = add_atom(f"O{o_counter}", "O", o_pos)
            add_bond(cid, oid)
            o_counter += 1

            # OH hydrogen added immediately to preserve CHHCHOHCHH ordering
            h_oh_pos = o_pos + oh * np.array([0.0, 0.326 * orient, 0.945])
            add_atom(f"H{h_counter}", "H", h_oh_pos)
            add_bond(oid, len(atoms))
            h_counter += 1

    # Now add the first CH3 group (will be part of last 8 atoms)
    first_ch3_pos = bpos[0]
    first_cid = add_atom(f"C{nC-1}", "C", first_ch3_pos)  # Use high number to place at end
    c_ids.append(first_cid)
    
    # Bond first CH3 to the first repeat unit carbon
    if len(c_ids) > 1:
        add_bond(first_cid, c_ids[0])
    
    # Add first CH3 hydrogens
    h1_pos = first_ch3_pos + ch * np.array([0.0, 0.945, 0.326])
    h2_pos = first_ch3_pos + ch * np.array([0.0, -0.473, -0.816])
    h3_pos = first_ch3_pos + ch * np.array([0.0, -0.473, 0.490])
    
    add_atom(f"H{h_counter}", "H", h1_pos)
    add_bond(first_cid, len(atoms))
    h_counter += 1
    
    add_atom(f"H{h_counter}", "H", h2_pos)
    add_bond(first_cid, len(atoms))
    h_counter += 1
    
    add_atom(f"H{h_counter}", "H", h3_pos)
    add_bond(first_cid, len(atoms))
    h_counter += 1

    # Now add the last CH3 group (the very last atoms)
    last_ch3_pos = bpos[-1]
    last_cid = add_atom(f"C{nC}", "C", last_ch3_pos)
    c_ids.append(last_cid)
    
    # Bond last CH3 to previous carbon
    add_bond(previous_cid, last_cid)
    
    # Add last CH3 hydrogens (these will be the very last 3 atoms)
    h1_pos = last_ch3_pos + ch * np.array([0.0, 0.945, 0.326])
    h2_pos = last_ch3_pos + ch * np.array([0.0, -0.473, -0.816])
    h3_pos = last_ch3_pos + ch * np.array([0.0, -0.473, 0.490])
    
    add_atom(f"H{h_counter}", "H", h1_pos)
    add_bond(last_cid, len(atoms))
    h_counter += 1
    
    add_atom(f"H{h_counter}", "H", h2_pos)
    add_bond(last_cid, len(atoms))
    h_counter += 1
    
    add_atom(f"H{h_counter}", "H", h3_pos)
    add_bond(last_cid, len(atoms))
    h_counter += 1

    if not cap:
        # Remove the terminal CH3 groups (last 8 atoms) and associated bonds.
        if len(atoms) < 8:
            raise ValueError(f"Cannot trim terminal CH3 groups: only {len(atoms)} atoms present.")
        keep_n = len(atoms) - 8
        atoms = atoms[:keep_n]
        trimmed_bonds: Dict[int, List[int]] = {}
        for atom_id, neighbors in bonds.items():
            if atom_id > keep_n:
                continue
            kept_neighbors = [nbr for nbr in neighbors if nbr <= keep_n]
            if kept_neighbors:
                trimmed_bonds[atom_id] = kept_neighbors
        bonds = trimmed_bonds

    bonds_1based = {i: sorted(v) for i, v in bonds.items()}

    if print_info:
        # Count groups based on the backbone pattern
        ch3_count = 0
        ch2_count = 0
        choh_count = 0
        ch3_positions = []
        ch2_positions = []
        choh_positions = []
        
        # Use the backbone types to determine group positions
        for i, t in enumerate(types):
            carbon_pos = i + 1  # 1-based carbon position
            if t == "CH3":
                ch3_count += 1
                ch3_positions.append(carbon_pos)
            elif t == "CH2":
                ch2_count += 1
                ch2_positions.append(carbon_pos)
            elif t == "CHOH":
                choh_count += 1
                choh_positions.append(carbon_pos)

        print(f"CH3 groups: {ch3_count} at positions {ch3_positions}")
        print(f"CH2 groups: {ch2_count} at positions {ch2_positions}")
        print(f"CHOH groups: {choh_count} at positions {choh_positions}")
        print(f"Total atoms: {len(atoms)}")
        print(f"Structure follows CH3-(CH2-CHOH-CH2)n-CH3 pattern with {n} repeat units")
        print(f"Terminal groups: {'CH3 at both ends' if cap else 'trimmed (no CH3 caps)'}")
        if cap:
            print("Last 8 atoms contain CH3 terminal groups in order CHHHCHHH")

    # Auto-write PDB file with custom or default filename
    if output_file is None:
        output_file = f"PVA{n}_trim.pdb" if not cap else f"PVA{n}.pdb"
    write_pdb_with_conect(atoms, bonds_1based, output_file)

    return atoms, bonds_1based


def write_pdb_with_conect(atoms: List[dict], bonds: Dict[int, List[int]], out_pdb: str):
    """
    Write atoms and bonds to PDB file with CONECT records.
    
    This function creates a standard PDB file containing atomic coordinates
    and CONECT records that define the bond connectivity between atoms.
    The CONECT records ensure that bond information is preserved when
    the PDB file is loaded into other molecular visualization or analysis tools.
    
    Args:
        atoms (List[dict]): List of atom dictionaries with keys:
                           - name: Atom name (e.g., "C1", "H1", "O1")
                           - element: Element symbol (e.g., "C", "H", "O")
                           - x, y, z: Cartesian coordinates in Angstroms
        bonds (Dict[int, List[int]]): Dictionary mapping atom IDs (1-based) 
                                     to lists of bonded atom IDs
        out_pdb (str): Output PDB filename
    
    Output Format:
        - HETATM records: Standard PDB atom coordinate format
        - CONECT records: Bond connectivity information
        - END record: File termination
        
    Example:
        >>> atoms = [{'name': 'C1', 'element': 'C', 'x': 0.0, 'y': 0.0, 'z': 0.0}]
        >>> bonds = {1: []}
        >>> write_pdb_with_conect(atoms, bonds, "molecule.pdb")
        # Creates molecule.pdb with atom and bond information
        
    Note:
        - PDB title reflects CH3-(CH2-CHOH-CH2)n-CH3 structure
        - Maintains proper atom ordering for CH3 terminal groups
    """
    with open(out_pdb, "w") as f:
        f.write("TITLE     PVA Polymer - CH3-(CH2-CHOH-CH2)n-CH3 Pattern\n")
        for i, a in enumerate(atoms, start=1):
            f.write(
                f"HETATM{i:5d} {a['name']:<4} UNK A   1    "
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                f"  1.00  0.00          {a['element']:>2}\n"
            )
        for i in range(1, len(atoms) + 1):
            neigh = bonds.get(i, [])
            if not neigh:
                f.write(f"CONECT{i:5d}\n")
            else:
                for k in range(0, len(neigh), 4):
                    f.write(
                        f"CONECT{i:5d}"
                        + "".join(f"{j:5d}" for j in neigh[k:k + 4])
                        + "\n"
                    )
        f.write("END\n")


def trim_terminal_ch3(pdb_in: str, pdb_out: str = None) -> str:
    """
    Remove the last 8 atoms (two terminal CH3 groups) from a PVA PDB while
    preserving the ordering of all other atoms and rewriting CONECT records.

    Args:
        pdb_in (str): Path to input PVA PDB (any chain length built by this module)
        pdb_out (str): Output PDB path. If None, uses \"<input>_trim.pdb\".

    Returns:
        str: Path to the written trimmed PDB file.

    Raises:
        ValueError: If the file has fewer than 8 atoms.
    """
    with open(pdb_in, "r") as f:
        lines = f.readlines()

    prefix = []
    atom_lines = []
    conect_lines = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            atom_lines.append(line)
        elif line.startswith("CONECT"):
            conect_lines.append(line)
        elif not line.startswith("END"):
            prefix.append(line)

    if len(atom_lines) < 8:
        raise ValueError(f"{pdb_in} has only {len(atom_lines)} atoms; need at least 8 to trim.")

    kept_atoms = atom_lines[:-8]
    # Map old serial -> new serial
    serial_map = {}
    rewritten_atoms = []
    for new_idx, line in enumerate(kept_atoms, start=1):
        old_serial = int(line[6:11])
        serial_map[old_serial] = new_idx
        rewritten_atoms.append(f"{line[:6]}{new_idx:5d}{line[11:]}")

    # Rewrite CONECT lines with updated serials, dropping references to trimmed atoms
    rewritten_conect = []
    for line in conect_lines:
        tokens = line.split()
        # tokens[0] == "CONECT"
        mapped = []
        for tok in tokens[1:]:
            try:
                sid = int(tok)
            except ValueError:
                continue
            if sid in serial_map:
                mapped.append(serial_map[sid])
        center_old = int(tokens[1])
        if center_old not in serial_map:
            continue
        center_new = serial_map[center_old]
        if not mapped:
            rewritten_conect.append(f"CONECT{center_new:5d}\n")
        else:
            # max 4 neighbors per line; split if needed
            for i in range(0, len(mapped), 4):
                chunk = mapped[i:i + 4]
                rewritten_conect.append(
                    f"CONECT{center_new:5d}" + "".join(f"{nbr:5d}" for nbr in chunk) + "\n"
                )

    if pdb_out is None:
        if pdb_in.lower().endswith(".pdb"):
            pdb_out = pdb_in[:-4] + "_trim.pdb"
        else:
            pdb_out = pdb_in + "_trim.pdb"

    with open(pdb_out, "w") as f:
        for line in prefix:
            f.write(line)
        for line in rewritten_atoms:
            f.write(line)
        for line in rewritten_conect:
            f.write(line)
        f.write("END\n")

    return pdb_out


def main():
    """
    Command line interface for PVA builder.
    
    This function provides a command-line interface for building PVA polymers
    with specified repeat units and hard-coded geometry.
    
    Usage:
        python pva_builder.py --n 7 --out my_pva.pdb
    
    Arguments:
        --n: Number of CH2-CHOH-CH2 repeat units (required)
        --out: Output PDB filename (optional, defaults to PVA{n}.pdb or PVA{n}_trim.pdb if --no-cap)
        --cap/--no-cap: Include or trim terminal CH3 groups (optional)
        --verbose: Print detailed polymer information (optional flag)
    
    Examples:
        # Build PVA with 7 repeat units
        python pva_builder.py --n 7
        
        # Build PVA with custom output filename
        python pva_builder.py --n 5 --out custom_pdb.pdb
        
        # Build PVA with verbose output
        python pva_builder.py --n 3 --verbose
        
    Output:
        - Generates PVA{n}.pdb file with CH3-(CH2-CHOH-CH2)n-CH3 structure
          (or PVA{n}_trim.pdb when --no-cap is used)
        - Prints group statistics and atom counts
        - Creates CONECT records for bond information
        - Ensures proper atom ordering with CH3 terminal groups in last 8 atoms
        
    Note:
        - The resulting polymer has CH3 terminal groups (not CH2)
        - Atom ordering follows CHHCHOHCHH per repeat unit
        - Last 8 atoms contain CH3 groups in CHHHCHHH order
        - Geometry is hard-coded, no external dependencies required
    """
    ap = argparse.ArgumentParser(description="Build PVA with CH3-(CH2-CHOH-CH2)n-CH3 pattern")
    ap.add_argument("--n", type=int, required=True, help="Number of CH2-CHOH-CH2 repeat units")
    ap.add_argument("--out", default=None, help="Output PDB file")
    ap.add_argument("--verbose", action="store_true", help="Print detailed information about the polymer")
    ap.add_argument(
        "-cap",
        "--cap",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include terminal CH3 caps (default: True). Use --no-cap to trim them.",
    )
    args = ap.parse_args()

    out = args.out or (f"PVA{args.n}.pdb" if args.cap else f"PVA{args.n}_trim.pdb")
    atoms, bonds = build_pva(
        n=args.n,
        print_info=args.verbose,
        output_file=out,
        cap=args.cap,
    )
    
    counts = Counter(a["element"] for a in atoms)
    print(f"Wrote {out}")
    print("Counts:", dict(counts), "Total:", len(atoms))


if __name__ == "__main__":
    main()
