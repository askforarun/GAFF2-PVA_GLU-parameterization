import subprocess
import os
import re
from collections import Counter
from pathlib import Path

import MDAnalysis as mda
import numpy as np

try:
    from src.system_constants import GLU_ATOMS_PER_MOLECULE, PVA_ATOMS_PER_MONOMER
except ImportError:
    from system_constants import GLU_ATOMS_PER_MOLECULE, PVA_ATOMS_PER_MONOMER

REFERENCE_DATA = Path(__file__).parent / "charge_data"


# ---------------------------------------------------------------------------
# External tool wrappers
# ---------------------------------------------------------------------------

def generate_monomer_from_smiles(smiles: str, output_name: str = "structure") -> str:
    """Generate a PDB file from a SMILES string using obabel."""
    output_pdb = f"{output_name}.pdb"
    result = subprocess.run(
        f"obabel -:'{smiles}' -opdb --gen3d -O {output_pdb}",
        shell=True, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"obabel failed: {result.stderr}")
    if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
        raise RuntimeError(f"obabel failed to generate PDB file: {output_pdb}")
    return output_pdb


def run_tleap_with_error_check(input_file: str = "tleap.in"):
    """Run tleap and raise RuntimeError if leap.log contains errors."""
    subprocess.run(f"tleap -f {input_file}", shell=True, capture_output=True, text=True)
    if os.path.exists("leap.log"):
        with open("leap.log") as f:
            log = f.read()
        match = re.search(r"Errors = (\d+)", log)
        if match and int(match.group(1)) > 0:
            raise RuntimeError(f"tleap failed with {match.group(1)} errors. Check leap.log.")
        if "FATAL" in log.upper():
            raise RuntimeError("tleap encountered a fatal error. Check leap.log.")


def getfiles(pdb_file: str):
    """Run antechamber + parmchk2 + tleap to produce mol2/frcmod/top/crd for a PDB."""
    base = pdb_file[:-4] if pdb_file.endswith(".pdb") else pdb_file
    ff = "gaff2"
    if os.path.exists("leap.log"):
        os.remove("leap.log")
    subprocess.run(
        f"antechamber -j 4 -at {ff} -dr no -fi pdb -fo mol2 -i {base}.pdb -o {base}.mol2",
        shell=True,
    )
    subprocess.run(f"parmchk2 -i {base}.mol2 -o {base}.frcmod -f mol2 -a Y", shell=True)
    with open("tleap.in", "w") as f:
        f.write(f"source leaprc.{ff}\n")
        f.write(f"SUS = loadmol2 {base}.mol2\n")  # SUS: arbitrary tleap handle name
        f.write("check SUS\n")
        f.write(f"loadamberparams {base}.frcmod\n")
        f.write(f"saveamberparm SUS {base}.top {base}.crd\n")
        f.write("quit")
    run_tleap_with_error_check()


# ---------------------------------------------------------------------------
# Charge loading
# ---------------------------------------------------------------------------

def _read_corrected_charges(filepath: str) -> list:
    """
    Parse corrected partial charges from a reference charge file.

    Supports two column layouts:
      - New format: Atom  Type  OrigCharge  CorrectedCharge[*]
      - Old format: Atom  OrigCharge  CorrectedCharge[*]

    Header/separator/summary lines are silently skipped.
    """
    charges = []
    with open(filepath) as f:
        for line in f:
            if (line.startswith(("#", "=", "-", "*"))
                    or not line.strip()
                    or any(kw in line for kw in ("GROUP", "SUM", "CORRECTED"))):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    float(parts[2])
                    charges.append(float(parts[3].rstrip("*")))
                except ValueError:
                    continue
            elif len(parts) >= 3:
                try:
                    float(parts[1])
                    charges.append(float(parts[2].rstrip("*")))
                except ValueError:
                    continue
    return charges


def load_system_charges(
    chain_length: int,
    n_pva: int,
    n_glu: int,
    pva_charge_file: str = None,
    glu_charge_file: str = None,
) -> np.ndarray:
    """
    Load and validate partial charges for a PVA-glutaraldehyde hydrogel system.

    Returns a 1-D numpy array:
        [ PVA_chain_1 | ... | PVA_chain_n_pva | GLU_1 | ... | GLU_n_glu ]

    Raises RuntimeError if a charge file is missing or empty.
    Raises ValueError if charge count or neutrality checks fail.
    """
    if pva_charge_file is None:
        pva_charge_file = str(REFERENCE_DATA / "PVA_monomercharges.txt")
    if glu_charge_file is None:
        glu_charge_file = str(REFERENCE_DATA / "glutaraldehyde_charges.txt")

    if not os.path.exists(pva_charge_file):
        raise RuntimeError(f"PVA charge file not found: {pva_charge_file}")
    monomer_charges = _read_corrected_charges(pva_charge_file)
    if len(monomer_charges) != PVA_ATOMS_PER_MONOMER:
        raise ValueError(
            f"Expected {PVA_ATOMS_PER_MONOMER} monomer charges, found {len(monomer_charges)}"
        )
    single_chain_charges = monomer_charges * chain_length
    total_pva = np.sum(single_chain_charges)
    print(f"Total PVA charge (1 chain, {chain_length} monomers): {total_pva:.10f}")
    if abs(total_pva) > 1e-8:
        raise ValueError(f"PVA chain is not neutral: total charge = {total_pva:.10f}")

    if not os.path.exists(glu_charge_file):
        raise RuntimeError(f"GLU charge file not found: {glu_charge_file}")
    single_glu_charges = _read_corrected_charges(glu_charge_file)
    if not single_glu_charges:
        raise RuntimeError(f"No valid GLU charges found in {glu_charge_file}")
    if len(single_glu_charges) != GLU_ATOMS_PER_MOLECULE:
        raise ValueError(
            f"Expected {GLU_ATOMS_PER_MOLECULE} GLU charges, found {len(single_glu_charges)}"
        )
    total_glu = np.sum(single_glu_charges)
    print(f"Total GLU charge (1 molecule): {total_glu:.10f}")
    if abs(total_glu) > 1e-8:
        raise ValueError(f"GLU molecule is not neutral: total charge = {total_glu:.10f}")

    all_charges = np.array(
        single_chain_charges * n_pva + single_glu_charges * n_glu, dtype=float
    )
    total = np.sum(all_charges)
    print(f"Total system charge ({n_pva} PVA + {n_glu} GLU): {total:.10f}")
    if abs(total) > 1e-6:
        print(f"WARNING: System is not neutral! Total charge = {total:.10f}")
    else:
        print("System is charge neutral")
    return all_charges


# ---------------------------------------------------------------------------
# Force-field parameter writing
# ---------------------------------------------------------------------------

def _parse_frcmod(path: str) -> dict:
    """Read an AMBER frcmod file **once** and return section data.

    Returns a dict with:
    - ``'BOND'``:     non-empty data lines from the BOND section.
    - ``'ANGLE'``:    non-empty data lines from the ANGLE section.
    - ``'DIHE_raw'``: all non-empty lines from the DIHE section onwards,
                      stopping just before NONBON.  The IMPR header line
                      ("IMPR") is intentionally *included* in this list so
                      that downstream code can detect it with an
                      ``'IMP' in line`` guard, exactly as the original
                      file-iteration logic did.
    """
    # IMPR is intentionally absent so its header line falls through into DIHE_raw.
    TOP_HEADERS = {'MASS', 'BOND', 'ANGLE', 'DIHE', 'NONBON'}
    result: dict = {'BOND': [], 'ANGLE': [], 'DIHE_raw': []}
    current = None
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            first = line.split()[0]
            if first in TOP_HEADERS:
                if first == 'NONBON':
                    break
                current = first
                continue  # skip the section-header line itself
            if current == 'BOND':
                result['BOND'].append(line)
            elif current == 'ANGLE':
                result['ANGLE'].append(line)
            elif current == 'DIHE':       # also captures IMPR header + IMPR data
                result['DIHE_raw'].append(line)
    return result


def print_ff_parameters_crosslinking(frcmod_file, parm_lammps_path):
    """
    Extract and append all force field parameters (bonds, angles, dihedrals)
    from an frcmod file to parm.lammps with proper section headers.

    The function reads the existing max bond/angle/dihedral coefficient IDs
    from parm_lammps_path to determine offsets, then appends the new
    crosslinking parameters. No external counts need to be passed in.

    Args:
        frcmod_file (str): Input frcmod file containing crosslinking parameters.
        parm_lammps_path (str): Path to parm.lammps (read for offsets, then appended to).

    Returns:
        tuple[int, int, int]: (n_new_bonds, n_new_angles, n_new_dihedrals)
            Number of new bond / angle / dihedral coefficient entries appended.
    """
    # ------------------------------------------------------------------ #
    # Read existing max coeff IDs from parm.lammps (determines offsets).  #
    # ------------------------------------------------------------------ #
    n_bonds = n_angles = n_dihedrals = 0
    with open(parm_lammps_path) as _f:
        for _line in _f:
            _parts = _line.split()
            if _parts and _parts[0] == 'bond_coeff':
                n_bonds = max(n_bonds, int(_parts[1]))
            elif _parts and _parts[0] == 'angle_coeff':
                n_angles = max(n_angles, int(_parts[1]))
            elif _parts and _parts[0] == 'dihedral_coeff':
                n_dihedrals = max(n_dihedrals, int(_parts[1]))

    # ------------------------------------------------------------------ #
    # Parse frcmod ONCE into section lists.                               #
    # ------------------------------------------------------------------ #
    parsed = _parse_frcmod(frcmod_file)
    bond_lines = parsed['BOND']
    angle_lines = parsed['ANGLE']
    dihe_raw = parsed['DIHE_raw']

    # dih_type_counts replicates the original ``dih_str`` list:
    # count every non-IMP line's first token (type string), including
    # IMPR-section data lines (they have no 'IMP' in the data, only the
    # 'IMPR' header does — and that line IS in dihe_raw but is filtered out).
    dih_type_counts: Counter = Counter(
        line.split()[0]
        for line in dihe_raw
        if line.split() and 'IMP' not in line
    )

    def _coeff(divisor: float, barrier: float):
        """Compute barrier/divisor; return 0 for zero divisor; None to skip."""
        if abs(divisor) > 1e-5:
            return barrier / divisor
        if int(divisor) == 0:
            return 0
        return None

    with open(parm_lammps_path, 'a') as out:

        # --- BONDS ---------------------------------------------------- #
        out.write('\n')
        out.write('# =============================================================================\n')
        out.write('# BOND PARAMETERS FOR CROSSLINKING\n')
        out.write('# Parameters for bonds created during crosslinking\n')
        out.write('# =============================================================================\n')
        out.write('\n')

        new_bond_id = 1
        for param_line in bond_lines:
            p = param_line.split()
            bond_id = n_bonds + new_bond_id
            out.write(f"bond_coeff {bond_id} {p[1]} {p[2]} # {p[0]}\n")
            new_bond_id += 1

        # --- ANGLES --------------------------------------------------- #
        out.write('\n')
        out.write('# =============================================================================\n')
        out.write('# ANGLE PARAMETERS FOR CROSSLINKING\n')
        out.write('# Parameters for angles created during crosslinking\n')
        out.write('# =============================================================================\n')
        out.write('\n')

        na = 1
        for angle_line in angle_lines:
            p = angle_line.split()
            out.write("angle_coeff {} {} {} # {}\n".format(
                n_angles + na, p[1], p[2], p[0]))
            na += 1

        # --- DIHEDRALS ------------------------------------------------ #
        out.write('\n')
        out.write('# =============================================================================\n')
        out.write('# DIHEDRAL PARAMETERS FOR CROSSLINKING\n')
        out.write('# Parameters for dihedrals created during crosslinking\n')
        out.write('# =============================================================================\n')
        out.write('\n')

        ndhid = 1
        # Use a shared iterator so the inner multi-term loop consumes
        # continuation lines from the same sequence — identical to the
        # original ``for line_2 in f`` advancing the open file handle.
        dihe_iter = iter(dihe_raw)
        for line_1 in dihe_iter:
            if 'IMP' in line_1:          # IMPR header stops the writing pass
                break
            p = line_1.split()
            if not p:
                continue
            dih_type = p[0]
            divisor   = float(p[1])
            barrier   = float(p[2])
            phase     = float(p[3])
            period    = float(p[4])

            if period > 0:
                coeff = _coeff(divisor, barrier)
                if coeff is not None:
                    out.write("dihedral_coeff {} 1 {} {} {} # {}\n".format(
                        n_dihedrals + ndhid, coeff,
                        abs(int(period)), phase, dih_type))
                    ndhid += 1

            elif period < 0:
                m = dih_type_counts[dih_type]
                coeff = _coeff(divisor, barrier)
                if coeff is not None:
                    out.write("dihedral_coeff {} {} {} {} {}".format(
                        n_dihedrals + ndhid, m, coeff,
                        abs(int(period)), phase))
                    ndhid += 1
                count = 1
                for line_2 in dihe_iter:   # consumes continuation lines
                    p2    = line_2.split()
                    coeff2 = _coeff(float(p2[1]), float(p2[2]))
                    if coeff2 is not None:
                        out.write(" {} {} {}".format(
                            coeff2,
                            abs(int(float(p2[4]))),
                            float(p2[3])))
                    count += 1
                    if int(count) == int(m):
                        break
                out.write(" # {}\n".format(dih_type))

        out.write('\n')

    return new_bond_id - 1, na - 1, ndhid - 1


def update_data_file_with_counts_and_charges(old_data_path, new_data_path,
                                             parm_lammps_path,
                                             n_new_bonds, n_new_angles, n_new_dihedrals,
                                             charges,
                                             n_dummy_atom_types=2,
                                             n_dummy_bond_types=1,
                                             n_dummy_angle_types=1,
                                             n_dummy_dihedral_types=1,
                                             dummy_atom_masses=(12.010, 12.010),
                                             dummy_pair_coeffs=(
                                                 (0.1094, 3.3996695084235338, "c3"),
                                                 (0.1094, 3.3996695084235338, "c6"),
                                             ),
                                             dummy_bond_coeff=(300.9000, 1.5380)):
    """
    Read old_data_path, parse type counts from its header, then:
      1. Update type-count headers explicitly as:
         - atom types     -> +n_dummy_atom_types
         - bond types     -> +(n_new_bonds + n_dummy_bond_types)
         - angle types    -> +(n_new_angles + n_dummy_angle_types)
         - dihedral types -> +(n_new_dihedrals + n_dummy_dihedral_types)
      2. Insert dummy MASS entries for the new atom types using
         `dummy_atom_masses`.
      3. Write corrected partial charges into the Atoms section.
      4. Append explicit dummy pair and bond coefficients to parm_lammps_path
         using `dummy_pair_coeffs` and `dummy_bond_coeff`, while keeping the
         current default dummy angle and dihedral coefficients.

    Returns (new_data_path, parm_lammps_path).
    """
    if n_dummy_atom_types != 2:
        raise ValueError(
            "update_data_file_with_counts_and_charges currently requires "
            "n_dummy_atom_types == 2 because exactly two dummy MASS entries "
            "and pair coefficients are written."
        )
    if len(dummy_atom_masses) != n_dummy_atom_types:
        raise ValueError(
            "dummy_atom_masses must provide one mass per dummy atom type."
        )
    if len(dummy_pair_coeffs) != n_dummy_atom_types:
        raise ValueError(
            "dummy_pair_coeffs must provide one pair coefficient per dummy atom type."
        )
    for coeff in dummy_pair_coeffs:
        if len(coeff) not in (2, 3):
            raise ValueError(
                "Each dummy pair coefficient must be (epsilon, sigma) or "
                "(epsilon, sigma, label)."
            )
    if len(dummy_bond_coeff) != 2:
        raise ValueError("dummy_bond_coeff must contain (k, r0).")

    def _rewrite_count(line, stripped, new_count):
        parts = stripped.split()
        if parts:
            parts[0] = str(new_count)
        return line[: len(line) - len(stripped)] + " ".join(parts) + "\n"

    with open(old_data_path) as f:
        lines = f.readlines()

    ntypes = nbondtypes = nangletypes = ndihtypes = 0
    for line in lines:
        s = line.lstrip()
        if "atom types" in s:
            ntypes = int(s.split()[0])
        elif "bond types" in s:
            nbondtypes = int(s.split()[0])
        elif "angle types" in s:
            nangletypes = int(s.split()[0])
        elif "dihedral types" in s:
            ndihtypes = int(s.split()[0])

    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if "atom types" in stripped:
            lines[i] = _rewrite_count(
                line, stripped, ntypes + n_dummy_atom_types
            )
        elif "bond types" in stripped:
            lines[i] = _rewrite_count(
                line,
                stripped,
                nbondtypes + n_new_bonds + n_dummy_bond_types,
            )
        elif "angle types" in stripped:
            lines[i] = _rewrite_count(
                line,
                stripped,
                nangletypes + n_new_angles + n_dummy_angle_types,
            )
        elif "dihedral types" in stripped:
            lines[i] = _rewrite_count(
                line,
                stripped,
                ndihtypes + n_new_dihedrals + n_dummy_dihedral_types,
            )

    # Locate the Masses section header dynamically — robust to any header length.
    masses_data_start = None
    for i, line in enumerate(lines):
        if line.strip().lower() == "masses":
            masses_data_start = i + 2  # skip header line + blank line
            break
    if masses_data_start is None:
        raise RuntimeError("Masses section not found in data file")

    dummy_masses = [
        f"{ntypes + idx + 1} {mass:.3f}\n"
        for idx, mass in enumerate(dummy_atom_masses)
    ]
    insert_at = masses_data_start + ntypes  # right after the last existing mass entry
    lines[insert_at:insert_at] = dummy_masses

    atoms_start = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("atoms"):
            atoms_start = i + 2
            break
    if atoms_start is None:
        raise RuntimeError("Atoms section not found in data file")
    for idx, charge in enumerate(charges):
        parts = lines[atoms_start + idx].split()
        if len(parts) >= 4:
            parts[3] = f"{charge:.10f}"
            lines[atoms_start + idx] = " ".join(parts) + "\n"

    with open(new_data_path, "w") as f:
        f.writelines(lines)

    with open(parm_lammps_path, "a") as pf:
        pf.write("\n# =============================================================================\n")
        pf.write("# DUMMY PARAMETERS FOR CROSSLINKING\n")
        pf.write("# =============================================================================\n\n")
        for idx, coeff in enumerate(dummy_pair_coeffs):
            epsilon = float(coeff[0])
            sigma = float(coeff[1])
            label = str(coeff[2]) if len(coeff) == 3 else ""
            type_id = ntypes + idx + 1
            line = f"\n pair_coeff {type_id}  {type_id} {epsilon:.4f} {sigma:.16f}"
            if label:
                line += f" # {label}"
            pf.write(line)
        pf.write(
            f"\n bond_coeff {nbondtypes + n_new_bonds + n_dummy_bond_types} "
            f"{dummy_bond_coeff[0]:.4f} {dummy_bond_coeff[1]:.4f}"
        )
        pf.write(
            f"\n angle_coeff {nangletypes + n_new_angles + n_dummy_angle_types} "
            "62.9000 111.5100"
        )
        pf.write(
            f"\n dihedral_coeff {ndihtypes + n_new_dihedrals + n_dummy_dihedral_types} "
            "1 0.0000 3 0.0000"
        )
        pf.write("\n")

    return new_data_path, parm_lammps_path
