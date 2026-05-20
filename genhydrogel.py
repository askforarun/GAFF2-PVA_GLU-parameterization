import os
import shutil
import secrets
import subprocess
from pathlib import Path

import numpy as np


from pva_builder import build_pva
from molecular_utils import (
    generate_monomer_from_smiles,
    getfiles,
    load_system_charges,
    run_tleap_with_error_check,
    print_ff_parameters_crosslinking,
    update_data_file_with_counts_and_charges,
)
from file_operations import cleanup_beginning, cleanup_end
from amber_to_lammps import amber2lammps
from system_constants import GLU_ATOMS_PER_MOLECULE, PVA_ATOMS_PER_MONOMER

REFERENCE_DATA  = Path(__file__).parent / "charge_data"
GLU_PDB         = str(REFERENCE_DATA / "glutaraldehyde.pdb")
CROSSLINK_SMILES = "CCCC1OCCC(CC(O)CC(O)CC(O)CC)O1"
FORCE_FIELD     = "gaff2"


# ---------------------------------------------------------------------------
# Step 2 — build PVA, minimise, trim, pack with GLU
# ---------------------------------------------------------------------------

def build_pva_system(
    n: int,
    n_pva: int,
    n_glu: int,
    init_box_size: float,
    packmol_seed: int | None,
):
    """Build the packed PVA/GLU starting structure with Packmol.

    Inputs
    ------
    n
        Number of PVA monomers per chain.
    n_pva
        Number of PVA chains to place in the box.
    n_glu
        Number of GLU molecules to place in the box.
    init_box_size
        Initial cubic Packmol box size in Angstrom.
    packmol_seed
        Optional Packmol seed. If None, Packmol will choose a random seed.

    Outputs
    -------
    tuple[str, str]
        `(polymer_pdb, crosslinker_pdb)` where:
        - `polymer_pdb` is the trimmed PVA structure file
        - `crosslinker_pdb` is the GLU reference PDB used for packing
    """
    polymer = f"PVA{n}_trim.pdb"
    build_pva(n, output_file=polymer, cap=False)
    if not os.path.exists(polymer):
        raise FileNotFoundError(f"'{polymer}' was not created by build_pva.")

    def _count_pdb_atoms(pdb_path: str) -> int:
        with open(pdb_path, "r") as handle:
            return sum(
                1
                for line in handle
                if line.startswith("ATOM  ") or line.startswith("HETATM")
            )

    expected_total_atoms = (
        n_pva * PVA_ATOMS_PER_MONOMER * n + n_glu * GLU_ATOMS_PER_MOLECULE
    )
    box_size = init_box_size
    maxit = 50
    for attempt in range(1, 6):
        if packmol_seed is None:
            # Packmol expects a positive integer seed. Use a large random value.
            packmol_seed = secrets.randbelow(2_147_483_647) + 1

        with open("packmol.inp", "w") as f:
            f.write(
                "tolerance 2.0\n"
                f"maxit {maxit}\n"
                "output packed_system.pdb\n"
                "filetype pdb\n"
                f"seed {packmol_seed}\n"
            )
            f.write(f"structure {polymer}\nnumber {n_pva}\n")
            f.write(
                f"inside box 0 0 0 {box_size} {box_size} {box_size}\n"
            )
            f.write("end structure\n")
            f.write(f"structure {GLU_PDB}\nnumber {n_glu}\n")
            f.write(
                f"inside box 0 0 0 {box_size} {box_size} {box_size}\n"
            )
            f.write("end structure\n")

        subprocess.run("packmol < packmol.inp > packmol.log", shell=True)
        with open("packmol.log") as f:
            if any("ERROR" in line for line in f):
                raise RuntimeError("Packmol failed. Check packmol.log.")
        if not os.path.exists("packed_system.pdb"):
            raise RuntimeError("packed_system.pdb not created — packmol produced no output.")

        found_atoms = _count_pdb_atoms("packed_system.pdb")
        if found_atoms == expected_total_atoms:
            break

        print(
            f"Packmol placed {found_atoms} / {expected_total_atoms} atoms "
            f"(attempt {attempt}). Increasing box size and retrying."
        )
        box_size = int(box_size * 1.1) + 10
        maxit = min(maxit * 2, 500)
        packmol_seed = None
    else:
        raise RuntimeError(
            "Packmol could not place all molecules after retries. "
            "Check packmol.log and consider increasing init_box_size."
        )

    return polymer, GLU_PDB


# ---------------------------------------------------------------------------
# Step 3 — independent FF parameterisation for PVA and GLU
# ---------------------------------------------------------------------------

def parameterise_molecules(polymer: str, crosslinker: str):
    """Parameterize the packed-molecule building blocks for AMBER.

    Inputs
    ------
    polymer
        Input PVA PDB file to parameterize.
    crosslinker
        Input GLU PDB file to parameterize.

    Outputs
    -------
    tuple[str, str]
        `(pva_mod, glu_mod)` base names without extension for the corrected
        PVA and GLU parameterized files.
    """
    pva_base = os.path.splitext(polymer)[0]
    getfiles(polymer)
    with open(f"{pva_base}.mol2") as f:
        content = f.read().replace("ha", "hc").replace("c2", "c3")
    pva_mod = f"{pva_base}_mod"
    with open(f"{pva_mod}.mol2", "w") as f:
        f.write(content)
    subprocess.run(
        f"parmchk2 -i {pva_mod}.mol2 -o {pva_mod}.frcmod -f mol2 -a Y",
        shell=True,
    )
    print(f"PVA FF parameters → {pva_mod}.frcmod")

    glu_src = Path(crosslinker)
    glu_pdb = glu_src.name
    shutil.copy2(glu_src, glu_pdb)
    glu_base = os.path.splitext(glu_pdb)[0]
    getfiles(glu_pdb)
    with open(f"{glu_base}.mol2") as f:
        content = f.read().replace("c2", "c6").replace("h4", "h1")
    glu_mod = f"{glu_base}_mod"
    with open(f"{glu_mod}.mol2", "w") as f:
        f.write(content)
    subprocess.run(
        f"parmchk2 -i {glu_mod}.mol2 -o {glu_mod}.frcmod -f mol2 -a Y",
        shell=True,
    )
    print(f"GLU FF parameters → {glu_mod}.frcmod")
    return pva_mod, glu_mod


# ---------------------------------------------------------------------------
# Step 4 — AMBER topology + AMBER→LAMMPS conversion
# ---------------------------------------------------------------------------

def convert_to_lammps(pva_mod: str, glu_mod: str, n_pva: int, n_glu: int):
    """Generate AMBER topologies and convert the packed system to LAMMPS.

    Inputs
    ------
    pva_mod
        Base name of the corrected PVA molecule files.
    glu_mod
        Base name of the corrected GLU molecule files.
    n_pva
        Number of PVA chains in the packed system.
    n_glu
        Number of GLU molecules in the packed system.

    Outputs
    -------
    None
        Writes `data.lammps` and `parm.lammps` in the working directory.
    """
    ff = FORCE_FIELD
    for mod, label in ((pva_mod, "PVA"), (glu_mod, "GLU")):
        tleap_in = f"tleap_{label.lower()}.in"
        with open(tleap_in, "w") as f:
            f.write(f"source leaprc.{ff}\n")
            f.write(f"MOL = loadmol2 {mod}.mol2\n")
            f.write("check MOL\n")
            f.write(f"loadamberparams {mod}.frcmod\n")
            f.write(f"saveamberparm MOL {mod}.top {mod}.crd\n")
            f.write("quit")
        run_tleap_with_error_check(tleap_in)
        print(f"{label} topology → {mod}.top")

    amber2lammps(
        data_file="data.lammps",
        param_file="parm.lammps",
        topologies=[f"{pva_mod}.top", f"{glu_mod}.top"],
        molecule_counts=[n_pva, n_glu],
        pdb_file="packed_system.pdb",
        charges_target=[0, 0],
        verbose=True,
    )


# ---------------------------------------------------------------------------
# Step 6 — crosslink reference structure from SMILES
# ---------------------------------------------------------------------------

def build_crosslink_reference(smiles: str = CROSSLINK_SMILES,
                               name: str = "crosslinked_structure") -> str:
    """Build and parameterize the crosslink reference fragment from SMILES.

    Inputs
    ------
    smiles
        SMILES string for the crosslink reference motif.
    name
        Base name for the generated reference files.

    Outputs
    -------
    str
        Base name of the minimized reference structure, for example
        `crosslinked_structure_min`.
    """
    generate_monomer_from_smiles(smiles, name)
    result = subprocess.run(
        f"obabel {name}.pdb -opdb -O {name}_min.pdb --minimize --ff mmff94",
        shell=True, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Crosslink structure minimization failed: {result.stderr}")
    cross_struct = f"{name}_min"
    getfiles(f"{cross_struct}.pdb")
    return cross_struct


# ---------------------------------------------------------------------------
# Step 7 — append FF parameters, add dummy atom types, and update data file
# ---------------------------------------------------------------------------

def prepare_lammps_params(cross_struct: str, charges: np.ndarray) -> str:
    """Finalize the LAMMPS data and parameter files for staged crosslinking.

    Inputs
    ------
    cross_struct
        Base name of the minimized crosslink reference structure.
    charges
        Array of per-atom charges to write into the LAMMPS `Atoms` section.

    Outputs
    -------
    str
        Path to the updated LAMMPS data file, currently
        `data_mod_with_charges.lammps`.

    Notes
    -----
    This step also:
    - appends crosslink force-field terms to `parm.lammps`
    - adds explicit dummy atom, pair, bond, angle, and dihedral type counts
    - inserts the dummy masses and dummy pair/bond coefficients used by the
      staged crosslink workflow
    """
    print("Appending crosslink force-field parameters to parm.lammps...")
    n_new_bonds, n_new_angles, n_new_dihedrals = print_ff_parameters_crosslinking(
        f"{cross_struct}.frcmod", "parm.lammps"
    )
    n_dummy_atom_types = 2
    n_dummy_bond_types = 1
    n_dummy_angle_types = 1
    n_dummy_dihedral_types = 1
    dummy_atom_masses = (12.010, 12.010)
    dummy_pair_coeffs = (
        (0.1094, 3.3996695084235338, "c3"),
        (0.1094, 3.3996695084235338, "c6"),
    )
    # Keep the `c3` / `c6` labels explicit here because update_angles_dihedrals.py
    # later reads these `pair_coeff ... # <label>` comments from parm.lammps to
    # recover atom-type names and apply the built-in c3=c6 chemical equivalence.
    dummy_bond_coeff = (300.9000, 1.5380)
    # Convert the raw AMBER-derived data file into the workflow-ready LAMMPS
    # data file with the following explicit type-count deltas:
    # - atom types     : +2 dummy atom types
    # - bond types     : +(n_new_bonds + 1 dummy bond type)
    # - angle types    : +(n_new_angles + 1 dummy angle type)
    # - dihedral types : +(n_new_dihedrals + 1 dummy dihedral type)
    updated_path, _ = update_data_file_with_counts_and_charges(
        "data.lammps",
        "data_mod_with_charges.lammps",
        "parm.lammps",
        n_new_bonds, n_new_angles, n_new_dihedrals,
        charges,
        n_dummy_atom_types=n_dummy_atom_types,
        n_dummy_bond_types=n_dummy_bond_types,
        n_dummy_angle_types=n_dummy_angle_types,
        n_dummy_dihedral_types=n_dummy_dihedral_types,
        dummy_atom_masses=dummy_atom_masses,
        dummy_pair_coeffs=dummy_pair_coeffs,
        dummy_bond_coeff=dummy_bond_coeff,
    )
    return updated_path


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def generatehydrogel(
    n: int,
    n_pva: int,
    n_glu: int,
    init_box_size: float,
    packmol_seed: int | None = None,
):
    """Generate a PVA-GLU hydrogel system ready for LAMMPS crosslinking.

    Inputs
    ------
    n
        Number of PVA monomers per chain.
    n_pva
        Total number of PVA chains in the system.
    n_glu
        Total number of GLU molecules in the system.
    init_box_size
        Initial cubic Packmol box side length in Angstrom.
    packmol_seed
        Optional Packmol seed. If None, Packmol will choose a random seed.

    Outputs
    -------
    None
        Produces the working-directory files needed by the downstream workflow,
        including `data_mod_with_charges.lammps` and `parm.lammps`.
    """
    cleanup_beginning()

    charges = load_system_charges(n, n_pva, n_glu)
    polymer, crosslinker = build_pva_system(
        n, n_pva, n_glu, init_box_size, packmol_seed
    )
    pva_mod, glu_mod     = parameterise_molecules(polymer, crosslinker)
    convert_to_lammps(pva_mod, glu_mod, n_pva, n_glu)
    cross_struct         = build_crosslink_reference()
    prepare_lammps_params(cross_struct, charges)

    cleanup_end()


if __name__ == "__main__":
    generatehydrogel(n=10, n_pva=20, n_glu=10, init_box_size=100)
