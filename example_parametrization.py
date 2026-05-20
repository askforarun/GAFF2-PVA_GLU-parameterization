#!/usr/bin/env python3
"""Example GAFF2 parametrization workflow for PVA-GLU systems.

The script is intentionally explicit so manuscript systems can be regenerated
with a chosen PVA chain length and molecule counts. It builds and parametrizes
PVA, parametrizes the GLU reference, writes corrected ``*_mod`` MOL2/FRCMOD/TOP
files, loads the pre-extracted charge array, and optionally converts a combined
PDB to LAMMPS format.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

# Add current directory to path for imports when run as a script.
REPO_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO_ROOT))

from molecular_utils import getfiles, load_system_charges, run_tleap_with_error_check
from pva_builder import build_pva
from system_constants import GLU_ATOMS_PER_MOLECULE, PVA_ATOMS_PER_MONOMER


def parse_args() -> argparse.Namespace:
    """Parse command-line options for one manuscript system."""
    parser = argparse.ArgumentParser(
        description="Build and parameterize one PVA-GLU system with GAFF2.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--chain-length",
        type=int,
        default=17,
        help="Number of PVA monomers per chain. Manuscript values: 17, 21, 25.",
    )
    parser.add_argument(
        "--n-pva",
        type=int,
        default=2,
        help="Number of PVA chains represented in the charge array / combined PDB.",
    )
    parser.add_argument(
        "--n-glu",
        type=int,
        default=1,
        help="Number of GLU molecules represented in the charge array / combined PDB.",
    )
    parser.add_argument(
        "--combined-pdb",
        type=Path,
        default=None,
        help="Optional combined PDB containing PVA and GLU molecules for LAMMPS conversion.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for LAMMPS outputs when --combined-pdb is provided.",
    )
    return parser.parse_args()


def run_checked(command: list[str]) -> None:
    """Run an external command and fail fast if it returns a non-zero status."""
    subprocess.run(command, check=True)


def apply_atom_type_corrections(input_mol2: Path, output_mol2: Path, replacements: dict[str, str]) -> None:
    """Write a corrected MOL2 file using the same replacements as genhydrogel.py."""
    content = input_mol2.read_text()
    for old, new in replacements.items():
        content = content.replace(old, new)
    output_mol2.write_text(content)


def run_parmchk2(mol2_file: Path, frcmod_file: Path) -> None:
    """Generate an FRCMOD file from a corrected MOL2 file."""
    run_checked(["parmchk2", "-i", str(mol2_file), "-o", str(frcmod_file), "-f", "mol2", "-a", "Y"])


def run_tleap_for_modified_mol2(base: Path) -> None:
    """Generate AMBER topology and coordinates for a corrected MOL2/FRCMOD pair."""
    tleap_input = Path(f"tleap_{base.name}.in")
    tleap_input.write_text(
        "\n".join(
            [
                "source leaprc.gaff2",
                f"MOL = loadmol2 {base}.mol2",
                "check MOL",
                f"loadamberparams {base}.frcmod",
                f"saveamberparm MOL {base}.top {base}.crd",
                "quit",
            ]
        )
    )
    run_tleap_with_error_check(str(tleap_input))


def parameterize_and_correct_pva(chain_length: int) -> tuple[Path, Path]:
    """Build, parameterize, and write corrected PVA files."""
    pva_pdb = Path(f"PVA{chain_length}_trim.pdb")
    atoms, bonds = build_pva(n=chain_length, output_file=str(pva_pdb), cap=False)
    print(f"Built {pva_pdb} with {len(atoms)} atoms and {sum(len(b) for b in bonds.values()) // 2} bonds")

    getfiles(str(pva_pdb))
    pva_base = pva_pdb.with_suffix("")
    pva_mod_base = Path(f"{pva_base}_mod")
    apply_atom_type_corrections(
        pva_base.with_suffix(".mol2"),
        pva_mod_base.with_suffix(".mol2"),
        {"ha": "hc", "c2": "c3"},
    )
    run_parmchk2(pva_mod_base.with_suffix(".mol2"), pva_mod_base.with_suffix(".frcmod"))
    run_tleap_for_modified_mol2(pva_mod_base)
    return pva_pdb, pva_mod_base


def parameterize_and_correct_glu() -> Path:
    """Parameterize the GLU reference and write corrected GLU files."""
    glu_pdb = Path("charge_data/glutaraldehyde.pdb")
    if not glu_pdb.exists():
        raise FileNotFoundError(f"Missing GLU reference structure: {glu_pdb}")

    getfiles(str(glu_pdb))
    glu_base = glu_pdb.with_suffix("")
    glu_mod_base = glu_base.with_name(f"{glu_base.name}_mod")
    apply_atom_type_corrections(
        glu_base.with_suffix(".mol2"),
        glu_mod_base.with_suffix(".mol2"),
        {"c2": "c6", "h4": "h1"},
    )
    run_parmchk2(glu_mod_base.with_suffix(".mol2"), glu_mod_base.with_suffix(".frcmod"))
    run_tleap_for_modified_mol2(glu_mod_base)
    return glu_mod_base


def maybe_convert_to_lammps(args: argparse.Namespace, pva_mod_base: Path, glu_mod_base: Path) -> None:
    """Convert to LAMMPS only when the user provides a combined PDB."""
    output_prefix = args.output_prefix or f"pva{args.chain_length}_glu"
    data_file = f"{output_prefix}.lammps"
    param_file = f"{output_prefix}_parm.lammps"

    if args.combined_pdb is None:
        print("\nLAMMPS conversion skipped: no --combined-pdb was provided.")
        print("To convert later, provide a combined PDB with molecules ordered as:")
        print(f"  {args.n_pva} x PVA{args.chain_length}_trim_mod, then {args.n_glu} x glutaraldehyde_mod")
        print("Example:")
        print(
            "  python example_parametrization.py "
            f"--chain-length {args.chain_length} --n-pva {args.n_pva} --n-glu {args.n_glu} "
            "--combined-pdb combined.pdb"
        )
        return

    if not args.combined_pdb.exists():
        raise FileNotFoundError(f"Combined PDB not found: {args.combined_pdb}")

    from amber_to_lammps import amber2lammps

    amber2lammps(
        data_file=data_file,
        param_file=param_file,
        topologies=[str(pva_mod_base.with_suffix(".top")), str(glu_mod_base.with_suffix(".top"))],
        molecule_counts=[args.n_pva, args.n_glu],
        pdb_file=str(args.combined_pdb),
        charges_target=[0, 0],
        verbose=True,
    )
    print(f"\nLAMMPS conversion complete: {data_file}, {param_file}")


def main() -> bool:
    """Run the parametrization workflow for one selected chain length."""
    args = parse_args()

    print("=" * 70)
    print("PVA-GLU GAFF2 Parametrization Example")
    print("=" * 70)
    print("\nConfiguration:")
    print(f"  Chain length:    {args.chain_length} monomers per PVA chain")
    print(f"  PVA chains:      {args.n_pva}")
    print(f"  GLU molecules:   {args.n_glu}")
    print(f"  Combined PDB:    {args.combined_pdb or 'not provided'}")

    try:
        print("\nSTEP 1-2: Build and parameterize corrected PVA")
        pva_pdb, pva_mod_base = parameterize_and_correct_pva(args.chain_length)

        print("\nSTEP 3: Parameterize corrected GLU reference")
        glu_mod_base = parameterize_and_correct_glu()

        print("\nSTEP 4: Load pre-extracted partial charges")
        charges = load_system_charges(
            chain_length=args.chain_length,
            n_pva=args.n_pva,
            n_glu=args.n_glu,
        )
        pva_total_atoms = args.chain_length * PVA_ATOMS_PER_MONOMER * args.n_pva
        glu_total_atoms = args.n_glu * GLU_ATOMS_PER_MOLECULE
        print(f"Loaded {len(charges)} charges")
        print(f"  PVA atoms: {pva_total_atoms}")
        print(f"  GLU atoms: {glu_total_atoms}")
        print(f"  Net charge: {sum(charges):.5f}")

        print("\nGenerated corrected topology inputs:")
        print(f"  PVA PDB:        {pva_pdb}")
        print(f"  PVA MOL2:       {pva_mod_base.with_suffix('.mol2')}")
        print(f"  PVA FRCMOD:     {pva_mod_base.with_suffix('.frcmod')}")
        print(f"  PVA topology:   {pva_mod_base.with_suffix('.top')}")
        print(f"  GLU MOL2:       {glu_mod_base.with_suffix('.mol2')}")
        print(f"  GLU FRCMOD:     {glu_mod_base.with_suffix('.frcmod')}")
        print(f"  GLU topology:   {glu_mod_base.with_suffix('.top')}")

        print("\nSTEP 5: Optional AMBER-to-LAMMPS conversion")
        maybe_convert_to_lammps(args, pva_mod_base, glu_mod_base)

    except Exception as exc:
        print(f"\nError: {exc}")
        return False

    print("\nPipeline complete.")
    return True


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
