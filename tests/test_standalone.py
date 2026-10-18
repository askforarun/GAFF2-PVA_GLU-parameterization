"""Standalone checks; no AmberTools, Packmol, or simulation jobs are run."""
import contextlib
import gzip
import io
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import example_parametrization as example
import genhydrogel
from molecular_utils import REFERENCE_DATA, load_system_charges
from system_constants import pva_atoms_per_chain


class StandaloneTests(unittest.TestCase):
    def test_local_reference_paths(self):
        self.assertEqual(REFERENCE_DATA, ROOT / "charge_data")
        self.assertEqual(genhydrogel.REFERENCE_DATA, ROOT / "charge_data")

    def test_charge_counts_and_cli_reporting(self):
        cases = [(7, 3, 1, 187), (17, 2700, 1350, 371250),
                 (21, 2100, 1050, 347550), (25, 1800, 900, 348300)]
        for n, n_pva, n_glu, count in cases:
            with self.subTest(n=n):
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    charges = load_system_charges(n, n_pva, n_glu)
                self.assertEqual(len(charges), count)
                self.assertAlmostEqual(float(sum(charges)), 0, places=8)
                argv = ["example_parametrization.py", "--chain-length", str(n),
                        "--n-pva", str(n_pva), "--n-glu", str(n_glu)]
                with patch.object(sys, "argv", argv), \
                     patch.object(example, "parameterize_and_correct_pva",
                                  return_value=(Path("PVA.pdb"), Path("PVA_mod"))), \
                     patch.object(example, "parameterize_and_correct_glu",
                                  return_value=Path("GLU_mod")), \
                     patch.object(example, "maybe_convert_to_lammps"), \
                     contextlib.redirect_stdout(output):
                    self.assertTrue(example.main())
                self.assertIn(f"PVA atoms: {n_pva*pva_atoms_per_chain(n)}", output.getvalue())

    def test_production_packed_examples(self):
        for n, n_glu in [(17, 1350), (21, 1050), (25, 900)]:
            path = ROOT / "examples" / f"packed_system_pva{n}_glu{n_glu}.pdb.gz"
            with self.subTest(n=n), gzip.open(path, "rt") as handle:
                records = [line for line in handle if line.startswith(("ATOM  ", "HETATM"))]
                split = 2*n_glu*pva_atoms_per_chain(n)
                self.assertEqual(len(records), split + 31*n_glu)
                def atom_names(source):
                    return [line[12:16].strip() for line in source.read_text().splitlines()
                            if line.startswith(("ATOM  ", "HETATM"))]
                pva_names = atom_names(ROOT / "examples" / f"PVA{n}_trim.pdb")
                glu_names = atom_names(ROOT / "charge_data" / "glutaraldehyde.pdb")
                self.assertEqual([line[12:16].strip() for line in records[:split]],
                                 pva_names * (2*n_glu))
                self.assertEqual([line[12:16].strip() for line in records[split:]],
                                 glu_names * n_glu)


if __name__ == "__main__":
    unittest.main()
