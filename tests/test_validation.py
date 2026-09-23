"""Adversarial checks for scientific validators and evidence preservation."""
from __future__ import annotations

import contextlib
import csv
import io
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

from scripts import final_scientific_validation as scientific
from scripts import validate_residue_results as residue

ROOT = Path(__file__).resolve().parents[1]


class ValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def copy_file(self, relative: str) -> Path:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, path)
        return path

    def table_fixture(self) -> Path:
        path = self.copy_file("data/table_i_peaks.csv")
        self.copy_file("outputs/fig2/fig2_table_with_page.csv")
        self.copy_file("outputs/fig2/fig2_linear_fit_summary.csv")
        return path

    def residue_fixture(self) -> None:
        for name in ("residue_matched_summary.csv", "residue_matched_samples.csv", "residue_entropy_bounds.csv"):
            self.copy_file("data/" + name)

    @staticmethod
    def mutate(path: Path, mutation) -> None:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            fields, rows = reader.fieldnames, list(reader)
        mutation(rows)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    def test_retained_table_and_derived_csv_agree(self) -> None:
        self.table_fixture()
        result = scientific.table_fit_checks(self.root)
        self.assertEqual(result["table_rows_checked"], 11)

    def test_changed_table_value_is_not_hidden_by_constants(self) -> None:
        path = self.table_fixture()
        self.mutate(path, lambda rows: rows[0].update(S_n="4.073"))
        with self.assertRaisesRegex(residue.ValidationError, "CSV and production defaults"):
            scientific.table_fit_checks(self.root)

    def test_changed_derived_fit_is_rejected(self) -> None:
        self.table_fixture()
        path = self.root / "outputs/fig2/fig2_linear_fit_summary.csv"
        self.mutate(path, lambda rows: rows[0].update(slope="0.75"))
        with self.assertRaisesRegex(residue.ValidationError, "stored log2_M_n fit slope"):
            scientific.table_fit_checks(self.root)

    def test_zero_cannot_replace_tiny_fit_p_value(self) -> None:
        self.table_fixture()
        path = self.root / "outputs/fig2/fig2_linear_fit_summary.csv"
        self.mutate(path, lambda rows: rows[0].update(p_value="0"))
        with self.assertRaisesRegex(residue.ValidationError, "p_value"):
            scientific.table_fit_checks(self.root)

    def test_exhaustive_check_calls_production_formula(self) -> None:
        result = scientific.exhaustive_checks((2,))
        self.assertEqual(result["supports_enumerated"], 15)
        with mock.patch.object(scientific.core, "fixed_cardinality_average_purity", return_value=0.5):
            with self.assertRaisesRegex(residue.ValidationError, "production average purity"):
                scientific.exhaustive_checks((2,))

    def test_nonfinite_production_output_is_rejected(self) -> None:
        with mock.patch.object(scientific.core, "entropy_from_support", return_value=float("nan")):
            with self.assertRaisesRegex(residue.ValidationError, "non-finite"):
                scientific.exhaustive_checks((2,))

    def test_all_recorded_residue_statistics_and_seed_prefix(self) -> None:
        self.residue_fixture()
        result = residue.residue_checks(self.root, replay_samples=1)
        self.assertEqual(result["raw_rows_checked"], 12000)
        self.assertEqual(result["seed_replayed_supports"], 12)

    def test_nonfinite_raw_entropy_is_rejected(self) -> None:
        self.residue_fixture()
        path = self.root / "data/residue_matched_samples.csv"
        self.mutate(path, lambda rows: rows[0].update(position_entropy="nan"))
        with self.assertRaisesRegex(residue.ValidationError, "non-finite"):
            residue.residue_checks(self.root, replay_samples=0)

    def test_duplicate_sample_cannot_silently_replace_a_missing_sample(self) -> None:
        self.residue_fixture()
        path = self.root / "data/residue_matched_samples.csv"
        self.mutate(path, lambda rows: rows[1].update(sample=rows[0]["sample"]))
        with self.assertRaisesRegex(residue.ValidationError, "duplicate or missing sample"):
            residue.residue_checks(self.root, replay_samples=0)

    def test_missing_group_cannot_shrink_validation_scope(self) -> None:
        self.residue_fixture()
        path = self.root / "data/residue_matched_summary.csv"
        self.mutate(path, lambda rows: rows.pop())
        with self.assertRaisesRegex(residue.ValidationError, "12 recorded groups"):
            residue.residue_checks(self.root, replay_samples=0)

    def test_summary_quantile_is_reconstructed_from_raw_samples(self) -> None:
        self.residue_fixture()
        path = self.root / "data/residue_matched_summary.csv"
        self.mutate(path, lambda rows: rows[0].update(null_position_q025="5.0"))
        with self.assertRaisesRegex(residue.ValidationError, "null_position_q025"):
            residue.residue_checks(self.root, replay_samples=0)

    def test_replay_uses_production_sampler(self) -> None:
        self.residue_fixture()
        with mock.patch.object(residue, "sample_matched_support", return_value=np.array([0])):
            with self.assertRaisesRegex(residue.ValidationError, "sampler cardinality"):
                residue.residue_checks(self.root, replay_samples=1)

    def test_optimized_cli_still_rejects_corrupt_data(self) -> None:
        path = self.table_fixture()
        self.mutate(path, lambda rows: rows[0].update(S_n="nan"))
        code = (
            "from pathlib import Path; import sys; "
            "from scripts import final_scientific_validation as m; "
            "original=m.table_fit_checks; root=Path(sys.argv[1]); "
            "m.table_fit_checks=lambda: original(root); sys.argv=['validator']; "
            "raise SystemExit(m.main())"
        )
        for flag in ("-O", "-OO"):
            with self.subTest(flag=flag):
                process = subprocess.run([sys.executable, flag, "-c", code, str(self.root)],
                    cwd=ROOT, env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
                    capture_output=True, text=True, check=False)
                self.assertEqual(process.returncode, 1, process.stderr)
                self.assertIn("Status: FAIL", process.stderr)
                self.assertIn("non-finite", process.stderr)

    def test_report_defaults_to_stdout_and_never_overwrites(self) -> None:
        report = "# Test report\n\nChecked.\n"
        with contextlib.redirect_stdout(io.StringIO()) as output:
            residue.write_report(report, None)
        self.assertEqual(output.getvalue().strip(), report.strip())
        path = self.root / "report.md"
        with contextlib.redirect_stdout(io.StringIO()):
            residue.write_report(report, path)
        self.assertEqual(path.read_text(), report)
        with self.assertRaisesRegex(residue.ValidationError, "overwrite"):
            residue.write_report("changed", path)
        self.assertEqual(path.read_text(), report)
        with self.assertRaisesRegex(residue.ValidationError, "Markdown"):
            residue.write_report(report, self.root / "report.txt")
        with self.assertRaisesRegex(residue.ValidationError, "evidence directories"):
            residue.write_report(report, ROOT / "data/new-validation-report.md")


if __name__ == "__main__":
    unittest.main()
