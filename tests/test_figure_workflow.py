"""Regression checks for evidence-safe figure reproduction."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "reproduce_publication_figures.py"
SPEC = importlib.util.spec_from_file_location("figure_workflow", SCRIPT)
workflow = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(workflow)


def option(command: list[str], name: str) -> Path:
    return Path(command[command.index(name) + 1])


class FigureWorkflowTests(unittest.TestCase):
    def test_every_mode_separates_outputs_and_fig5_consumes_its_new_data(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "new"
            for mode in ("redraw", "smoke", "full"):
                with self.subTest(mode=mode):
                    commands = workflow.build_commands(mode, destination)
                    for command in commands:
                        flag = "--out" if "--out" in command else "--outdir"
                        self.assertTrue(option(command, flag).is_relative_to(destination))
                    if mode == "redraw":
                        self.assertEqual(len(commands), 7)
                        self.assertEqual(len({Path(command[1]).stem for command in commands}), 7)
                        self.assertTrue(all("from_csv" in command[1] or "fig5_qft_residue_controls" in command[1] for command in commands))
                    else:
                        by_script = {Path(command[1]).name: command for command in commands}
                        plot = by_script["fig5_qft_residue_controls.py"]
                        for producer in ("regenerate_fig5_baseline.py", "run_residue_controls.py"):
                            self.assertEqual(option(by_script[producer], "--outdir"), option(plot, "--data-dir"))
                            self.assertLess(commands.index(by_script[producer]), commands.index(plot))

    def test_existing_pngs_do_not_hide_a_missing_csv(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            released, data, out = root / "released", root / "data", root / "new"
            for number, names in workflow.CSV_INPUTS.items():
                folder = data if number == 5 else released / f"fig{number}"
                folder.mkdir(parents=True, exist_ok=True)
                for name in names:
                    (folder / name).write_text("placeholder\n")
            missing = released / "fig4" / "fig4_zoom_summary.csv"
            missing.unlink()
            for filename in workflow.FIGURE_FILES.values():
                path = out / filename
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"old figure")
            result = subprocess.run([sys.executable, str(SCRIPT), "--released-dir", str(released), "--data-dir", str(data), "--outdir", str(out)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(str(missing), result.stderr)
            self.assertNotIn("+ ", result.stdout)
            self.assertTrue(all((out / name).read_bytes() == b"old figure" for name in workflow.FIGURE_FILES.values()))

    def test_evidence_destination_and_symlink_alias_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            evidence = root / "released"
            evidence.mkdir()
            alias = root / "alias"
            alias.symlink_to(evidence, target_is_directory=True)
            for destination in (root, evidence, evidence / "fig1", alias / "fig1"):
                with self.subTest(destination=destination), self.assertRaises(ValueError):
                    workflow.validate_output_root(destination, evidence)
            workflow.validate_output_root(root / "generated", evidence)

    def test_figure_two_redraw_keeps_input_directory_read_only(self):
        from subset_states.csv_plotting import plot_fig2

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "released"
            source.mkdir()
            name = "fig2_table_with_page.csv"
            shutil.copy2(ROOT / "outputs" / "fig2" / name, source / name)
            before = (source / name).read_bytes()
            destination = root / "generated" / "figure.png"
            plot_fig2(source, destination)
            self.assertEqual(list(source.iterdir()), [source / name])
            self.assertEqual((source / name).read_bytes(), before)
            self.assertTrue(destination.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))


if __name__ == "__main__":
    unittest.main()
