# Reproduce the checks and figures

Run commands from the repository root. Python code and CSV data remain executable and machine-readable; all narrative documentation and reports use Markdown.

## Environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The cleanup was checked with Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, and Matplotlib 3.10.8. For those numerical package versions use `python -m pip install -r requirements-tested.txt`. That file pins the direct dependencies, not every transitive dependency or BLAS implementation. CI checks the supported Python 3.10 and 3.12 configurations against the dependency ranges in `requirements.txt`.

Matplotlib uses its own math-text renderer for plot labels, with external TeX disabled. No document compiler is required.

## Scientific validation

```bash
python scripts/run_smoke_tests.py
python scripts/final_scientific_validation.py
python scripts/validate_residue_results.py
python scripts/check_repository.py
```

The test suite covers exact identities, API input rejection, basis ordering, numerical entropy behavior, validator failure cases, and figure-output isolation. The comprehensive validator independently enumerates all nonempty supports for $`n = 2`$ and $`n = 4`$, compares against production routines, reads the actual peak table, and checks released residue summaries against the 12,000 observations. The residue command is a focused view of those checks. Default validators print Markdown and do not modify released evidence.

To save a new report:

```bash
python scripts/final_scientific_validation.py --output generated/scientific-validation.md
```

Choose a new Markdown path; report writing rejects existing files and protected evidence destinations. Acceptance checks remain active under `python -O`. A passing report is limited to the checks it lists: it does not reconstruct the high-$`n`$ global searches or establish novelty.

The repository check verifies frozen evidence hashes, local Markdown link paths, and the absence of standalone TeX/PDF artifacts. It does not fetch external links or inspect missing historical records.

## Three figure workflows

| Command | Work performed | Output directory |
|---|---|---|
| `python scripts/reproduce_publication_figures.py` | Redraw all seven figures from released CSVs; no simulations | `generated/figures/` |
| `python scripts/reproduce_publication_figures.py --smoke` | Small fresh numerical runs, including Figure 5 baselines and matched controls; Figure 2 refits the retained table | `generated/smoke/` |
| `python scripts/reproduce_publication_figures.py --full` | Original-size figure computations; expensive; Figure 2 still uses the historical table | `generated/full/` |

The driver rejects destinations that overlap the released `data/` or `outputs/` trees. Missing input CSVs cause failure even if old images are present. Individual scripts also default to generated directories. Explicit custom destinations on individual scripts remain the caller's responsibility.

All outputs are PNG by default. Pixel identity with historical images is not guaranteed across library/font versions. Fresh simulations may differ in floating-point details even with the same seed. CSV redraw demonstrates that a figure can be reconstructed from its supplied inputs; it does not independently regenerate those inputs.

## Local peak check

A modest example is:

```bash
python scripts/peak_scaling_verification.py --n-values 10 --samples 20 --points 7
```

This samples a neighborhood around a supplied Table I value, then fits a quadratic in $`\log_2 M`$. It does not search the entire support-size range. The original global searches for the historical peaks are not part of `--full` either.

Inspect a workload without executing it:

```bash
python scripts/peak_scaling_verification.py --all-table --dry-run
```

The default schedule reaches dimensions that can be expensive with dense exact entropy. A $`d \times d`$ real coefficient matrix needs $`8N`$ bytes before copies, reduced states, or SVD workspace; a complex state vector needs $`16N`$ bytes. At $`n = 30`$ these alone are 8 GiB and 16 GiB. Dense SVD time grows roughly as $`N^{3/2}`$, multiplied by support sizes and sample counts. The September cleanup did not run those high-$`n`$ searches or the full figure simulations.

## Evidence protection

[The manifest](../validation/evidence_sha256.json) records all 50 released CSV/PNG files and the original license. `generated/` is ignored by Git. The manifest's hashes were computed from the pre-cleanup commit and checked against the working tree; they must not be silently regenerated to hide a mismatch.

See [provenance](../PROVENANCE.md) for original-data limitations and [the sanity check](SANITY_CHECK.md) for corrections made during this cleanup.
