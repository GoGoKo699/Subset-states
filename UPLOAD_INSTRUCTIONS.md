# Upload instructions

This package is designed as an **overlay** for the existing
`GoGoKo699/Subset-states` repository.

1. Create a branch such as `publication-2026`.
2. Upload the contents of this folder into the repository root, allowing files
   with the same names to be replaced.
3. Do **not** delete the existing repository first; this preserves historical
   raw output CSVs that are not duplicated in the overlay.
4. Remove or relocate the paths listed in `FILES_TO_REMOVE.txt`.
5. Run the smoke tests and validation shown in `README.md`.
6. Inspect the branch, then merge it into `main` and create the release tag.

The separate snapshot ZIP is a clean self-contained publication candidate, but
the overlay workflow is preferred for preserving repository history and older
raw outputs.
