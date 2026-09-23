#!/usr/bin/env python3
"""Check the Markdown-only documentation policy and frozen evidence checksums."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

from check_manuscript import validate as validate_manuscript

ROOT = Path(__file__).resolve().parents[1]
SKIP = {'.git', '.venv', 'venv', 'env', 'generated', '__pycache__', 'build', 'dist'}
FORBIDDEN = {'.tex', '.bib', '.sty', '.cls', '.aux', '.bbl', '.blg', '.pdf'}


def main() -> None:
    errors = validate_manuscript()
    manifest = json.loads((ROOT / 'validation/evidence_sha256.json').read_text())
    for name, expected in manifest['files'].items():
        path = ROOT / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            errors.append(f'released evidence changed or missing: {name}')
    files = [p for p in ROOT.rglob('*') if p.is_file()
             and not any(part in SKIP or part.endswith('.egg-info') for part in p.relative_to(ROOT).parts)]
    for path in files:
        rel = path.relative_to(ROOT)
        if path.suffix.lower() in FORBIDDEN:
            errors.append(f'unsupported document/build artifact: {rel}')
        if path.suffix == '.txt' and path.name not in {'requirements.txt', 'requirements-tested.txt'}:
            errors.append(f'narrative text should be Markdown: {rel}')
        if path.suffix != '.md':
            continue
        content = path.read_text()
        if re.search(r'\\(?:begin|end|frac|mathrm|mathbb|gamma|lambda)\b|\\[\[\]()]|\$|```math\b', content):
            errors.append(f'use plain-text or Unicode mathematics in Markdown: {rel}')
        # This repository uses inline links without spaces or optional titles.
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', content):
            url = urlsplit(target)
            if url.scheme or not url.path:
                continue
            if not (path.parent / unquote(url.path)).exists():
                errors.append(f'broken local link in {rel}: {target}')
    if errors:
        print('\n'.join(errors), file=sys.stderr)
        raise SystemExit(1)
    docs = sum(p.suffix == '.md' for p in files)
    print(f'PASS: {len(manifest["files"])} evidence checksums; {docs} Markdown documents; '
          'manuscript inventory, figures, table and fragments; local link paths; no TeX/PDF artifacts.')


if __name__ == '__main__':
    main()
