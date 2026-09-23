#!/usr/bin/env python3
"""Check manuscript objects and links; semantic review remains a separate task."""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]


def anchors(content: str) -> set[str]:
    """Explicit IDs plus GitHub-style heading IDs used by these documents."""
    prose = re.sub(r'^```.*?^```\s*$', '', content, flags=re.M | re.S)
    result = set(re.findall(r'<a id="([^"]+)"></a>', prose))
    counts: Counter[str] = Counter()
    for heading in re.findall(r'^#{1,6}\s+(.+)$', prose, re.M):
        slug = re.sub(r'[^\w\- ]', '', heading.lower()).replace(' ', '-')
        count = counts[slug]
        result.add(f'{slug}-{count}' if count else slug)
        counts[slug] += 1
    return result


def validate(root: Path = ROOT) -> list[str]:
    errors: list[str] = []
    manifest = json.loads((root / 'validation/manuscript_manifest.json').read_text())
    paper = root / manifest['document']
    if not paper.is_file():
        return ['complete manuscript is missing: PAPER.md']
    content = paper.read_text()
    ids = re.findall(r'<a id="([^"]+)"></a>', content)
    for name, count in Counter(ids).items():
        if count != 1:
            errors.append(f'duplicate manuscript anchor: {name}')
    for name in manifest['required_anchors']:
        if name not in ids:
            errors.append(f'missing manuscript object: {name}')

    counts = manifest['counts']
    if len(re.findall(r'^### ', content, re.M)) != counts['subsections']:
        errors.append('manuscript subsection count differs from the source inventory')
    statements = re.findall(r'^\*\*(Proposition|Theorem|Corollary|Remark) ([\d.]+)', content, re.M)
    expected_statements = [(s['kind'], s['number']) for s in manifest['statements']]
    if [(kind, number.rstrip('.')) for kind, number in statements] != expected_statements:
        errors.append('manuscript numbered statements differ from the source inventory')
    for heading in ['Data availability', 'Acknowledgements', 'Conflict of interest']:
        if f'\n## {heading}\n' not in content:
            errors.append(f'missing manuscript declaration: {heading}')
    equations = re.findall(r'^\*\*\((\d+)\)\*\*$', content, re.M)
    if equations != [str(n) for n in range(1, counts['numbered_equations'] + 1)]:
        errors.append('manuscript equation markers must appear once each in source order')
    references = re.findall(r'^\*\*\[(\d+)\]\*\*', content, re.M)
    if references != [str(n) for n in range(1, counts['references'] + 1)]:
        errors.append('manuscript bibliography entries must appear once each in source order')
    for number in range(1, counts['references'] + 1):
        if f'](#ref-{number})' not in content:
            errors.append(f'uncited manuscript bibliography entry: {number}')
    captions = re.findall(r'^\*\*Figure (\d+)\.\*\*', content, re.M)
    if captions != [str(n) for n in range(1, counts['figures'] + 1)]:
        errors.append('manuscript figure captions must appear once each in source order')
    for figure in manifest['figures']:
        path = root / figure['path']
        if f']({figure["path"]})' not in content:
            errors.append(f'manuscript figure is not embedded: {figure["path"]}')
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != figure['sha256']:
            errors.append(f'imported manuscript figure changed or missing: {figure["path"]}')

    table_rows = []
    for line in content.splitlines():
        match = re.fullmatch(r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|', line)
        if match:
            table_rows.append((int(match[1]), int(match[2]), float(match[3])))
    with (root / 'data/table_i_peaks.csv').open(newline='') as handle:
        expected_rows = [(int(r['n']), int(r['M_n']), float(r['S_n']))
                         for r in csv.DictReader(handle)]
    if table_rows != expected_rows or len(table_rows) != counts['table_data_rows']:
        errors.append('manuscript peak table differs from released CSV or source row count')

    # Check fragment targets in the paper and its supporting repository prose.
    skip = {'.git', '.venv', 'venv', 'env', 'generated', 'build', 'dist'}
    for path in root.rglob('*.md'):
        if any(part in skip or part.endswith('.egg-info') for part in path.relative_to(root).parts):
            continue
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            url = urlsplit(target)
            if url.scheme or not url.fragment:
                continue
            destination = path.parent / unquote(url.path) if url.path else path
            if not destination.is_file() or destination.suffix != '.md':
                continue  # File existence is covered by check_repository.py.
            if unquote(url.fragment) not in anchors(destination.read_text()):
                errors.append(f'broken local fragment in {path.relative_to(root)}: {target}')
    return errors


def main() -> None:
    errors = validate()
    if errors:
        print('\n'.join(errors), file=sys.stderr)
        raise SystemExit(1)
    print('PASS: manuscript inventory, 88 equations, 33 references, 7 figures, '
          '11 table rows, and local fragment links.')


if __name__ == '__main__':
    main()
