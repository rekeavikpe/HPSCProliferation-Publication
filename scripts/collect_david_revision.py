#!/usr/bin/env python3
"""Collect writing tables and figures without optimizer checkpoints or Julia binaries."""
from pathlib import Path
import re
import sys
import zipfile
root = Path(sys.argv[1] if len(sys.argv) > 1 else 'results/david_revision_20261014').resolve()
if not (root / 'COMPLETE').is_file():
    raise SystemExit('The complete suite has not finished. Check progress.log and the current stage log.')
archive = root.with_name(root.name + '_writing.zip')
allowed = {'.csv', '.toml', '.txt', '.md', '.html', '.png', '.pdf', '.svg'}
with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as bundle:
    for file in sorted(root.rglob('*')):
        if not file.is_file():
            continue
        relative = file.relative_to(root)
        if any(re.fullmatch(r'run\d+', part) for part in relative.parts) or 'logs' in relative.parts:
            continue
        if file.suffix.lower() in allowed or file.name == 'COMPLETE':
            bundle.write(file, Path(root.name) / relative)
print(archive)
print(f'{archive.stat().st_size / 1024**2:.1f} MB')
