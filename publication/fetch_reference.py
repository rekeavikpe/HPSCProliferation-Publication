"""Fetch Neto's original diffusivity table at a pinned commit and verify its checksum."""
from pathlib import Path
import json,hashlib,urllib.request
root=Path(__file__).resolve().parent.parent/'data/reference/neto2022'
item=next(x for x in json.loads((root/'provenance.json').read_text())['files'] if x['file']=='2PDabPlot.txt')
path=root/item['file']
if path.exists():
 data=path.read_bytes()
else:
 print("Downloading the pinned Neto reference table from the authors' repository.",flush=True)
 data=urllib.request.urlopen(item['url'],timeout=60).read()
if hashlib.sha256(data).hexdigest()!=item['sha256']:raise SystemExit('Neto reference checksum does not match; refusing to use changed data.')
if not path.exists():path.write_bytes(data)
print('PASS: Neto reference table matches the recorded source checksum.',flush=True)
