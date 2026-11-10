"""Extract labeled media-exchange replicates from the collaborator workbook."""
import csv
import hashlib
import json
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
NAME = 'mTesR1_VWBR_100mL_vessel_Data_for_Amon_Lab_March_10_2025.xlsx'
source = ROOT / 'data/source' / NAME
sheet_name = 'Media Exchange Schedule'
sheet = openpyxl.load_workbook(source, data_only=True)[sheet_name]
rows = []
for protocol, columns in ((20, ('S', 'T', 'U')), (50, ('M', 'N', 'O'))):
    for column in columns:
        replicate = sheet[f'{column}3'].value
        assert all(sheet[f'{column}{r}'].value == replicate for r in (23, 33, 43))
        for state, first_row in (('cells', 4), ('glucose', 24), ('lactate', 34), ('aggregate', 44)):
            for day in range(7):
                cell = sheet[f'{column}{first_row + day}']
                value = cell.value
                status = 'measured'
                if state == 'aggregate' and day == 0:
                    status = 'day_zero_placeholder_not_used'
                elif value is None or value == 'No data':
                    status = 'missing'
                else:
                    assert isinstance(value, (float, int)) and value >= 0, cell.coordinate
                rows.append(dict(protocol_percent=protocol, replicate=replicate, day=day,
                                 state=state, value=value if isinstance(value, (float, int)) else '',
                                 status=status, source_file=NAME, source_sheet=sheet_name,
                                 source_cell=cell.coordinate))
out = ROOT / 'data/processed/media_exchange_measurements.csv'
with out.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
manifest = {'source_file': NAME, 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'source_sheet': sheet_name, 'output_sha256': hashlib.sha256(out.read_bytes()).hexdigest(),
            'note': 'Day-zero aggregate zeros are not measured diameters and are excluded; simulation retains the calibration initial-size assumption.'}
(ROOT / 'data/processed/media_exchange_source.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(f'Extracted {len(rows)} records from six labeled bioreactors.')
