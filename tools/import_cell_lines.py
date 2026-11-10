"""Extract the no-DS cell-line cohorts used by the archived transfer analysis."""
import csv
import hashlib
import json
from pathlib import Path
import openpyxl
ROOT=Path(__file__).resolve().parents[1]
name='585A1 and SCTi003-A Cell Counts VWBR- 100mL_For Reke March 10 2025.xlsx'
source=ROOT/'data/source'/name
w=openpyxl.load_workbook(source,data_only=True)
rows=[]
for line,cols,days,rpm,starts,labels in [
    ('585A1',('Q','R','S'),4,50,(6,23,31,39),(5,22,30,38)),
    ('SCTi003-A',('C','D','E'),5,40,(4,23,32,40),(3,22,31,39))]:
    s=w[line]
    assert 'no' in str(s['Q3' if line=='585A1' else 'C2'].value).lower()
    for col in cols:
        replicate=s[f'{col}{labels[0]}'].value
        assert all(s[f'{col}{row}'].value==replicate for row in labels)
        for state,start in zip(('cells','glucose','lactate','aggregate'),starts):
            for day in range(days+1):
                cell=s[f'{col}{start+day}']; value=cell.value
                status='missing' if value is None or value=='No data' else 'measured'
                if status=='measured': assert isinstance(value,(float,int)) and value>=0
                rows.append(dict(cell_line=line,replicate=replicate,volume_ml=100,rpm=rpm,dextran_sulfate='absent',day=day,state=state,value=value if status=='measured' else '',status=status,source_file=name,source_sheet=line,source_cell=cell.coordinate))
out=ROOT/'data/processed/cell_line_measurements.csv'
with out.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
(ROOT/'data/processed/cell_line_source.json').write_text(json.dumps(dict(source_file=name,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),cohort='Workbook no-DS groups, matching original analysis; manuscript Methods +DS statement remains inconsistent'),indent=2)+'\n')
print(f'Extracted {len(rows)} records from six no-DS cell-line vessels.')
