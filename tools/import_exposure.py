"""Extract size-resolved P50/P75 summaries from the verified CFD workbook."""
import csv
import math
from import_data import ROOT, verify_sources
import openpyxl


def main():
    verify_sources(ROOT)
    name = 'Particles_REBUILT_ALLNONAN_fromTXT.xlsx'
    sheet = openpyxl.load_workbook(ROOT / 'data/source' / name, data_only=False)['Particles_REBUILT_ALLNONAN']
    expected = {1:'Volume_mL', 2:'RPM', 3:'Particle_Size_um', 6:'Shear Median',
                8:'Shear 75th Percentile', 12:'EDR Median', 14:'EDR 75th Percentile'}
    for column, label in expected.items():
        if sheet.cell(1,column).value != label:
            raise ValueError(f'Unexpected CFD header in column {column}')
    records=[]
    for row in range(2, sheet.max_row+1):
        volume,rpm,size=(sheet.cell(row,c).value for c in (1,2,3))
        if not all(isinstance(x,(int,float)) for x in (volume,rpm,size)):
            raise ValueError(f'Invalid CFD key at row {row}')
        for percentile,shear_col,edr_col in [(50,6,12),(75,8,14)]:
            tau,edr=sheet.cell(row,shear_col),sheet.cell(row,edr_col)
            if not all(isinstance(x,(int,float)) and math.isfinite(x) and x >= 0 for x in (tau.value,edr.value)):
                raise ValueError(f'Invalid CFD exposure at row {row}')
            records.append(dict(volume_ml=volume,rpm=rpm,size_um=size,percentile=percentile,
                shear_pa=tau.value,edr_m2_s3=edr.value,source_file=name,
                source_sheet=sheet.title,shear_cell=tau.coordinate,edr_cell=edr.coordinate))
    keys=[tuple(r[k] for k in ('volume_ml','rpm','size_um','percentile')) for r in records]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate CFD keys')
    with (ROOT/'data/processed/lagrangian_exposure.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=records[0]);writer.writeheader();writer.writerows(records)
    print(f'Imported {len(records)} CFD percentile records.')


if __name__ == '__main__':
    main()
