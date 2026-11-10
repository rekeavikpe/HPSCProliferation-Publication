"""Read supplied CFD summaries without changing the source workbooks."""
import csv, math, json
from pathlib import Path
import openpyxl
import hashlib
ROOT=Path(__file__).resolve().parents[1]

OUT=ROOT/'data/processed'

def main():
    OUT.mkdir(exist_ok=True)
    records=[]
    name='Particles_REBUILT_ALLNONAN_fromTXT.xlsx'
    s=openpyxl.load_workbook(ROOT/'data/source'/name,data_only=True)['Particles_REBUILT_ALLNONAN']
    columns={'lag_mean':(4,10),'lag_max':(5,11),'lag_p50':(6,12),'lag_p25':(7,13),'lag_p75':(8,14),'lag_p90':(9,15)}
    expected=('Shear Mean','Shear Maximum','Shear Median','Shear 25th Percentile','Shear 75th Percentile','Shear 90th Percentile','EDR Mean','EDR Maximum','EDR Median','EDR 25th Percentile','EDR 75th Percentile','EDR 90th Percentile')
    assert tuple(s.cell(1,i).value for i in range(4,16))==expected
    for row in range(2,52):
        volume,rpm,size=[s.cell(row,c).value for c in (1,2,3)]
        for scenario,(sc,ec) in columns.items():
            records.append(dict(scenario=scenario,volume_ml=volume,rpm=rpm,size_um=size,shear_pa=s.cell(row,sc).value,edr_m2_s3=s.cell(row,ec).value,source_file=name,source_sheet=s.title,shear_cell=s.cell(row,sc).coordinate,edr_cell=s.cell(row,ec).coordinate))
    name='Summary_of_FluidDomain_Particle_EDR_ShearStress.xlsx';s=openpyxl.load_workbook(ROOT/'data/source'/name,data_only=True)['Fluid_Domain']
    cols={'eul_va':(3,9),'eul_max':(4,10),'eul_p50':(5,11),'eul_p75':(6,12),'eul_p90':(7,13)}
    assert [s.cell(2,c).value for c in range(3,8)]==['Volume Average','Maximum','Median','75th Percentile','90th Percentile']
    volume=None
    for row in range(3,13):
        if s.cell(row,1).value is not None:volume=s.cell(row,1).value
        rpm=s.cell(row,2).value
        for scenario,(sc,ec) in cols.items():
            records.append(dict(scenario=scenario,volume_ml=volume,rpm=rpm,size_um='',shear_pa=s.cell(row,sc).value,edr_m2_s3=s.cell(row,ec).value,source_file=name,source_sheet=s.title,shear_cell=s.cell(row,sc).coordinate,edr_cell=s.cell(row,ec).coordinate))
    for r in records:
        assert all(isinstance(r[k],(int,float)) and math.isfinite(r[k]) and r[k]>=0 for k in ('volume_ml','rpm','shear_pa','edr_m2_s3'))
    keys=[(r['scenario'],r['volume_ml'],r['rpm'],r['size_um']) for r in records];assert len(set(keys))==len(keys)
    previous=list(csv.DictReader((ROOT/'data/processed/lagrangian_exposure.csv').open()))
    for r in previous:
        match=next(x for x in records if x['scenario']==f"lag_p{r['percentile']}" and x['volume_ml']==float(r['volume_ml']) and x['rpm']==float(r['rpm']) and x['size_um']==float(r['size_um']))
        assert match['shear_pa']==float(r['shear_pa']) and match['edr_m2_s3']==float(r['edr_m2_s3'])
    s=openpyxl.load_workbook(ROOT/'data/source/Summary_of_VA_EDR_ShearStress.xlsx',data_only=True)['Sheet1'];volume=None
    for row in range(2,12):
        if s.cell(row,1).value is not None:volume=s.cell(row,1).value
        match=next(x for x in records if x['scenario']=='eul_va' and x['volume_ml']==volume and x['rpm']==s.cell(row,2).value)
        assert math.isclose(match['shear_pa'],s.cell(row,3).value,rel_tol=1e-12)
        assert math.isclose(match['edr_m2_s3'],s.cell(row,4).value,rel_tol=1e-12)
    with (OUT/'exposure_summary_inputs.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    (OUT/'exposure_summary_sources.json').write_text(json.dumps({'source_sha256':{name:hashlib.sha256((ROOT/'data/source'/name).read_bytes()).hexdigest() for name in ('Particles_REBUILT_ALLNONAN_fromTXT.xlsx','Summary_of_FluidDomain_Particle_EDR_ShearStress.xlsx','Summary_of_VA_EDR_ShearStress.xlsx')},'records':len(records),'scenarios':list(columns)+list(cols),'source_policy':'Use rebuilt particle workbook; do not mix older particle summaries into the same percentile series. Eulerian values come from Fluid_Domain; VA checked against separate VA workbook.','excluded_rows':'Fluid_Domain row20 has no volume/RPM identifiers and is not imported. Notes/blank rows excluded.','interpretation':'Workbook-provided marginal shear and EDR summaries; paired maxima need not occur simultaneously. No raw particle trajectories or CFD weighting scheme reconstructed.','lagrangian_size_policy':'Same linear interpolation and endpoint clamping as the main fit.'},indent=2)+'\n')
    print('Verified and imported',len(records),'source-cell-traced records across',len(columns)+len(cols),'scenarios')

if __name__=='__main__':main()
