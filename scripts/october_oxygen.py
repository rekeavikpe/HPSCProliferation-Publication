"""Explicit conditional scaling, not a measured 500 mL coefficient."""
import csv,math,pathlib,sys,json,collections,statistics
source,out=map(pathlib.Path,sys.argv[1:3]);out.mkdir(parents=True,exist_ok=True)
rows=list(csv.DictReader((source/'condition_summary.csv').open()))
anchor=float(next(r['conditional_kla_h'] for r in rows if float(r['volume_ml'])==100 and float(r['rpm'])==30))
ratio=5**(-1/6);scaled=anchor*ratio
q0=17.7;cstar=.21;D=2e-9;packing=1e9 # mol/m3 numerically equals mM
aggregates=collections.defaultdict(list)
for r in csv.DictReader(pathlib.Path('data/processed/measurements.csv').open()):
    if r['state']=='aggregate' and r['status']=='measured': aggregates[(r['condition'],r['day'])].append(float(r['value']))
actual_diameter=max(statistics.mean(v) for v in aggregates.values())
output=[];core=[]
for r in rows:
    volume=float(r['volume_ml']);kla=scaled if volume==500 else float(r['conditional_kla_h'])
    req=float(r['required_kla_h_at_q_reference']);da=req/kla
    output.append(dict(condition=r['source_condition'],volume_ml=volume,rpm=r['rpm'],peak_cells_ml=r['peak_condition_mean_cells_ml'],kla_h=kla,Da=da,kla_sensitivity_min_h=.44,kla_sensitivity_max_h=9.,q_min_amol=8.85,q_max_amol=35.4,Da_sensitivity_min=req*.5/9,Da_sensitivity_max=req*2/.44,source='geometric similarity and constant power number assumed' if volume==500 else 'Neto correlation at 37 C'))
    for qm in [.5,1.,2.]:
        # Packing cells/mL -> cells/m3. Constant consumption sphere, no external film.
        Q=q0*qm*1e-18*packing*1e6
        for diameter,basis in [(336.,'David draft value'),(actual_diameter,'maximum condition-day mean in current processed measurements')]:
            bulk=cstar*(1-da*qm);drop=Q*(diameter*1e-6/2)**2/(6*D)
            core.append(dict(condition=r['source_condition'],q_multiplier=qm,diameter_um=diameter,diameter_basis=basis,packing_cells_ml=packing,diffusivity_m2_s=D,bulk_mM=bulk,drop_mM=drop,center_mM=bulk-drop,anoxic_diameter_um=2*math.sqrt(max(bulk,0)*6*D/Q)*1e6,physical_steady_state=bulk>=0 and bulk-drop>=0))
def write(name,values):
    with (out/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
write('condition_summary.csv',output);write('aggregate_core_sensitivity.csv',core)
exposures=list(csv.DictReader(pathlib.Path('data/processed/exposure_summary_inputs.csv').open()))
eps={int(float(r['volume_ml'])):float(r['edr_m2_s3']) for r in exposures if r['scenario']=='eul_va' and float(r['rpm'])==30 and float(r['volume_ml']) in [100,500]}
write('scaling_scenarios.csv',[
    dict(scenario='constant power number geometric scaling',edr_ratio=5**(2/3),area_per_volume_ratio=5**(-1/3),kla_500_h=scaled,qualification='surface EDR assumed to follow volume scaling'),
    dict(scenario='volume-average CFD EDR proxy sensitivity',edr_ratio=eps[500]/eps[100],area_per_volume_ratio=5**(-1/3),kla_500_h=anchor*(eps[500]/eps[100])**.25*5**(-1/3),qualification='available volume-average EDR is not surface EDR; area scaling still assumed')])

(out/'scaling_assumptions.json').write_text(json.dumps(dict(anchor_100mL_30rpm_h=anchor,volume_ratio=5,kla_ratio=ratio,scaled_500mL_h=scaled,formula='kLa2/kLa1=(epsilon_surface2/epsilon_surface1)^0.25 * (a2/a1)',assumptions=['geometrically similar vessels','same RPM, fluid and gas boundary conditions','constant power number','surface dissipation ratio follows volume-average dissipation scaling V^(2/3)','interfacial area ratio follows V^(2/3)'],verified_actual_500mL_geometry=False,verified_surface_CFD=False,sources=['https://doi.org/10.1002/aic.690160403','https://doi.org/10.3390/bioengineering9080386']),indent=2))
(out/'interpretation.txt').write_text('The 500 mL result is a scaling scenario anchored to the Neto 100 mL correlation. It is not a Neto correlation prediction for 500 mL or an experimental measurement. Actual surface CFD and vessel geometry have not verified these scaling assumptions. The intra-aggregate calculation uses assumed packing and diffusivity, uniform consumption and surface concentration equal to bulk. Negative concentrations invalidate that steady-state regime. A 336 um mean diameter is not the maximum individual diameter. Do not claim that hypoxia is excluded.\n')
(out/'COMPLETE').write_text('Conditional calculations completed. Assumptions remain explicit.\n')
