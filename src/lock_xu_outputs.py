from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd

ROOT=Path(r'D:\hill')
OUT=ROOT/'paper2'/'data'/'locked_xu2026'
frames=[]; summaries=[]
for year in (2000,2020):
 for season in ('summer','winter'):
  frames.append(pd.read_csv(OUT/f'patch_thermal_{year}_{season}.csv'))
  summaries.append(pd.read_csv(OUT/f'summary_{year}_{season}.csv'))
patch=pd.concat(frames,ignore_index=True)
summary=pd.concat(summaries,ignore_index=True)
patch.to_csv(OUT/'patch_thermal_all.csv',index=False)
summary.to_csv(OUT/'group_summary_all.csv',index=False)

# Endpoint changes within each season and seasonal contrasts within each year.
wide=summary.pivot_table(index=['dimension','group'],columns=['year','season'],values=['dT_mean_C','dT_median_C','positive_fraction','n'])
rows=[]
for idx,row in wide.iterrows():
 dimension,group=idx
 for season in ('summer','winter'):
  def val(metric,year):
   key=(metric,year,season); return row[key] if key in row.index else np.nan
  rows.append({'comparison':'endpoint_change','dimension':dimension,'group':group,'season':season,
               'value_2000':val('dT_mean_C',2000),'value_2020':val('dT_mean_C',2020),
               'change_C':val('dT_mean_C',2020)-val('dT_mean_C',2000),
               'median_change_C':val('dT_median_C',2020)-val('dT_median_C',2000),
               'positive_fraction_change':val('positive_fraction',2020)-val('positive_fraction',2000)})
 for year in (2000,2020):
  sk=('dT_mean_C',year,'summer'); wk=('dT_mean_C',year,'winter')
  if sk in row.index and wk in row.index:
   rows.append({'comparison':'seasonal_contrast','dimension':dimension,'group':group,'year':year,
                'summer_mean_C':row[sk],'winter_mean_C':row[wk],
                'summer_minus_winter_C':row[sk]-row[wk]})
comparison=pd.DataFrame(rows)
comparison.to_csv(OUT/'locked_comparisons.csv',index=False)

audit=[]
for (year,season),d in patch.groupby(['year','season']):
 valid=d.dT_C.notna()
 audit.append({'year':year,'season':season,'n_patches':len(d),'n_valid':int(valid.sum()),
               'success_pct':100*valid.mean(),'mean_T1_valid_fraction':d.T1_valid_fraction.mean(),
               'pct_T1_full_support':100*(d.T1_valid_fraction>=.999999).mean(),
               'median_T0_pixels':d.T0_pixels.median(),'min_T0_pixels':d.T0_pixels.min(),
               'mean_dT_C':d.dT_C.mean(),'median_dT_C':d.dT_C.median()})
audit=pd.DataFrame(audit)
audit.to_csv(OUT/'quality_audit.csv',index=False)

inputs=[]
for year in (2000,2020):
 for season in ('summer','winter'):
  for role,path in [('patches',ROOT/'数据输出'/'03_Results'/'patches'/f'is_patches_{year}.shp'),('lst',ROOT/'数据输出'/'04_LST'/f'LST_{year}_{season}_albers.tif'),('dem',ROOT/'数据输出'/'01_DEM_Slope'/'QTP_DEM_30m.tif'),('impervious',ROOT/'数据输出'/'02_LandCover'/f'QTP_IS_{year}.tif')]:
   stat=path.stat(); inputs.append({'year':year,'season':season,'role':role,'path':str(path),'size_bytes':stat.st_size,'modified':stat.st_mtime})
pd.DataFrame(inputs).drop_duplicates(['year','season','role']).to_csv(OUT/'input_manifest.csv',index=False)
config={'reference':'Xu et al. (2026), Habitat International 171, 103778','years':[2000,2020],
        'seasons':['summer','winter'],'slope_classes':['0-5','5-10','10-15','15-20','20-30','30-plus'],
        'hillside_threshold_deg':5,'ring_count':10,'ring_area_ratio':0.5,'t1_area_multiplier':6.0,
        'elevation_tolerance_m':0,'lst_support_m':100,'dem_support_m':30,
        't0_domain':'QTP-wide same-elevation other land after excluding impervious pixels and all T1 footprints',
        'paper1_primary_outputs_used':False}
(OUT/'analysis_lock.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
print('QUALITY')
print(audit.to_string(index=False))
print('\nELEVATION SUMMARY')
print(summary[summary.dimension=='elevation'].sort_values(['group','year','season']).to_string(index=False))
print('\nELEVATION COMPARISONS')
print(comparison[comparison.dimension=='elevation'].to_string(index=False))
print('\nSLOPE COMPARISONS')
print(comparison[comparison.dimension=='slope'].to_string(index=False))
