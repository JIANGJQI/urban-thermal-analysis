from pathlib import Path
import pickle,time
import numpy as np
import pandas as pd
import shapely
from scipy.stats import pearsonr
ROOT=Path(r'D:\hill'); OUT=ROOT/'paper2'/'data'/'locked_xu2026'
patch=pd.read_csv(OUT/'patch_thermal_all.csv')
land=pd.read_csv(OUT/'landscape_indices_2000_2020.csv')
label_map={'[0-5)':'0-5','[5-10)':'5-10','[10-15)':'10-15','[15-20)':'15-20','[20-30)':'20-30','[30+]':'30-plus'}
land['slope_label']=land.slope_class.map(label_map)
metrics=[]
for year in (2000,2020):
 with (OUT/f'footprints_{year}.pkl').open('rb') as f: footprints=pickle.load(f)
 for season in ('summer','winter'):
  d=patch[(patch.year==year)&(patch.season==season)&patch.dT_C.notna()].copy()
  for code,label in enumerate(('0-5','5-10','10-15','15-20','20-30','30-plus')):
   q=d[d.slope_class==code]
   positive=q[q.dT_C>0]
   geoms=[footprints[int(i)] for i in positive.source_index]
   started=time.perf_counter()
   area=shapely.area(shapely.union_all(geoms))/1_000_000 if geoms else 0.0
   metrics.append({'year':year,'season':season,'slope_label':label,'n':len(q),
    'TMax_C':q.dT_C.max(),'TMean_C':q.dT_C.mean(),'TArea_union_km2':area,
    'n_positive':len(positive),'positive_fraction':(q.dT_C>0).mean()})
   print(year,season,label,len(geoms),f'{area:.3f}',f'{time.perf_counter()-started:.1f}s',flush=True)
thermal=pd.DataFrame(metrics); thermal.to_csv(OUT/'warming_metrics_by_slope.csv',index=False)
# Xu-style static and temporal-change correlations across six slope classes.
rows=[]
indices=['PLAND','NP','AREA_MN','LSI','COHESION','AI']; outcomes=['TMax_C','TMean_C','TArea_union_km2']
for season in ('summer','winter'):
 for year in (2000,2020):
  merged=land[land.year==year].merge(thermal[(thermal.year==year)&(thermal.season==season)],on='slope_label')
  for index in indices:
   for outcome in outcomes:
    r,p=pearsonr(merged[index],merged[outcome]); rows.append({'analysis':'static','season':season,'year':year,'index':index,'outcome':outcome,'r':r,'R2':r*r,'p':p,'n':len(merged)})
 # Endpoint change correlations, mirroring the source article.
 l0=land[land.year==2000].set_index('slope_label'); l1=land[land.year==2020].set_index('slope_label')
 t0=thermal[(thermal.year==2000)&(thermal.season==season)].set_index('slope_label'); t1=thermal[(thermal.year==2020)&(thermal.season==season)].set_index('slope_label')
 for index in indices:
  x=l1[index]-l0[index]
  for outcome in outcomes:
   y=t1[outcome]-t0[outcome]; r,p=pearsonr(x,y); rows.append({'analysis':'endpoint_change','season':season,'year':np.nan,'index':index,'outcome':outcome,'r':r,'R2':r*r,'p':p,'n':len(x)})
pd.DataFrame(rows).to_csv(OUT/'pearson_landscape_warming.csv',index=False)
print('\nTHERMAL'); print(thermal.to_string(index=False))
print('\nP<.05'); print(pd.DataFrame(rows).query('p < .05').to_string(index=False))
