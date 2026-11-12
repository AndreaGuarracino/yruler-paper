#!/usr/bin/env python3
"""Summarise repeated CV as one mean absolute error per individual."""
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
def read(path):
    with (ROOT / path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))
def write(path, rows, fields):
    with (ROOT / path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
def summaries(data):
    out=[]
    for c in sorted(data):
        by=defaultdict(list)
        for (_,s),v in data[c].items(): by[s].append(abs(v))
        assert by and {len(v) for v in by.values()}=={10}, (c,{len(v) for v in by.values()})
        x=np.array([np.mean(v) for v in by.values()]); out.append({'class':c,'n_males':len(x),'mean_abs':f'{np.mean(x):.5f}','median_abs':f'{np.median(x):.5f}','p90_abs':f'{np.percentile(x,90):.5f}','max_individual_mean_abs':f'{max(x):.5f}','max_prediction_abs':f'{max(abs(v) for v in data[c].values()):.5f}'})
    return out
D=defaultdict(dict)
for r in read('results/model/cv_units.tsv'): D[r['class']][(r['repeat'],r['male'])]=float(r['error_pct'])
D['dyz1:dyz2']={k:100*((1+v/100)/(1+D['dyz2'][k]/100)-1) for k,v in D['dyz1'].items() if k in D['dyz2']}
classes={'dyz3','yq12','dyz1','dyz2','dyz18s','dyz19','dyz17','dyz1:dyz2'}
fields=['class','n_males','mean_abs','median_abs','p90_abs','max_individual_mean_abs','max_prediction_abs']
cv=summaries({k:v for k,v in D.items() if k in classes}); write('results/model/cv_individual_summary.tsv',cv,fields)
ONT=defaultdict(dict)
for r in read('results/model/ont_cv_units.tsv'):
    key=(r['repeat'],r['male'])
    assert key not in ONT[r['class']], (r['class'],key)
    ONT[r['class']][key]=float(r['error_pct'])
ont=summaries({k:v for k,v in ONT.items() if k in classes}); write('results/model/ont_cv_individual_summary.tsv',ont,fields)
out=[]
for depth in ('d1x','d2x','d5x','d10x','d15x','d20x','d25x','d30x','full'):
    path='results/model/depth/depth_'+('f1' if depth=='full' else 'f'+depth)+'.tsv'; dd=defaultdict(dict)
    for r in read(path):
        if r['class'] in classes: dd[r['class']][(r['repeat'],r['sample'])]=float(r['error_pct'])
    for row in summaries(dd): out.append({'depth':depth,**row})
write('results/model/depth_individual_summary.tsv',out,['depth']+fields)
for row in cv: print(row['class'],row['median_abs'],row['p90_abs'])
print('depth d1x:',[(r['class'],r['median_abs']) for r in out if r['depth']=='d1x'])
