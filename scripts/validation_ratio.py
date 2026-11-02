"""Paired DYZ1:DYZ2 validation; no independent model is fitted for the ratio."""
from contextlib import contextmanager
import numpy as np

KEY='dyz1_dyz2_ratio'
LABEL='DYZ1:DYZ2'
COLOR='#303d43'


def paired(a,b):
    left={r['sample']:r for r in a};right={r['sample']:r for r in b}
    result=[]
    for sample in sorted(left.keys()&right.keys()):   # males with both references (HG01261 has no DYZ2 reference)
        x,y=left[sample],right[sample]
        if 'fold' in x:assert x['fold']==y['fold']
        ref=x['reference_bp']/y['reference_bp']
        factor=(1+x['signed']/100)/(1+y['signed']/100)
        assert np.isfinite(factor) and factor>0
        r=dict(sample=sample,signed=100*(factor-1),reference_ratio=ref,estimate_ratio=ref*factor)
        if 'fold' in x:r['fold']=x['fold']
        if 'chr2_depth' in x:
            assert np.isclose(x['chr2_depth'],y['chr2_depth']);r['chr2_depth']=x['chr2_depth']
        result.append(r)
    return result


def augment(cv,lo,depth):
    cv={**cv,KEY:{r:paired(cv['dyz1'][r],cv['dyz2'][r]) for r in cv['dyz1']}}
    lo={**lo,KEY:{g:paired(lo['dyz1'][g],lo['dyz2'][g]) for g in lo['dyz1']}}
    depth={f:{**groups,KEY:{r:paired(groups['dyz1'][r],groups['dyz2'][r]) for r in groups['dyz1']}}
           for f,groups in depth.items()}
    assert len(cv[KEY])==10 and len({len(rows) for rows in cv[KEY].values()})==1 and len(cv[KEY][0])>0   # same paired men in every repeat
    return cv,lo,depth


@contextmanager
def view(data):
    old=(data.KEYS,data.LABELS,data.COLORS,data.CV,data.LO,data.DEP)
    data.KEYS=[*data.KEYS,KEY];data.LABELS={**data.LABELS,KEY:LABEL};data.COLORS=[*data.COLORS,COLOR]
    data.CV,data.LO,data.DEP=augment(data.CV,data.LO,data.DEP)
    try:yield
    finally:data.KEYS,data.LABELS,data.COLORS,data.CV,data.LO,data.DEP=old


def summary(cv,lo,depth,coverage):
    cv,lo,depth=augment(cv,lo,depth)
    med=[np.median(np.abs([r['signed'] for r in rows])) for rows in cv[KEY].values()]
    pooled=np.abs([r['signed'] for rows in cv[KEY].values() for r in rows])
    errors=np.abs([r['signed'] for rows in lo[KEY].values() for r in rows])
    s=dict(target=KEY,label=LABEL,cv_n=len({r['sample'] for rows in cv[KEY].values() for r in rows}),cv_repeats=10,cv_mean_median_absolute_pct=float(np.mean(med)),
           cv_min_median_absolute_pct=float(min(med)),cv_max_median_absolute_pct=float(max(med)),
           lolo_n=len(errors),lolo_median_absolute_pct=float(np.median(errors)),
           cv_mean_absolute_pct=float(np.mean(pooled)),cv_p90_absolute_pct=float(np.percentile(pooled,90)),
           cv_max_absolute_pct=float(max(pooled)))
    ds=[]
    for fraction,groups in depth.items():
        m=[np.median(np.abs([r['signed'] for r in rows])) for rows in groups[KEY].values()]
        ds.append(dict(fraction=fraction,median_autosomal_depth=coverage[fraction],target=KEY,label=LABEL,
                       repeats=10,mean_median_absolute_pct=float(np.mean(m)),
                       min_median_absolute_pct=float(min(m)),max_median_absolute_pct=float(max(m))))
    return s,ds
