# Truth lengths per class and male; families
import os, collections
def load(O):
    T = collections.defaultdict(dict); A12 = {}
    for l in open(f'{O}/truth_regions.tsv').read().split('\n')[1:]:   # DYZ3 and Yq12 truth arrays
        if l: m, c, src, a, ctg, s, e = l.split('\t'); T[c][m] = float(int(e) - int(s)); A12.update({m: a} if c == 'yq12' else {})
    yq = set(T['yq12'])
    h = None
    for l in open(f'{O}/indep_setB.tsv'):   # Hallast annotation table (um_halregions.py)
        x = l.rstrip('\n').split('\t')
        if h is None: h = x; continue
        r = dict(zip(h, x)); m = r['assembly'].replace('_R2', '')
        if m not in yq or r['assembly'] != A12[m]: continue   # same assembly as the Yq12 truth
        T['dyz1'][m] = float(r['dyz1_units']); T['dyz2'][m] = float(r['dyz2_units'])
        T['dyz18s'][m] = float(r['dyz18_units']) + float(r['q3k1_units']); T['yq2k7'][m] = float(r['q2k7_units'])
        T['dyz18sub'][m] = T['dyz18s'][m] + T['yq2k7'][m]   # whole DYZ18 subregion
        if r.get('q2k7a_units'): T['yq2k7a'][m] = float(r['q2k7a_units']); T['dyz18suba'][m] = T['dyz18s'][m] + T['yq2k7a'][m]   # amount variants
        T['dyz18'][m] = float(r['dyz18_units']); T['yq3k1'][m] = float(r['q3k1_units'])
    h = None
    for l in open(f'{O}/dyz19_clean.tsv'):   # DYZ19 arrays passing flank, gap and flag checks
        x = l.rstrip('\n').split('\t')
        if h is None: h = x; continue
        r = dict(zip(h, x))
        if r['clean'] == 'True': T['dyz19'][r['assembly'].replace('_R2', '')] = float(r['unit_bp_array'])
    f = os.environ.get('D17T', f'{os.path.dirname(O)}/newclass/truth_all.tsv'); c17 = os.environ.get('D17C', 'dyz17_span')   # DYZ17: D17T table (release: um_halregions.py)
    if os.path.exists(f):
        h = None
        for l in open(f):
            x = l.rstrip('\n').split('\t')
            if h is None: h = x; continue
            r = dict(zip(h, x)); T['dyz17'][r['sample'].replace('_R2', '')] = float(r[c17])
            if r.get('dyz17_amount') and c17 != 'dyz17_amount': T['dyz17a'][r['sample'].replace('_R2', '')] = float(r['dyz17_amount'])   # amount variant
    return T
def families(ped):   # 1kGP pedigree plus the GIAB trios
    F = {x.split()[1]: x.split()[0] for x in open(ped) if x.strip() and not x.startswith('FamilyID')}
    return F | {s: 'GIAB_AJ' for s in ('HG002', 'HG003', 'HG004')} | {s: 'GIAB_HAN' for s in ('HG005', 'HG006', 'HG007')}
