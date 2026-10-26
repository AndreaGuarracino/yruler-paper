# Truth lengths per class and male; families
import collections
def load(O):
    T = collections.defaultdict(dict)
    for l in open(f'{O}/truth_regions.tsv').read().split('\n')[1:]:   # DYZ3 and Yq12 truth arrays
        if l: m, c, src, a, ctg, s, e = l.split('\t'); T[c][m] = float(int(e) - int(s))
    yq = set(T['yq12'])
    h = None
    for l in open(f'{O}/indep_setB.tsv'):   # unit bp from our profile HMMs, assemblies with a truth Yq12
        x = l.rstrip('\n').split('\t')
        if h is None: h = x; continue
        r = dict(zip(h, x)); m = r['assembly'].replace('_R2', '')
        if m not in yq: continue
        T['dyz1'][m] = float(r['dyz1_units']); T['dyz2'][m] = float(r['dyz2_units'])
        T['dyz18s'][m] = float(r['dyz18_units']) + float(r['q3k1_units']); T['yq2k7'][m] = float(r['q2k7_units'])
    h = None
    for l in open(f'{O}/dyz19_clean.tsv'):   # DYZ19 arrays passing flank, gap and flag checks
        x = l.rstrip('\n').split('\t')
        if h is None: h = x; continue
        r = dict(zip(h, x))
        if r['clean'] == 'True': T['dyz19'][r['assembly'].replace('_R2', '')] = float(r['unit_bp_array'])
    return T
def families(ped):
    return {x.split()[1]: x.split()[0] for x in open(ped) if x.strip() and not x.startswith('FamilyID')}
