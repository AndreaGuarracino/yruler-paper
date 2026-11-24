# Summaries of the CV, lineage-out, batch holdout and depth runs, and the tables for the figure scripts
# usage: (TRUTHX as the model) python3 summary.py cv|lolo|holdout <base folds dir> <floor folds dir or none> <out.tsv> (dirs relative to $P/$UM)
#        python3 summary.py depth <dir with cv_depth_<label>.tsv> <out.tsv> <label> ... | export <v1 dir> <out dir> (CHR8DEPTH=<chr8 window 31-mer depth per level and sample>)
# reported classes: floors for dyz2 and yq2k7, Yq12 = DYZ1 + DYZ2 + DYZ18S + rest; errors in % of the assembly length
import os, sys, glob, collections, numpy as np
if sys.argv[1] == 'export':   # per prediction: repeat, fold, class, sample, choice, markers, error, truth bp, estimate bp, chr8 window depth of the test run; level 1 = the main CV (cv/), other levels = depth/
    import csv
    V, OD = sys.argv[2:4]; os.makedirs(OD, exist_ok=True); CL = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17']; LV = {'1': 'full', 'd30x': 'd30x', 'd25x': 'd25x', 'd20x': 'd20x', 'd15x': 'd15x', 'd10x': 'd10x', 'd5x': 'd5x', 'd2x': 'd2x', 'd1x': 'd1x'}
    def log(*a): print(*a, flush=True)
    def units(f):
        R = {}
        for l in list(open(f))[1:]:
            c, r, m, e = l.rstrip('\n').split('\t')
            if c in CL: R[(c, r, m)] = e
        return R
    def rows(pat, key):   # (class, key fields, sample) -> row of the per-class fold files; DYZ2 from the G+C floor run
        R = {}
        for c in CL:
            d = pat(c == 'dyz2'); fs = glob.glob(f'{d}/*_c{c}.tsv')
            if not fs: sys.exit(f'ERROR no files {d}/*_c{c}.tsv')
            for f in fs:
                for l in open(f):
                    x = l.rstrip('\n').split('\t'); R[(c, key(x), x[3] if key(x) is not None else x[2])] = x
        return R
    DEP = {(r['level'], r['sample']): r['chr8_depth_mean'] for r in csv.DictReader(open(os.environ['CHR8DEPTH']), delimiter='\t')}
    def dep(lv, s):   # chr8 window 31-mer depth of the test run; full depth when the sample has no run at that level (the CV used its full-depth counts)
        if (lv, s) not in DEP: log(f'{s}: no {lv} chr8 depth, full depth used')
        return DEP[(lv, s)] if (lv, s) in DEP else DEP[('full', s)]
    for f, lv in LV.items():
        E = units(f'{V}/cv/cv_units.tsv' if f == '1' else f'{V}/depth/cv_depth_{lv}_units.tsv')
        X = rows(lambda fl: f'{V}/cv/folds{"_dyz2_gcfloor" * fl}' if f == '1' else f'{V}/depth/folds{"_dyz2_gcfloor" * fl}/{lv}', lambda x: x[0])
        n = collections.Counter(); miss = 0
        with open(f'{OD}/depth_f{f}.tsv', 'w') as o:
            o.write('repeat\tfold\tclass\tsample\tchoice_presence_outside_ruler_filter\tmarkers\terror_pct\treference_bp\testimate_bp\tchr8_depth_test\n')
            for (c, r, m), e in sorted(E.items(), key=lambda z: (CL.index(z[0][0]), int(z[0][1]), z[0][2])):
                x = X.get((c, r, m))
                if x is None: miss += 1; continue
                t = float(x[7]); o.write('\t'.join([r, x[1], c, m, x[4], x[5], e, f'{t:.1f}', f'{t * (1 + float(e) / 100):.0f}', dep(lv, m)]) + '\n'); n[c] += 1
        if miss: sys.exit(f'ERROR level {lv}: {miss} predictions without a CV row')
        log(f'level {f} ({lv}): ' + ', '.join(f'{c} {n[c]}' for c in CL) + f'; men {len({k[2] for k in E})}')
    E = units(f'{V}/lineage_out/lolo_units.tsv'); X = rows(lambda fl: f'{V}/lineage_out/folds{"_dyz2_gcfloor" * fl}', lambda x: None)
    with open(f'{OD}/lineage_out.tsv', 'w') as o:
        o.write('lineage\tclass\tsample\tchoice_presence_outside_ruler_filter\tmarkers\terror_pct\treference_bp\testimate_bp\n'); k = 0
        for (c, r, m), e in sorted(E.items(), key=lambda z: (CL.index(z[0][0]), z[0][2])):
            x = X.get((c, None, m))
            if x is None: sys.exit(f'ERROR lineage-out {c} {m}: no row')
            t = float(x[6]); o.write('\t'.join([x[0], c, m, x[3], x[4], e, f'{t:.1f}', f'{t * (1 + float(e) / 100):.0f}']) + '\n'); k += 1
    log(f'lineage_out: {k} rows, lineages {sorted({x[0] for x in X.values()})}')
    sys.exit(0)
if sys.argv[1] == 'depth':   # depth <dir> <out.tsv> <label> ...: one table over the depth levels
    M = sys.argv[2]
    with open(sys.argv[3], 'w') as o:
        o.write('depth\tclass\tn_predictions\tn_males\tmean_abs\tmedian_abs\tp90_abs\tmax_abs\tbias_median\trepeat_median_mean\tfemale_median_bp\tfemale_max_bp\n')
        for l in sys.argv[4:]:
            R = [x.rstrip('\n').split('\t') for x in open(f'{M}/cv_depth_{l}.tsv')]; F = {x.split('\t')[0]: x.rstrip('\n').split('\t') for x in list(open(f'{M}/cv_depth_{l}_females.tsv'))[1:]}
            for r in R[1:]: o.write('\t'.join([l] + r[:8] + [r[9], F.get(r[0], ['.'] * 5)[2], F.get(r[0], ['.'] * 5)[4]]) + '\n')
    print(open(sys.argv[3]).read(), flush=True); sys.exit(0)
P = os.environ['P']; O = f'{P}/{os.environ["UM"]}'; IN = os.environ['IN']; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import truth
K, B, F, OUT = sys.argv[1], f'{O}/{sys.argv[2]}', f'{O}/{sys.argv[3]}', sys.argv[4]
def log(*a): print(*a, flush=True)
T = truth.load(f'{IN}/truth')
if os.environ.get('TRUTHX'):
    for l in open(os.environ['TRUTHX']): c, s = l.split()[:2]; T[c].pop(s, None)
G = {}; TRU = {}; RULE = os.environ.get('RULE', '')   # (class, unit) -> training group: (repeat, fold) for cv, lineage for lolo; RULE=r<presence>_<outside>: rows of one fixed rule (cv_eval.py allrules)
def read(d):   # (class, unit) -> % error; unit = (repeat, male) for cv, male otherwise
    E = {}
    for fn in sorted(glob.glob(f'{d}/*.tsv')):
        if fn.endswith('_fem.tsv'): continue
        for l in open(fn):
            x = l.rstrip('\n').split('\t')
            if K == 'cv': rep, k, c, s, e = x[0], x[1], x[2], x[3], x[6]; u = (rep, s)
            elif K == 'lolo': c, s, e = x[1], x[2], x[5]; u = s
            else: c, s, e = x[:3]; u = s   # holdout: class, male, error[, truth]
            if (c.split('|')[1] if '|' in c else '') != RULE: continue
            c = c.split('|')[0]; E[(c, u)] = float(e); G[(c, u)] = (x[0], x[1]) if K == 'cv' else x[0] if K == 'lolo' else '.'
            if len(x) > {'cv': 7, 'lolo': 6}.get(K, 3): TRU[(c, u)] = float(x[-1])   # truth used by cv_eval.py (per-fold all-copies truth)
    return E
EB, EF = read(B), read(F) if os.path.isdir(F) else {}
if os.environ.get('TRUTHNEW'):   # same predictions scored against another truth
    TN = {(l.split()[0], l.split()[1]): float(l.split()[2]) for l in open(os.environ['TRUTHNEW'])}
    for E in (EB, EF):
        for (c, u), e in E.items():
            s = u[1] if K == 'cv' else u
            if (c, s) in TN and s in T[c]: E[(c, u)] = 100 * (T[c][s] * (1 + e / 100) / TN[(c, s)] - 1)
    for (c, s), v in TN.items():
        if s in T[c]: T[c][s] = v
    log('TRUTHNEW', os.environ['TRUTHNEW'], len(TN), 'truths')
log(K, 'base rows', len(EB), 'floor rows', len(EF))
R = collections.defaultdict(dict)
FC = sorted({c for c, u in EF})   # classes with G+C floors: those in the floor dir
for (c, u), e in EB.items():
    if c != 'yq12' and c not in FC: R[c][u] = e
for c in FC: R[c] = {u: e for (cc, u), e in EF.items() if cc == c}
EX = {l.split()[0] for l in open(os.environ['EXCL'])} if os.environ.get('EXCL') else set()
RS = {m: T['yq12'][m] - T['dyz1'][m] - T['dyz2'][m] - T['dyz18s'][m] for m in T['yq12'] if all(m in T[c] for c in ('dyz1', 'dyz2', 'dyz18s'))}   # Yq12 rest per truth male
D2 = EF if EF else EB; R['yq12'] = {}
for (c, u), e1 in EB.items():   # Yq12 = DYZ1 + DYZ2 + DYZ18S + rest (class yq12rest2), the only definition
    s = u[1] if K == 'cv' else u
    if c != 'dyz1' or ('dyz2', u) not in D2 or ('yq12rest2', u) not in EB or ('dyz18s', u) not in EB or ('yq12', u) not in EB: continue
    t12, t1, t2, t18, tr = (TRU[(c0, u)] if (c0, u) in TRU else T[c0][s] for c0 in ('yq12', 'dyz1', 'dyz2', 'dyz18s', 'yq12rest2')); R['yq12'][u] = 100 * ((t1 * (1 + e1 / 100) + t2 * (1 + D2[('dyz2', u)] / 100) + t18 * (1 + EB[('dyz18s', u)] / 100) + tr * (1 + EB[('yq12rest2', u)] / 100)) / t12 - 1)
log(f'Yq12 = DYZ1 + DYZ2 + DYZ18S + rest: {len(R["yq12"])} predictions, rest median of all truth males {np.median(list(RS.values())):.0f} bp')
CL = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'yq2k7', 'dyz18sub', 'dyz19', 'dyz17', 'yq2k7a', 'dyz18suba', 'dyz17a']   # last 3: amount variants
with open(OUT[:-4] + '_units.tsv', 'w') as o:   # per prediction, for paired comparisons
    o.write('class\trepeat\tmale\terror_pct\n')
    for c in CL:
        for u, e in sorted(R[c].items()): o.write(f'{c}\t{u[0] if K == "cv" else "."}\t{u[1] if K == "cv" else u}\t{e:.4f}\n')
with open(OUT, 'w') as o:
    o.write('class\tn_predictions\tn_males\tmean_abs\tmedian_abs\tp90_abs\tmax_abs\tbias_median\tworst_male' + ('\trepeat_median_mean\trepeat_median_min\trepeat_median_max' if K == 'cv' else '') + '\n')
    for c in CL:
        if not R[c]: log(c, 'no rows'); continue
        a = np.abs(np.array(list(R[c].values()))); us = list(R[c]); w = us[int(a.argmax())]; ms = {u[1] if K == 'cv' else u for u in us}
        row = [c, len(a), len(ms), f'{a.mean():.2f}', f'{np.median(a):.2f}', f'{np.percentile(a, 90):.2f}', f'{a.max():.2f}', f'{np.median(list(R[c].values())):+.2f}', w[1] if K == 'cv' else w]
        if K == 'cv':
            if len({u[0] for u in us}) < 10: log(f'ERROR {c}: {len({u[0] for u in us})} of 10 repeats')
            rm = [np.median([abs(e) for u, e in R[c].items() if u[0] == r]) for r in sorted({u[0] for u in us})]; row += [f'{np.mean(rm):.2f}', f'{min(rm):.2f}', f'{max(rm):.2f}']
        o.write('\t'.join(map(str, row)) + '\n'); log('\t'.join(map(str, row)))
if K == 'cv':   # test females: bp, truth 0
    med = {c: np.median(list(T[c].values())) for c in CL if T[c]}; FB = collections.defaultdict(dict)   # class -> (repeat, fold, female) -> bp
    for d, cs in ((B, None), (F, tuple(FC))) if EF else ((B, ()),):
        for fn in glob.glob(f'{d}/*_fem.tsv'):
            for l in open(fn):
                x = l.rstrip('\n').split('\t'); c = x[2]
                if (c.split('|')[1] if '|' in c else '') != RULE: continue
                c = c.split('|')[0]
                if (cs and c not in cs) or (cs is None and c in FC): continue
                FB[c][(x[0], x[1], x[3])] = float(x[6])
    FB['yq12'] = {u: FB['dyz1'][u] + FB['dyz2'][u] + FB['dyz18s'][u] + FB['yq12rest2'][u] for u in FB['dyz1'] if all(u in FB[c] for c in ('dyz2', 'dyz18s', 'yq12rest2'))}   # Yq12 = DYZ1 + DYZ2 + DYZ18S + rest2, same female, repeat and fold
    log('female Yq12 predictions', len(FB['yq12']), 'from DYZ1', len(FB['dyz1']), 'DYZ2', len(FB['dyz2']), 'DYZ18S', len(FB['dyz18s']), 'rest2', len(FB['yq12rest2']))
    with open(OUT[:-4] + '_females.tsv', 'w') as o:
        o.write('class\tn_predictions\tmedian_bp\tp90_bp\tmax_bp\tmedian_pct_of_median_truth\tmax_pct_of_median_truth\n')
        for c in CL:
            if not FB[c]: continue
            v = np.array(list(FB[c].values())); row = [c, len(v), f'{np.median(v):.0f}', f'{np.percentile(v, 90):.0f}', f'{v.max():.0f}', f'{100 * np.median(v) / med[c]:.4f}', f'{100 * v.max() / med[c]:.4f}']
            o.write('\t'.join(map(str, row)) + '\n'); log('females', '\t'.join(map(str, row)))
log('written', OUT)
