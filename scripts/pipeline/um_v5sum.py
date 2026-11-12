# Summaries of v5 CV, lineage-out, holdout
# usage: UM=markers_v5 (TRUTHX, TRUTHSET as the model) [RULE=r0.5_0.1] python3 um_v5sum.py cv|lolo|holdout <base dir> <floor dir or none> <out.tsv>; also rules, depth and pair modes (below)
# reported classes: floors for dyz2 and yq2k7, Yq12 = DYZ1 + DYZ2 + rest (median Yq12 - DYZ1 - DYZ2 of the training truth males); errors in % of the assembly length
import os, sys, glob, collections, numpy as np
if sys.argv[1] == 'floors':   # floors <out tsv> <g> ...: per class and G+C floor, median |SGDP - HGDP| of the same individuals (model_gcf<g>/same_male_pairs.tsv); no truth; SGDP vs 1kGP NYGC as a control
    P = os.environ['P']; R = collections.defaultdict(lambda: collections.defaultdict(list))
    for g in sys.argv[3:]:
        fn = f'{P}/model_gcf{g}/same_male_pairs.tsv'
        if not os.path.exists(fn): print(f'ERROR floor {g}: no {fn}', flush=True); continue
        for l in list(open(fn))[1:]:
            a, b, pr, c, x, y, d = l.rstrip('\n').split('\t')
            if d not in ('', 'nan') and pr in ('SGDP vs HGDP', 'SGDP vs 1kGP NYGC'): R[(c, pr)][int(g)].append(abs(float(d)))
    with open(sys.argv[2], 'w') as o:
        o.write('class\tpairs\tfloor\tn\tmedian_abs_diff\tmax_abs_diff\n')
        for (c, pr), F in sorted(R.items()):
            for g in sorted(F): v = np.array(F[g]); o.write(f'{c}\t{pr}\t{g}\t{len(v)}\t{np.median(v):.2f}\t{v.max():.2f}\n')
    for c in sorted({c for c, pr in R}):
        F = R[(c, 'SGDP vs HGDP')]
        if not F: continue
        gb = min(F, key=lambda g: (np.median(F[g]), g)); N = R[(c, 'SGDP vs 1kGP NYGC')]
        print(f'{c}: chosen floor {gb} (SGDP vs HGDP median |diff| {np.median(F[gb]):.2f}% over {len(F[gb])} pairs; no floor {np.median(F[0]):.2f}%' + (f'; control SGDP vs NYGC {np.median(N[gb]):.2f}% vs no floor {np.median(N[0]):.2f}%' if 0 in N and gb in N else '') + ')', flush=True)
    sys.exit(0)
if sys.argv[1] == 'rules':   # rules <prefix> <rule> ...: one row per fixed rule from <prefix>_<rule>.tsv, sorted by the mean of class medians
    print('rule\tmean_of_medians\tmax_of_max\t' + '\t'.join(f'{c}_median/p90/max' for c in ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'yq2k7', 'dyz18sub', 'dyz19', 'dyz17']), flush=True); rows = []
    for r in sys.argv[3:]:
        X = {x.split('\t')[0]: x.rstrip('\n').split('\t') for x in list(open(f'{sys.argv[2]}_{r}.tsv'))[1:]}
        rows.append((np.mean([float(v[4]) for v in X.values()]), max(float(v[6]) for v in X.values()), r, '\t'.join(f'{X[c][4]}/{X[c][5]}/{X[c][6]}' if c in X else '.' for c in ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'yq2k7', 'dyz18sub', 'dyz19', 'dyz17'])))
    for m, x, r, t in sorted(rows): print(f'{r}\t{m:.3f}\t{x:.2f}\t{t}', flush=True)
    sys.exit(0)
if sys.argv[1] in ('depth', 'pair'):   # depth table: depth <out.tsv> <label> ...; paired runs: pair <units A> <units B> <name A> <name B> <out.tsv>
    M = f"{os.environ['P']}/{os.environ.get('MDIR', 'model_v5_sum')}"
    if sys.argv[1] == 'depth':
        with open(sys.argv[2], 'w') as o:
            o.write('depth\tclass\tn_predictions\tn_males\tmean_abs\tmedian_abs\tp90_abs\tmax_abs\tbias_median\trepeat_median_mean\tfemale_median_bp\tfemale_max_bp\n')
            for l in sys.argv[3:]:
                R = [x.rstrip('\n').split('\t') for x in open(f'{M}/cv_depth_{l}.tsv')]; F = {x.split('\t')[0]: x.rstrip('\n').split('\t') for x in list(open(f'{M}/cv_depth_{l}_females.tsv'))[1:]}
                for r in R[1:]: o.write('\t'.join([l] + r[:8] + [r[9], F.get(r[0], ['.'] * 5)[2], F.get(r[0], ['.'] * 5)[4]]) + '\n')
        print(open(sys.argv[2]).read(), flush=True); sys.exit(0)
    A, B, NA, NB, OUT = sys.argv[2:7]
    def rd(f):
        D = collections.defaultdict(lambda: collections.defaultdict(list))
        for l in list(open(f))[1:]: c, r, m, e = l.rstrip('\n').split('\t'); D[c][m].append(abs(float(e)))
        return {c: {m: np.mean(v) for m, v in d.items()} for c, d in D.items()}
    X, Y = rd(A), rd(B)
    with open(OUT, 'w') as o:
        o.write(f'class\tcommon_males\tonly_{NA}\tonly_{NB}\tmedian_{NA}\tmedian_{NB}\tp90_{NA}\tp90_{NB}\tmax_{NA}\tmax_{NB}\tmedian_paired_diff\t{NB}_better\t{NA}_better\n')
        for c in X:
            cm = sorted(set(X[c]) & set(Y.get(c, {}))); a = np.array([X[c][m] for m in cm]); b = np.array([Y[c][m] for m in cm]); d = b - a
            row = [c, len(cm), len(set(X[c]) - set(cm)), len(set(Y.get(c, {})) - set(cm)), f'{np.median(a):.2f}', f'{np.median(b):.2f}', f'{np.percentile(a, 90):.2f}', f'{np.percentile(b, 90):.2f}', f'{a.max():.2f}', f'{b.max():.2f}', f'{np.median(d):+.3f}', int((d < 0).sum()), int((d > 0).sum())]
            o.write('\t'.join(map(str, row)) + '\n'); print('\t'.join(map(str, row)), flush=True)
    sys.exit(0)
P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_v5")}'; sys.path.insert(0, f'{P}/scripts'); import truth_cl
K, B, F, OUT = sys.argv[1], f'{O}/{sys.argv[2]}', f'{O}/{sys.argv[3]}', sys.argv[4]
def log(*a): print(*a, flush=True)
T = truth_cl.load(O)
if os.environ.get('TRUTHX'):
    for l in open(os.environ['TRUTHX']): c, s = l.split()[:2]; T[c].pop(s, None)
if os.environ.get('TRUTHSET'):
    for l in open(os.environ['TRUTHSET']):
        c, s, v = l.split()[:3]
        if s in T[c]: T[c][s] = float(v)
G = {}; TRU = {}; RULE = os.environ.get('RULE', '')   # (class, unit) -> training group: (repeat, fold) for cv, lineage for lolo; RULE=r<presence>_<outside>: rows of one fixed rule (um_eval.py ALLRULES)
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
            if len(x) > {'cv': 7, 'lolo': 6}.get(K, 3): TRU[(c, u)] = float(x[-1])   # truth used by um_eval.py (per-fold all-copies truth)
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
