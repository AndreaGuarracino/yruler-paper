# BED vs whole: length ratio spread, reads, time
import sys, glob, os, re, statistics as st
P = '/path/to/project'; MODES = sys.argv[1:] or ['bedv2', 'bedv2w']; NT = 30   # first 30 males chose v2w windows
def summ(d):
    f = f'{d}/summary.tsv'
    if not os.path.exists(f): return None
    r = [l.rstrip('\n').split('\t') for l in open(f) if not l.startswith(('#', 'sample'))]
    v = {x[2]: float(x[7]) for x in r}; y = float([x for x in r if x[2] == 'ysinglecopy'][0][4]); return {c: v[c] / y for c in v if c not in ('ysinglecopy', 'chr2depth')}
def run(d):
    t = open(f'{d}/run.log').read(); n = re.search(r'records read \(primary\) (\d+)', t); s = re.findall(r'^([\d.]+) s (\d+) kB', t, re.M)
    return (int(n.group(1)) if n else None, float(s[-1][0]) if s else None)
males = sorted(os.path.basename(os.path.dirname(d)) for d in glob.glob(f'{P}/yruler_dragen/*/whole') if summ(d))
print(f'males with whole-file summary {len(males)}; v2w windows chosen on the first {NT}')
for m in MODES:
    ok = [s for s in males if summ(f'{P}/yruler_dragen/{s}/{m}')]
    for lab, grp in (('all', ok), (f'first {NT}', [s for s in ok if s in males[:NT]]), ('other (holdout)', [s for s in ok if s not in males[:NT]])):
        if not grp: continue
        print(f'== {m} {lab}: {len(grp)} males')
        if lab == 'all':
            rw = [run(f'{P}/yruler_dragen/{s}/{m}') for s in grp]; print(f'  records mean {st.mean(x[0] for x in rw if x[0]):.0f}, time mean {st.mean(x[1] for x in rw if x[1]):.2f} s')
        for c in summ(f'{P}/yruler_dragen/{grp[0]}/whole'):
            r = [summ(f'{P}/yruler_dragen/{s}/{m}')[c] / summ(f'{P}/yruler_dragen/{s}/whole')[c] for s in grp if summ(f'{P}/yruler_dragen/{s}/whole')[c] > 0]
            sp = 100 * (max(r) - min(r)); print(f'  {c:6s} n {len(r)} ratio median {st.median(r):.4f} min {min(r):.4f} max {max(r):.4f} spread {sp:.2f} points {"PASS" if sp < 0.5 else "FAIL"}')
