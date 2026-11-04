#!/usr/bin/env python3
"""Reference chrY contigs with the arrays that pass QC, grouped by which arrays pass, rows ordered by similarity.

usage: python3 scripts/figure_array_structure.py <clusters_all.tsv> <out.pdf> <truth_regions.tsv> <dyz19_clean.tsv> <dyz17_clean.tsv> <rule_w4.txt>
rows: assemblies of the reference panels after QC (DYZ3, Yq12, DYZ19, DYZ17), the contig that carries their arrays, forward as assembled;
drawn: the yruler classes only. DYZ3, DYZ17 and DYZ19 = arrays that pass QC; Yq12 = black line over the Yq12 array that passes QC,
with DYZ1, DYZ2 and DYZ18S (DYZ18 + 3.1 kb) clusters >= 10 kb inside it (um_arrstruct.py: HMM hits, clusters within 5 kb);
contigs not named *_chrY (fragments) are shifted so that their first array starts at the median start of that class on the *_chrY contigs;
order in each group: average linkage on log array lengths and log bp and cluster counts of DYZ1, DYZ2, DYZ18S in Yq12 (Euclidean on z-scores)
"""
import sys, csv, lzma, collections as co, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
IN, OUT, TR, D19, D17, RX = sys.argv[1:7]; MINBP = 10000
from palette import CLASS, SOURCE, SOURCE_LABEL, TITLE, NAMES, FS, save, ORDER, CONTIG
COL = {k: CLASS[k] for k in ('DYZ3', 'DYZ17', 'DYZ19', 'DYZ1', 'DYZ2', 'DYZ18S')}
PANEL = ['dyz3', 'dyz17', 'dyz19', 'yq12']; NAME = {'dyz3': 'DYZ3', 'dyz17': 'DYZ17', 'dyz19': 'DYZ19', 'yq12': 'Yq12'}
rd = lambda p: csv.DictReader(lzma.open(p, 'rt') if p.endswith('.xz') else open(p), delimiter='\t')   # inputs may be xz
# arrays passing QC per assembly
A = co.defaultdict(dict)
for r in rd(TR): A[r['asm']][r['class']] = (r['contig'], int(r['start']), int(r['end']))
for c, fn in (('dyz19', D19), ('dyz17', D17)):
    for r in rd(fn):
        if r['clean'] == 'True': A[r['assembly']][c] = (r['contig'], int(r['start']), int(r['end']))
X = [l.split()[:2] for l in open(RX) if l.strip()]   # truth exclusions (class, male), as in the model
for c, m in X:
    for a, d in A.items():
        if a.replace('_R2', '') == m and c in d: del d[c]; print('excluded by rule', c, a, flush=True)
A = {a: d for a, d in A.items() if d}
bad = [a for a, d in A.items() if len({v[0] for v in d.values()}) > 1]
if bad: sys.exit(f'ERROR arrays on more than one contig: {bad}')
CTG = {a: next(iter(d.values()))[0] for a, d in A.items()}; print(len(A), 'assemblies in the panels', {c: sum(c in d for d in A.values()) for c in PANEL}, flush=True)
# clusters of DYZ1, DYZ2, DYZ18S inside Yq12
KQ = {'DYZ1': 'DYZ1', 'DYZ2': 'DYZ2', 'DYZ18': 'DYZ18S', 'Yqhet_3k1bp': 'DYZ18S'}; Q = co.defaultdict(list); L = {}
for r in rd(IN):
    a = r['assembly']
    if a not in A or r['contig'] != CTG[a]: continue
    L[a] = int(r['contig_length']); y = A[a].get('yq12')
    if r['class'] in KQ and y and int(r['annotated_bp']) >= MINBP and int(r['start']) >= y[1] and int(r['end']) <= y[2]: Q[a].append((KQ[r['class']], int(r['start']), int(r['end']), int(r['annotated_bp'])))
miss = [a for a in A if a not in L]
if miss: sys.exit(f'ERROR no clusters for {miss}')
whole = {a for a in A if CTG[a].endswith('_chrY')}
MED = {c: np.median([A[a][c][1] for a in whole if c in A[a]]) for c in PANEL}
SH = {a: 0 if a in whole else MED[min(A[a], key=lambda c: A[a][c][1])] - min(v[1] for v in A[a].values()) for a in A}
print(len(whole), 'assemblies with a *_chrY contig,', len(A) - len(whole), 'fragments (shifted)', flush=True)
for a, d in A.items():
    ks = [c for c in PANEL if c in d]
    if any(d[x][1] > d[y][1] for x, y in zip(ks, ks[1:])): print('WARNING order differs from DYZ3, DYZ17, DYZ19, Yq12:', a, {c: (d[c][1], d[c][2]) for c in ks}, flush=True)
# groups by arrays passing QC
G = co.defaultdict(list)
for a, d in A.items(): G[tuple(c for c in PANEL if c in d)].append(a)
groups = sorted(G, key=lambda k: (-len(k), -len(G[k])))
def feats(a):
    bp, nc = co.Counter(), co.Counter()
    for k, s, e, u in Q[a]: bp[k] += u; nc[k] += 1
    return [np.log10(1 + (A[a][c][2] - A[a][c][1] if c in A[a] else 0)) for c in PANEL] + [np.log10(1 + bp[k]) for k in ('DYZ1', 'DYZ2', 'DYZ18S')] + [np.log10(1 + nc[k]) for k in ('DYZ1', 'DYZ2', 'DYZ18S')]
rows = []   # (assembly or None for a group header, label)
for g in groups:
    m = sorted(G[g])
    if len(m) > 2:
        F = np.array([feats(a) for a in m]); Z = (F - F.mean(0)) / np.where(F.std(0) > 0, F.std(0), 1); m = [m[i] for i in leaves_list(linkage(Z, 'average', metric='euclidean'))]
    rows += [(None, '')] if rows else []; rows.append((None, f'{", ".join(NAME[c] for c in g)} pass QC (n = {len(m)})')); rows += [(a, a.replace('_R2', '')) for a in m]
    print(f'group {g}: {len(m)} assemblies, {sum(a in whole for a in m)} *_chrY contigs', flush=True)
# two columns, 180 mm wide, text 5 to 7 pt
cut = next(i for i, (a, l) in enumerate(rows) if i > 0 and a is None and l == ''); cols = [rows[:cut], rows[cut + 1:]]
n = max(len(c) for c in cols); MM = 1 / 25.4; ROW = 2.0 * MM; fig = plt.figure(figsize=(180 * MM, n * ROW + 20 * MM))
xmax = max(L[a] + SH[a] for a in A) / 1e6; SRC = SOURCE; H = n * ROW / fig.get_figheight(); B = 9 * MM / fig.get_figheight()
for j, rr in enumerate(cols):
    ax = fig.add_axes([0.075 + 0.5 * j, B, 0.415, H])
    for i, (a, lab) in enumerate(rr):
        y = n - 1 - i
        if a is None: ax.text(-0.115 * xmax, y + 0.15, lab, fontsize=FS['tick'], fontweight='bold', va='bottom', clip_on=False); continue
        ax.add_patch(plt.Rectangle((-0.03 * xmax, y + 0.1), 0.018 * xmax, 0.8, fc=SRC['HPRC R2' if a.endswith('_R2') else 'Hallast 2026'], lw=0, clip_on=False))
        s0 = SH[a] / 1e6; ax.add_patch(plt.Rectangle((s0, y + 0.1), L[a] / 1e6, 0.8, fc=CONTIG['face'], ec=CONTIG['edge'], lw=0.25))
        for c in ('dyz3', 'dyz17', 'dyz19'):
            if c in A[a]: _, s, e = A[a][c]; ax.add_patch(plt.Rectangle((s0 + s / 1e6, y + 0.1), max((e - s) / 1e6, 0.1), 0.8, fc=COL[NAME[c]], lw=0))
        for k, s, e, _ in sorted(Q[a], key=lambda x: x[0] == 'DYZ18S'): ax.add_patch(plt.Rectangle((s0 + s / 1e6, y + 0.1), max((e - s) / 1e6, 0.05), 0.8, fc=COL[k], lw=0))
        if 'yq12' in A[a]: _, s, e = A[a]['yq12']; ax.plot([s0 + s / 1e6, s0 + e / 1e6], [y + 0.97, y + 0.97], color='black', lw=0.4, solid_capstyle='butt')
    ax.set_xlim(0, xmax * 1.01); ax.set_ylim(n - len(rr), n); ax.set_yticks([n - 1 - i + 0.5 for i, (a, _) in enumerate(rr) if a])
    ax.set_yticklabels([l for a, l in rr if a]); ax.tick_params(axis='y', length=0, pad=0.035 * 72 * 0.415 * fig.get_figwidth())
    ax.set_xlabel('Position along the contig (Mb)', labelpad=1)
    for sp in ('top', 'right', 'left'): ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.4)
y0 = B + H + 1.5 * MM / fig.get_figheight(); kw = dict(loc='lower center', title_fontproperties=TITLE(FS['legend']))
K = [k for k in ORDER if k in COL or k == 'Yq12']; lc = fig.legend([plt.Line2D([0], [0], color='black', lw=0.8) if k == 'Yq12' else plt.Rectangle((0, 0), 1, 1, fc=COL[k]) for k in K], K, ncol=7, title=NAMES['class'], bbox_to_anchor=(0.30, y0), **kw)
lp = fig.legend([plt.Rectangle((0, 0), 1, 1, fc=c) for c in SRC.values()], [SOURCE_LABEL[k] for k in SRC], ncol=2, title=NAMES['source'], bbox_to_anchor=(0.78, y0), **kw)
save(fig, OUT); print(sum(1 for a, _ in rows if a), 'assemblies,', len(groups), 'groups', flush=True)
for k in ('DYZ1', 'DYZ2', 'DYZ18S'):
    c = [sum(x[0] == k for x in Q[a]) for a in A if 'yq12' in A[a]]; print(k, 'clusters >= 10 kb in Yq12 per assembly: median', np.median(c), 'min', min(c), 'max', max(c))
