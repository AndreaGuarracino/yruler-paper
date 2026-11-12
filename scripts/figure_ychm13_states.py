# Supplementary figure: chrY dose of unusual males
# usage: python3 scripts/figure_ychm13_states.py <dir with ychm13_flagged_100kb.tsv, .classes.tsv> <out.pdf> [unusual|pairs] [kmer <X/PAR summary.tsv> [Y-specific profile.tsv [CRAM-subset profile.tsv:sample,...]]]; 16:9, rows clustered within groups
# kmer: alignment-free only (marker k-mer dose on Y, no read depth; X dose from um_sexprof.py; PAR and XTR windows: Y = 2 (X + Y) - X from <summary>.chrY_XY.tsv)
from palette import NAMES, TITLE, MM, save
import palette as PL
import sys, csv, textwrap, warnings, collections as co, numpy as np, matplotlib
warnings.filterwarnings('ignore', 'Mean of empty slice')
matplotlib.use('Agg'); import matplotlib.pyplot as plt; plt.rcParams['hatch.color'] = PL.STAIN['gpos50']; plt.rcParams['hatch.linewidth'] = 0.4
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering
D, OUT = sys.argv[1:3]; MODE = sys.argv[3] if len(sys.argv) > 3 else 'unusual'; KM = sys.argv[5] if len(sys.argv) > 5 and sys.argv[4] == 'kmer' else None; MB = 1e6; INK = PL.INK   # MODE pairs: father-son pairs only
BANDS = [('p11.32', 0, 127522, 'gneg'), ('p11.31', 127522, 454425, 'gpos50'), ('p11.2', 454425, 10568667, 'gneg'), ('p11.1', 10568667, 10703039, 'acen'), ('cen', 10703039, 10726286, 'acen'),   # KaryoScope-ISCN chm13v2.0 cytobands.bed
         ('q11.1', 10726286, 10883195, 'acen'), ('q11.21', 10883195, 13307633, 'gneg'), ('q11.221', 13307633, 18006534, 'gpos50'), ('q11.222', 18006534, 20506564, 'gneg'),
         ('q11.223', 20506564, 25062557, 'gpos50'), ('q11.23', 25062557, 27449937, 'gneg'), ('q12', 27449937, 62460029, 'gvar')]
STAIN = PL.STAIN
GROUPS = [('Yp \u2260 Yq (gain and loss)', 'partial'), ('Yq loss', 'yq_loss'), ('Yp depletion', 'yp_deletion'), ('Low Y dosage', 'mosaic_loss'), ('Intermediate Y dosage', 'mosaic_gain'), ('Y gain', 'Y_gain'), ('PAR2 gain', 'par2_gain'), ('No Y call', 'no_Y')]
EXTRA = [('X gain', ['HG01683']), ('Unusual array estimates', ['HG01527', 'HG01529', 'NA18910', 'NA18911', 'HG04006', 'HG02804', 'HG02684', 'HG02686', 'NA18952']), ('Recorded females with Y signal', ['HGDP01239', 'HGDP01273', 'HG02772'])]   # females: user 2026-11-11
NOTE = {'HGDP01239': 'recorded female; XY', 'HGDP01273': 'recorded female; XY', 'HG02772': 'mother of HG02773; Yq12 and PAR2 without single-copy Y', 'NA19226': 'Y lost; offspring XY', 'HG02966': 'offspring; parent HG02965 without change', 'HG01529': 'DYZ18S', 'HG01527': 'DYZ18S', 'NA18911': 'DYZ19', 'NA18910': 'DYZ19', 'HG04006': 'DYZ3', 'HG02804': 'DYZ2', 'NA18952': 'DYZ3', 'HG02684': 'H3b2a', 'HG02686': 'H3b2a', 'NA11930': 'PAR2 1.43-1.51 in 3 libraries', 'HG02895': 'PAR2 1.45, 1.47 in 2 ONT runs', 'HG03685': 'PAR2 1.58 in ONT; Yq12-PAR2 on 15p', 'HG02773': 'PAR2 1.50 in Polaris'}
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
if not KM: sys.exit('ERROR alignment-free only: give kmer <X/PAR summary.tsv>')
YP = sys.argv[6] if KM and len(sys.argv) > 6 else None   # um_yprof2.py windows: Y-specific 31-mers
SUB = sys.argv[7].split(':') if YP and len(sys.argv) > 7 else None; SUBS = set(SUB[1].split(',')) if SUB else set()   # samples whose Y comes from the CRAM subset
K = [r for r in rd(YP) if r['level'] == 'window' and r['role'] != 'reference' and r['sample'] not in SUBS] + ([r for r in rd(SUB[0]) if r['chrom'] == 'chrY' and r['sample'] in SUBS] if SUB else []) if YP else [r for r in rd(f'{D}/ychm13_flagged_100kb.tsv') if r['level'] == 'window']
W = sorted({int(r['start']) for r in K}); step = W[1] - W[0]; W = list(range(0, 62460029, step)) if KM else W; nw = len(W); X = {}   # kmer: whole chrY, PAR2 included
for r in K:   # k-mer dose where >= 50 marker positions
    X.setdefault(r['sample'], np.full(nw, np.nan))[int(r['start']) // step] = float(r['ratio']) if int(r['key_positions']) >= 50 else np.nan
s0 = K[0]['sample']; kw = {int(r['start']) // step for r in K if r['sample'] == s0 and int(r['key_positions']) >= 50}
print(f'windows: {nw}; with k-mer dose {len(kw)}, no data {nw - len(kw)}')
E = {r['sample']: r for r in rd('results/estimates_1kgp.tsv')}; FEM = set(EXTRA[-1][1]); CT = {r['sample']: r for f in ('cohort_1kgp', 'cohort_fem', 'cohort_hgdp_sgdp') for r in rd(f'results/model/{f}.tsv')}
X = {s: v for s, v in X.items() if s in FEM or (s in E and E[s]['sex_ped'] == '1')}   # males, and the recorded females with a Y signal
print(f'recorded females with Y signal in the profile: {sorted(FEM & set(X))}', flush=True)
XP = {r['sample']: r for r in rd(KM)} if KM else {}; fx = lambda s, c: float(XP[s][c]) if s in XP else np.nan   # alignment-free X and PAR doses
import os; FX = {r['sample']: r for r in rd('results/checks/ed3_female_doses.tsv')} if os.path.exists('results/checks/ed3_female_doses.tsv') else {}   # females outside the cohort tables: doses from their recount (um_apply59.py)
CT.update({s: r for s, r in FX.items() if s not in CT}); PR = {r['sample']: r for r in K}; yfem = lambda s: float(CT[s]['y_dose']) if s in CT else float(PR[s]['y_dose_saved'])   # females outside the cohort tables: dose of the profile reads
info = {s: (CT[s]['y_state'] if s in CT else 'no_Y' if yfem(s) <= 0.025 else 'normal', yfem(s), fx(s, 'x_kmer') if s in XP else float(CT[s]['x_dose']) if s in CT else np.nan) if s in FEM else (E[s]['y_state'], float(E[s]['y_dose']), fx(s, 'x_kmer') if KM else float(E[s]['x_dose'])) for s in X}
XYS = [(int(x[1]), int(x[2])) for x in (l.split('\t') for l in open('data/chm13_annot/chm13v2.0_chrXY_sequence_class_v1.bed')) if x[0] == 'chrY' and x[3] in ('PAR1', 'PAR2', 'XTR1', 'XTR2')]   # PAR and XTR windows from the X + Y rulers
if KM:   # PAR and XTR windows: Y dose from the 31-mers shared by X and Y, Y = 2 (X + Y) - X
    nxy = 0
    for r in rd(KM.replace('.summary.tsv', '.chrY_XY.tsv')):
        j, a = int(r['start']) // step, int(r['start'])
        if r['sample'] in X and j < nw and any(p <= a + step / 2 < e for p, e in XYS) and int(r['key_positions']) >= 50 and r['ratio'] != 'nan': X[r['sample']][j] = 2 * float(r['ratio']) - info[r['sample']][2]; nxy += r['sample'] == s0
    print(f'PAR and XTR windows from X + Y keys in {s0}: {nxy}; windows still without data: {int(np.isnan(X[s0]).sum())} of {nw}')
def order(S):   # rows of a group by profile similarity
    if len(S) < 3: return S
    A = np.nan_to_num(np.clip(np.array([X[s] for s in S]), 0, 2), nan=1.0); return [S[i] for i in leaves_list(optimal_leaf_ordering(linkage(A, 'average'), A))]
ORD = __import__('os').environ.get('ORDER', 'window')   # window: 100-kb profiles; band: mean per cytoband
BI = [ix for ix in ([j for j in range(nw) if a <= j * step < b] for _, a, b, _ in BANDS) if ix]
feat = (lambda v: np.array([np.nanmean(v[ix]) if np.isfinite(v[ix]).any() else np.nan for ix in BI])) if ORD == 'band' else (lambda v: v)
print(f'row order within groups: {ORD} ({len(BI) if ORD == "band" else nw} features)', flush=True)
FS = [(r['parent'], r['child']) for r in rd('results/estimates_1kgp.relationships.tsv') if r['relationship'] == 'parent_offspring' and r['parent'] in X and r['child'] in X]   # father-son pairs shown
grp = {s: st for s in X for lab, st in GROUPS if info[s][0] == st}; grp.update({s: lab for lab, S in EXTRA for s in S if s in X})
REL = co.defaultdict(list)
for p, c in FS: REL[p].append((c, 'parent of')); REL[c].append((p, 'offspring of'))
rows, gl, relrow = [], [], set()
for lab, key in [(l, st) for l, st in GROUPS] + [(l, l) for l, S in EXTRA]:
    M = sorted(s for s in X if grp.get(s) == key); units, seen = [], set()
    for s in M:   # members linked by father-son pairs stay together; relatives outside every group join their unit
        if s in seen: continue
        u, todo = [], [s]
        while todo:
            x = todo.pop()
            if x in seen: continue
            seen.add(x); u.append(x); todo += [y for y, _ in REL[x] if y not in seen and grp.get(y) == key]
        units.append(sorted(u, key=lambda x: sum(1 for y, r in REL[x] if r == 'offspring of' and y in u)))   # fathers above sons
    if len(units) > 2:
        A = np.nan_to_num(np.clip(np.array([feat(np.nanmean([X[x] for x in u if x in M], 0)) for u in units]), 0, 2), nan=1.0); units = [units[i] for i in leaves_list(optimal_leaf_ordering(linkage(A, 'average'), A))]
    S = [x for u in units for x in u]; relrow |= {x for x in S if x not in grp}
    if S: gl.append((lab, len(rows), len(rows) + len(S))); rows.extend(S)
BASE = dict(NOTE)
for x in rows:   # father/son notes between shown males
    for y, r in REL[x]:
        if y in rows: NOTE[x] = (NOTE[x] + '; ' if x in NOTE else '') + f'{r} {y}'
if MODE == 'pairs':   # every father-son pair with an unusual member, father above son
    lab = {v: k for k, v in GROUPS}; lab.update({k: k for k, _ in EXTRA}); lab['Unusual array estimates'] = 'Unusual array'; gname = lambda x: lab[grp[x]] if x in grp else None
    def plab(p, c):
        a, b = gname(p), gname(c)
        return f'{a}: both' if a and a == b else '; '.join(f'{g}: {r}' for r, g in (('parent', a), ('offspring', b)) if g)
    P = sorted((pc for pc in FS if pc[0] in grp or pc[1] in grp), key=lambda pc: plab(*pc)); rows, gl = [], []
    for p, c in P: gl.append((plab(p, c), len(rows), len(rows) + 2)); rows += [p, c]
    BASE.update({'NA19226': 'Y lost in the sample'}); BASE.pop('HG02966', None)
    NOTE = {x: ('parent' if i % 2 == 0 else 'offspring') + (f'; {BASE[x]}' if x in BASE else '') for i, x in enumerate(rows)}
print(f'{len(rows)} males in {len(gl)} groups; missing {sorted(set(X) - set(rows))}')
print('row notes not drawn (for the caption):', {x: NOTE[x] for x in rows if x in NOTE}, flush=True)
C = co.defaultdict(lambda: np.zeros(nw))
for r in rd(f'{D}/ychm13_flagged_100kb.classes.tsv'): C[r['class']][int(r['start']) // step] = int(r['positions'])
BOUND = [0, 0.25, 0.75, 0.9, 1.1, 1.25, 1.75, 2.5]   # 0, about 0.5, slight loss, 1, slight gain, about 1.5, 2
cm = ListedColormap(PL.HEAT); cm.set_bad('white'); nrm = BoundaryNorm(BOUND, cm.N)
n = len(rows); ROW, HEAD, BOT, TOP = 2.5, 25, 8, 1   # mm: row pitch, tracks, bottom, top margins
FLAB = {'HGDP01239': 'XY profile', 'HGDP01273': 'XY profile', 'HG02772': 'Yq12 only'}   # notes per female (user 2026-11-11)
runs = lambda v: [(v[i], i, j) for i, j in zip([0] + [k for k in range(1, len(v)) if v[k] != v[k - 1]], [k for k in range(1, len(v)) if v[k] != v[k - 1]] + [len(v)])]   # consecutive equal labels
FG = [i for i, (lab, a, b) in enumerate(gl) if lab == EXTRA[-1][0]]; nf = gl[FG[0]][2] - gl[FG[0]][1] if FG else 0; GAP = 7 if nf else 0   # females: own block below the males
FW, FH = 180, TOP + HEAD + 6 + n * ROW + GAP + BOT; fig = plt.figure(figsize=(FW * MM, FH * MM))
ax_mm = lambda x0, x1, y0, y1: fig.add_axes([x0 / FW, y0 / FH, (x1 - x0) / FW, (y1 - y0) / FH])
hb = BOT + n * ROW + GAP   # heatmap top
ax0, cax = ax_mm(14, 126, hb + 6, hb + 6 + HEAD), ax_mm(168.5, 170.5, BOT, hb)
BL = [(rows[:n - nf], [g for g in gl if g[0] != EXTRA[-1][0]], BOT + nf * ROW + GAP, hb)] + ([(rows[n - nf:], [(lab, i, j) for lab, i, j in runs([FLAB[s] for s in rows[n - nf:]])], BOT, BOT + nf * ROW)] if nf else [])
xmax = (W[-1] + step) / MB
A = 'data/chm13_annot'   # T2T CHM13v2 annotations (github.com/marbl/CHM13)
SC = [(int(x[1]), int(x[2]), x[3], '#%02x%02x%02x' % tuple(int(v) for v in x[8].split(','))) for x in (l.split('\t') for l in open(f'{A}/chm13v2.0_chrXY_sequence_class_v1.bed')) if x[0] == 'chrY']
PAL = co.defaultdict(list)
for l in open(f'{A}/chm13v2.0Y_inverted_repeats_v1.bed'):
    x = l.split('\t')
    if x[0] == 'chrY' and x[3].startswith('P'): PAL[x[3].strip()].append((int(x[1]), int(x[2])))
NMC = {'PAR1': 'PAR', 'PAR2': 'PAR', 'XTR1': 'XTR', 'XTR2': 'XTR', 'X-DEG': 'X-degenerate', 'AMPL': 'Ampliconic', 'SAT': 'Satellite', 'DYZ17': 'Satellite', 'DYZ19': 'Satellite', 'CEN': 'Centromere', 'HET': 'Heterochromatin', 'OTHER': 'Other'}
for nm, a, b, st in BANDS:   # tracks in mm: cytobands 10.8-13.8, names above on two rows
    ax0.add_patch(plt.Rectangle((a / MB, 10.8), (b - a) / MB, 3, fc=STAIN[st], ec=PL.STAIN['acen'], lw=0.5, hatch='////' if st == 'gvar' else None))
for i, (nm, a, b, st) in enumerate(B for B in BANDS if B[2] - B[1] > 1e6):
    y = 14.4 + 2.8 * (i % 2); ax0.plot([(a + b) / 2 / MB] * 2, [13.8, y - 0.2], color=INK, lw=0.4); ax0.text((a + b) / 2 / MB, y, nm, ha='center', va='bottom', fontsize=PL.FS["label"], color=INK)
for lab, a, b in (('Yp', 0, 10.71), ('Yq', 10.71, 62.46)):   # arms
    ax0.plot([a + 0.15, b - 0.15], [21.3, 21.3], color=INK, lw=0.8); ax0.text((a + b) / 2, 21.6, lab, ha='center', va='bottom', fontsize=PL.FS["label"], color=INK, fontweight='bold')
CC = {}   # sequence classes 7.4-9.8, legend at the top right
for a, b, c, col in SC:
    ax0.add_patch(plt.Rectangle((a / MB, 7.4), (b - a) / MB, 2.4, fc=col, ec='white', lw=0.3)); CC.setdefault(NMC[c], col)
lax = ax_mm(128, 180, hb + 6, hb + 6 + HEAD); lax.axis('off')   # legend in its own header cell
lax.legend([plt.Rectangle((0, 0), 1, 1, fc=c) for c in CC.values()], list(CC), loc='center', ncol=2, title=NAMES['annotation'], title_fontproperties=TITLE(PL.FS['legend']), labelspacing=0.35)
for i, (p, arms) in enumerate(sorted(PAL.items(), key=lambda kv: np.mean(kv[1]))):   # palindromes 5.2-6.4, labels on two rows
    for a, b in arms: ax0.add_patch(plt.Rectangle((a / MB, 5.2), (b - a) / MB, 1.2, fc=INK, lw=0))
    ax0.text(np.mean([v for ab in arms for v in ab]) / MB, 4.9 - 2.4 * (i % 2), p, ha='center', va='top', fontsize=PL.FS["point"], color=INK)
for y, lab in ((12.3, 'Cytoband'), (8.6, NAMES['annotation']), (5.8, 'Palindromes')): ax0.text(-0.5, y, lab, ha='right', va='center', fontsize=PL.FS["tick"], color=INK, linespacing=1)
ax0.set_xlim(0, xmax); ax0.set_ylim(0, HEAD); ax0.axis('off')
nd_all = 0
for k, (R, G, y0, y1) in enumerate(BL):   # one block per sex: heatmap, dose bars, group labels
    m = len(R); ax2, ay, ax, axg = ax_mm(14, 126, y0, y1), ax_mm(128, 137, y0, y1), ax_mm(139, 148, y0, y1), ax_mm(149.5, 167, y0, y1)
    M = np.array([X[s] for s in R]); im = ax2.imshow(np.clip(M, 0, 2.49), aspect='auto', cmap=cm, norm=nrm, interpolation='nearest', extent=(0, xmax, m, 0))
    ax2.hlines(np.arange(1, m), 0, xmax, color='white', lw=0.5)   # thin gap between rows
    nd = np.isnan(M).astype(float); nd_all += int(nd.sum())
    if nd.any(): ax2.contourf(np.array(W) / MB + step / 2 / MB, np.arange(m) + 0.5, nd, levels=[0.5, 1.5], colors='none', hatches=['////'])   # no data: hatched white
    ax2.set_yticks(np.arange(m) + 0.5); ax2.set_yticklabels(R)
    for t, s in zip(ax2.get_yticklabels(), R): t.set_color(PL.GREY if s in relrow else INK)   # relatives outside the groups in grey
    ax2.set_xlim(0, xmax); last = k == len(BL) - 1
    if last: ax2.set_xlabel('CHM13v2 chrY position (Mb)')
    else: ax2.set_xticklabels([])
    if nf: ax2.set_title('Recorded females' if k else 'Recorded males', loc='left', fontsize=PL.FS['label'], fontweight='bold', pad=2)
    for a_, lab, col in ((ay, 'Y dose', 1), (ax, 'X dose', 2)):   # dose bars with numbers
        v = [info[s][col] for s in R]; a_.barh(np.arange(m) + 0.5, np.nan_to_num(v), height=0.62, color=PL.GREY, lw=0)
        for i, x in enumerate(v):
            if x == x: a_.text(min(x, 2.25) + 0.05, i + 0.5, f'{x:.2f}', va='center', ha='left', fontsize=PL.FS["point"], color=INK)
        for g in (1, 2): a_.axvline(g, color=PL.GREY, lw=0.5, ls=':')
        a_.set_xlim(0, 2.9); a_.set_ylim(m, 0); a_.set_yticks([]); a_.set_xticks([0, 1, 2]); a_.set_xlabel(lab if last else '')
        if not last: a_.set_xticklabels([])
        for sp_ in ('top', 'right'): a_.spines[sp_].set_visible(False)
    axg.set_ylim(m, 0); axg.set_xlim(0, 1); axg.set_autoscale_on(False); axg.axis('off')
    for lab, a, b in G:
        for a_ in (ax2, ay, ax): a_.axhline(a, color=INK, lw=0.6)
        axg.plot([0.08, 0.08], [a + 0.15, b - 0.15], color=INK, lw=0.6); axg.text(0.2, (a + b) / 2, textwrap.fill(lab, 12), va='center', ha='left', fontsize=PL.FS["tick"], color=INK)
print(f'no-data cells in the heatmap: {nd_all} of {n * nw}', flush=True); nd = np.array([nd_all])
cb = fig.colorbar(im, cax=cax, boundaries=BOUND, ticks=[(a + b) / 2 for a, b in zip(BOUND[:-1], BOUND[1:])]); cb.ax.set_yticklabels(['0', '0.5', 'loss', '1', 'gain', '1.5', '2'])
cb.set_label('Dose relative to 100 XY males without flagged Y profiles' + ('; hatched: no data' if nd.any() else ''))
save(fig, OUT)
