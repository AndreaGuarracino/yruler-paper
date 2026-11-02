# Supplementary figure: chrY dose of unusual males
# usage: python3 scripts/figure_ychm13_states.py <dir with ychm13_flagged_100kb.tsv, .classes.tsv, ylift_flagged_100kb.tsv> <out.pdf> [unusual|pairs] [kmer <X/PAR summary.tsv> [Y-specific profile.tsv [CRAM-subset profile.tsv:sample,...]]]; 16:9, rows clustered within groups
# kmer: alignment-free only (marker k-mer dose on Y, no read depth; X dose from um_sexprof.py; PAR and XTR windows: Y = 2 (X + Y) - X from <summary>.chrY_XY.tsv)
import sys, csv, textwrap, warnings, collections as co, numpy as np, matplotlib
warnings.filterwarnings('ignore', 'Mean of empty slice')
matplotlib.use('Agg'); import matplotlib.pyplot as plt; plt.rcParams['hatch.color'] = '#a9a8a4'; plt.rcParams['hatch.linewidth'] = 0.4
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering
D, OUT = sys.argv[1:3]; MODE = sys.argv[3] if len(sys.argv) > 3 else 'unusual'; KM = sys.argv[5] if len(sys.argv) > 5 and sys.argv[4] == 'kmer' else None; MB = 1e6; INK = '#3d3c39'   # MODE pairs: father-son pairs only
BANDS = [('p11.32', 0, 127522, 'gneg'), ('p11.31', 127522, 454425, 'gpos50'), ('p11.2', 454425, 10568667, 'gneg'), ('p11.1', 10568667, 10703039, 'acen'), ('cen', 10703039, 10726286, 'acen'),   # KaryoScope-ISCN chm13v2.0 cytobands.bed
         ('q11.1', 10726286, 10883195, 'acen'), ('q11.21', 10883195, 13307633, 'gneg'), ('q11.221', 13307633, 18006534, 'gpos50'), ('q11.222', 18006534, 20506564, 'gneg'),
         ('q11.223', 20506564, 25062557, 'gpos50'), ('q11.23', 25062557, 27449937, 'gneg'), ('q12', 27449937, 62460029, 'gvar')]
STAIN = {'gneg': '#fcfcfb', 'gpos50': '#a9a8a4', 'acen': '#5d5c58', 'gvar': '#d6d5d1'}
GROUPS = [('Yp \u2260 Yq (gain and loss)', 'partial'), ('Yq loss', 'yq_loss'), ('Yp deletion', 'yp_deletion'), ('Mosaic loss', 'mosaic_loss'), ('Mosaic gain', 'mosaic_gain'), ('Y gain (XYY)', 'Y_gain'), ('PAR2 gain', 'par2_gain'), ('No Y call', 'no_Y')]
EXTRA = [('X gain (XXY)', ['HG01683']), ('Unusual array estimates', ['HG01527', 'HG01529', 'NA18910', 'NA18911', 'HG04006', 'HG02804', 'HG02684', 'HG02686', 'NA18952'])]
NOTE = {'HG03006': 'CRAM', 'NA19226': 'Y lost; offspring XY', 'HG02966': 'parent without change', 'HG01529': 'DYZ18S', 'HG01527': 'DYZ18S', 'NA18911': 'DYZ19', 'NA18910': 'DYZ19', 'HG04006': 'DYZ3', 'HG02804': 'DYZ2', 'NA18952': 'DYZ3', 'HG02684': 'H3b2a', 'HG02686': 'H3b2a', 'NA11930': 'PAR2 1.43-1.51 in 3 libraries', 'HG02895': 'PAR2 1.47 in ONT', 'HG03685': 'PAR2 1.58 in ONT'}
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
YP = sys.argv[6] if KM and len(sys.argv) > 6 else None   # um_yprof2.py windows: Y-specific 31-mers
SUB = sys.argv[7].split(':') if YP and len(sys.argv) > 7 else None; SUBS = set(SUB[1].split(',')) if SUB else set()   # samples whose Y comes from the CRAM subset
K = [r for r in rd(YP) if r['level'] == 'window' and r['role'] != 'reference' and r['sample'] not in SUBS] + [r for r in rd(SUB[0]) if r['chrom'] == 'chrY' and r['sample'] in SUBS] if YP else [r for r in rd(f'{D}/ychm13_flagged_100kb.tsv') if r['level'] == 'window']; L = [] if KM else rd(f'{D}/ylift_flagged_100kb.tsv')
W = sorted({int(r['start']) for r in K}); step = W[1] - W[0]; W = list(range(0, 62460029, step)) if KM else W; nw = len(W); X = {}   # kmer: whole chrY, PAR2 included
for r in K:   # k-mer dose where >= 50 marker positions
    X.setdefault(r['sample'], np.full(nw, np.nan))[int(r['start']) // step] = float(r['ratio']) if int(r['key_positions']) >= 50 else np.nan
for r in L:   # read depth where normal males have >= 50 unique reads per million chr2 reads
    j = int(r['start']) // step
    if r['sample'] in X and j < nw and float(r['ref_reads_per_million_chr2']) >= 50 and r['ratio'] != 'nan': X[r['sample']][j] = float(r['ratio'])
s0 = K[0]['sample']; kw = {int(r['start']) // step for r in K if r['sample'] == s0 and int(r['key_positions']) >= 50}
lw = {int(r['start']) // step for r in L if r['sample'] == s0 and float(r['ref_reads_per_million_chr2']) >= 50 and int(r['start']) // step < nw}
print(f'windows: {nw}; read depth {len(lw)}, k-mers only {len(kw - lw)}, no data {nw - len(lw | kw)}')
E = {r['sample']: r for r in rd('results/estimates_1kgp.tsv')}; X = {s: v for s, v in X.items() if E[s]['sex_ped'] == '1'}   # males only
XP = {r['sample']: r for r in rd(KM)} if KM else {}; fx = lambda s, c: float(XP[s][c]) if s in XP else np.nan   # alignment-free X and PAR doses
info = {s: (E[s]['y_state'], float(E[s]['y_dose']), fx(s, 'x_kmer') if KM else float(E[s]['x_dose'])) for s in X}
XYS = [(int(x[1]), int(x[2])) for x in (l.split('\t') for l in open('data/chm13_annot/chm13v2.0_chrXY_sequence_class_v1.bed')) if x[0] == 'chrY' and x[3] in (('PAR1', 'PAR2') if YP else ('PAR1', 'PAR2', 'XTR1', 'XTR2'))]   # with Y-specific keys only the PARs
if KM:   # PAR and XTR windows: Y dose from the 31-mers shared by X and Y, Y = 2 (X + Y) - X
    nxy = 0
    for r in rd(KM.replace('.summary.tsv', '.chrY_XY.tsv')):
        j, a = int(r['start']) // step, int(r['start'])
        if r['sample'] in X and j < nw and any(p <= a + step / 2 < e for p, e in XYS) and int(r['key_positions']) >= 50 and r['ratio'] != 'nan': X[r['sample']][j] = 2 * float(r['ratio']) - info[r['sample']][2]; nxy += r['sample'] == s0
    print(f'PAR and XTR windows from X + Y keys in {s0}: {nxy}; windows still without data: {int(np.isnan(X[s0]).sum())} of {nw}')
def order(S):   # rows of a group by profile similarity
    if len(S) < 3: return S
    A = np.nan_to_num(np.clip(np.array([X[s] for s in S]), 0, 2), nan=1.0); return [S[i] for i in leaves_list(optimal_leaf_ordering(linkage(A, 'average'), A))]
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
        A = np.nan_to_num(np.clip(np.array([np.nanmean([X[x] for x in u if x in M], 0) for u in units]), 0, 2), nan=1.0); units = [units[i] for i in leaves_list(optimal_leaf_ordering(linkage(A, 'average'), A))]
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
C = co.defaultdict(lambda: np.zeros(nw))
for r in rd(f'{D}/ychm13_flagged_100kb.classes.tsv'): C[r['class']][int(r['start']) // step] = int(r['positions'])
BOUND = [0, 0.25, 0.75, 0.9, 1.1, 1.25, 1.75, 2.5]   # 0, about 0.5, slight loss, 1, slight gain, about 1.5, 2
cm = ListedColormap(['#123c73', '#2a78d6', '#9cc3ee', '#d4d3ce', '#f2a6a5', '#e34948', '#8f1f1f']); cm.set_bad('#ffffff'); nrm = BoundaryNorm(BOUND, cm.N)
fig = plt.figure(figsize=(13.33, 7.5)); n = len(rows)
gs = fig.add_gridspec(2, 5, height_ratios=[0.2, 0.8], width_ratios=[1, 0.06, 0.06, 0.1, 0.012], hspace=0.04, wspace=0.03, left=0.185, right=0.935, top=0.97, bottom=0.215)
ax0, ax2, ay, ax, axg, cax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2]), fig.add_subplot(gs[1, 3]), fig.add_subplot(gs[1, 4])
xmax = (W[-1] + step) / MB
A = 'data/chm13_annot'   # T2T CHM13v2 annotations (github.com/marbl/CHM13)
SC = [(int(x[1]), int(x[2]), x[3], '#%02x%02x%02x' % tuple(int(v) for v in x[8].split(','))) for x in (l.split('\t') for l in open(f'{A}/chm13v2.0_chrXY_sequence_class_v1.bed')) if x[0] == 'chrY']
PAL = co.defaultdict(list)
for l in open(f'{A}/chm13v2.0Y_inverted_repeats_v1.bed'):
    x = l.split('\t')
    if x[0] == 'chrY' and x[3].startswith('P'): PAL[x[3].strip()].append((int(x[1]), int(x[2])))
NMC = {'PAR1': 'PAR', 'PAR2': 'PAR', 'XTR1': 'XTR', 'XTR2': 'XTR', 'X-DEG': 'X-degenerate', 'AMPL': 'Ampliconic', 'SAT': 'Satellite', 'DYZ17': 'Satellite', 'DYZ19': 'Satellite', 'CEN': 'Centromere', 'HET': 'Heterochromatin', 'OTHER': 'Other'}
for nm, a, b, st in BANDS:   # cytobands (y 1.9-2.6), band names above
    ax0.add_patch(plt.Rectangle((a / MB, 1.9), (b - a) / MB, 0.7, fc=STAIN[st], ec='#5d5c58', lw=0.5, hatch='////' if st == 'gvar' else None))
for nm, a, b, st in (B for B in BANDS if B[2] - B[1] > 1e6): ax0.text((a + b) / 2 / MB, 2.68, nm, ha='center', va='bottom', fontsize=6.5, color=INK, rotation=90)
for lab, a, b in (('Yp', 0, 10.71), ('Yq', 10.71, 62.46)):   # arms
    ax0.plot([a + 0.15, b - 0.15], [4.0, 4.0], color=INK, lw=0.8); ax0.text((a + b) / 2, 4.1, lab, ha='center', va='bottom', fontsize=8, color=INK, fontweight='bold')
CC = {}   # sequence classes (y 1.0-1.6), one label per class on its widest block, legend at the top right
for a, b, c, col in SC:
    ax0.add_patch(plt.Rectangle((a / MB, 1.0), (b - a) / MB, 0.6, fc=col, ec='white', lw=0.3)); CC.setdefault(NMC[c], col)
fig.legend([plt.Rectangle((0, 0), 1, 1, fc=c) for c in CC.values()], list(CC), loc='upper right', bbox_to_anchor=(0.997, 0.99), ncol=2, fontsize=7.5, frameon=False, title='Sequence class', title_fontsize=8, handlelength=1.6, handleheight=1.1, columnspacing=1.0, labelspacing=0.35)
for i, (p, arms) in enumerate(sorted(PAL.items(), key=lambda kv: np.mean(kv[1]))):   # palindromes (y 0.45-0.7), labels on two rows
    for a, b in arms: ax0.add_patch(plt.Rectangle((a / MB, 0.45), (b - a) / MB, 0.25, fc=INK, lw=0))
    ax0.text(np.mean([v for ab in arms for v in ab]) / MB, 0.32 - 0.32 * (i % 2), p, ha='center', va='top', fontsize=5.5, color=INK)
if KM:   # windows whose Y dose comes from the X + Y keys
    for a, b in XYS: ax0.add_patch(plt.Rectangle((a / MB, -0.62), (b - a) / MB, 0.16, fc=INK, lw=0))
for y, lab in ((2.25, 'Cytoband'), (1.3, 'Sequence class'), (0.55, 'Palindromes')) + (((-0.54, 'Y = 2(X+Y) \u2212 X'),) if KM else ()): ax0.text(-0.5, y, lab, ha='right', va='center', fontsize=6.5, color=INK)
ax0.set_xlim(0, xmax); ax0.set_ylim(-0.75 if KM else -0.45, 4.6); ax0.axis('off')
M = np.array([X[s] for s in rows]); im = ax2.imshow(np.clip(M, 0, 2.49), aspect='auto', cmap=cm, norm=nrm, interpolation='nearest', extent=(0, xmax, n, 0))
nd = np.isnan(M).astype(float); ax2.contourf(np.array(W) / MB + step / 2 / MB, np.arange(n) + 0.5, nd, levels=[0.5, 1.5], colors='none', hatches=['////'])   # no data: hatched white
ax2.set_yticks(np.arange(n) + 0.5); ax2.set_yticklabels([s + (f'  ({NOTE[s]})' if s in NOTE else '') for s in rows], fontsize=6.5)
for t, s in zip(ax2.get_yticklabels(), rows): t.set_color('#8a8984' if s in relrow else INK)   # relatives outside the groups in grey
ax2.set_xlabel('CHM13v2 chrY position (Mb)', fontsize=7.5); ax2.tick_params(labelsize=6.5, length=2); ax2.set_xlim(0, xmax)
for a_, lab, col in ((ay, 'Y dose', 1), (ax, 'X dose', 2)):   # dose bars with numbers
    v = [info[s][col] for s in rows]; a_.barh(np.arange(n) + 0.5, v, height=0.62, color='#8a8984', lw=0)
    for i, x in enumerate(v): a_.text(min(x, 2.25) + 0.05, i + 0.5, f'{x:.2f}', va='center', ha='left', fontsize=5.5, color=INK)
    for g in (1, 2): a_.axvline(g, color='#b9b8b3', lw=0.5, ls=':')
    a_.set_xlim(0, 2.9); a_.set_ylim(n, 0); a_.set_yticks([]); a_.set_xticks([0, 1, 2]); a_.tick_params(labelsize=6, length=2); a_.set_xlabel(lab, fontsize=7)
    for sp_ in ('top', 'right'): a_.spines[sp_].set_visible(False)
axg.set_ylim(n, 0); axg.set_xlim(0, 1); axg.set_autoscale_on(False); axg.axis('off')
for lab, a, b in gl:
    for a_ in (ax2, ay, ax): a_.axhline(a, color=INK, lw=0.6)
    axg.plot([0.08, 0.08], [a + 0.15, b - 0.15], color=INK, lw=0.6); axg.text(0.2, (a + b) / 2, textwrap.fill(lab, 14), va='center', ha='left', fontsize=6.5, color=INK)
cb = fig.colorbar(im, cax=cax, boundaries=BOUND, ticks=[(a + b) / 2 for a, b in zip(BOUND[:-1], BOUND[1:])]); cb.ax.set_yticklabels(['0', '0.5', 'loss', '1', 'gain', '1.5', '2'])
cb.ax.tick_params(labelsize=6); cb.set_label('Copies vs 100 reference XY males; hatched white: no data', fontsize=7)
RULES = ('Heatmap: each 100 kb window of CHM13v2 chrY shows the coverage over the median of 100 reference XY males (Y state normal; 1 = as in the reference; 2 = twice). Where reads map uniquely, coverage is MQ20 reads per chromosome 2 read '
         '(GRCh38 alignments of the 1kGP CRAMs, moved to CHM13v2 with the UCSC hg38ToHs1 chain); in the arrays (centromere, DYZ19, Yq12), where reads do not map uniquely, it is the median count of the yruler marker k-mers that lie in the window. '
         'In Yq12 the value is the array amount relative to the median male, not a chromosome copy number. Colour steps: 0 (< 0.25), 0.5 (0.25-0.75), loss (0.75-0.9), 1 (0.9-1.1), gain (1.1-1.25), 1.5 (1.25-1.75), 2 (>= 1.75). Hatched white: no unique reads and fewer than 50 marker k-mers.\n'
         'Y dose = 2 x mean count of the Y single-copy ruler k-mers / mean count of the chromosome 2 k-mers (Y against an autosome). X dose = chrX MQ20 reads outside the pseudoautosomal regions / chr2 MQ20 reads, over the median of males with Y dose 0.95-1.05.\n'
         'Y state, first rule that applies: no Y call, Y dose <= 0.0222 (highest of 65 heldout females); Yp deletion, Yp dose < 0.5 and Yq dose within 0.1 of 1; Yq loss, Yp dose within 0.1 of 1 and Yq dose more than 0.1 below Yp; Yp \u2260 Yq (state partial in the tables), Yp and Yq doses differ by more than 0.1, so part of the Y is gained and part lost; PAR2 gain, PAR2 dose >= 1.25 with PAR1 within 0.1 of 1; Y gain, Y dose >= 1.75; '
         'mosaic gain, 1.25 to 1.75; mosaic loss, < 0.9. XYY and XXY need both doses: XYY is Y gain with X dose 1, XXY is X dose 2 with Y dose 1 (X dose in 1,567 XY males with Y state normal: 0.981-1.022). '
         'Unusual array estimates: males noted in the cohort analysis. Rows are ordered by average-linkage clustering within each group; parents and offspring that are both in the figure are drawn next to each other.')
if KM:   # alignment-free wording
    rx = [float(r['x_kmer']) for r in XP.values() if r['role'] == 'reference']; nref = len(rx)
    cen = [float(r['ratio']) for r in rd(YP) if r['role'] == 'reference' and r['level'] == 'class' and r['name'] == 'CEN'] if YP else [np.nan]; print(f'centromere class in reference males: {len(cen)}, {min(cen):.3f}-{max(cen):.3f}')
    RULES = RULES.replace(RULES.split('\n')[0], ('Heatmap: each 100 kb window of CHM13v2 chrY shows the median dose of the CHM13v2 chrY 31-mers absent from every other CHM13 sequence (Y-specific, multi-copy ones included) that lie in the window, counted with yruler in the reads that carry a yruler marker k-mer (picked without alignment)' if YP else 'Heatmap: each 100 kb window of CHM13v2 chrY shows the median dose of the yruler marker k-mers that lie in the window, counted in all reads of the 1kGP CRAMs without alignment') + ', over the median of 100 reference XY males (Y state normal; 1 = as in the reference; 2 = twice). '
        + (f'In the centromere (DYZ3), DYZ19 and Yq12 the value is the array amount relative to the median male and depends on the Y lineage (centromere in the 100 reference XY males: {min(cen):.3f}-{max(cen):.3f}), not a chromosome copy number. ' if YP else 'In Yq12 the value is the array amount relative to the median male, not a chromosome copy number. ') + 'Colour steps: 0 (< 0.25), 0.5 (0.25-0.75), loss (0.75-0.9), 1 (0.9-1.1), gain (1.1-1.25), 1.5 (1.25-1.75), 2 (>= 1.75). Hatched white: fewer than 50 k-mer positions in the window (segmental duplications).')
    FEM = [float(r['y_dose']) for r in rd('results/model_final_2026-11-02/cohort_fem.tsv') if r['cohort'] == '1kGP' and r['y_dose'] != 'nan']   # no_Y cut: highest 1kGP female
    RULES = RULES.replace('Y dose <= 0.0222 (highest of 65 heldout females)', f'Y dose <= {max(FEM):.4f} (highest of {len(FEM)} 1kGP females)')
    RULES = RULES.replace(RULES.split('\n')[1], 'Y dose = 2 x mean count of the Y single-copy ruler k-mers / mean count of the chromosome 8 k-mers in chr8:91-101 Mb (Y against an autosome), over the median of the 1kGP males. ' f'X dose = median dose of CHM13 X-specific 31-mers (absent from chrY, autosomes and chrM; 1 in 32 kept), counted with yruler in the chrX, chrY, chromosome 2 and unmapped reads of the CRAM, over the median of {nref} reference XY males. Window, X and PAR doses use the chromosome 8 count of the Y dose as denominator. '
        + ('PAR1 and PAR2' if YP else 'PAR1, PAR2 and the X-transposed region') + ' (black bar under the palindromes): Y dose inferred from CHM13 31-mers present on chrX and chrY and absent from autosomes, Y = 2 x (X + Y dose) - X dose, same reads and reference; a copy change on the X copy of these regions shows as a Y change.')
    if SUBS: RULES = RULES.replace(RULES.split('\n')[0], RULES.split('\n')[0] + f' {", ".join(sorted(SUBS))}: Y from the 31-mer counts of the chrX, chrY, chromosome 2 and unmapped reads of the CRAM (saved reads written before the whole-Y k-mers were added).')
    RULES = RULES.replace('(X dose in 1,567 XY males with Y state normal: 0.981-1.022)', f'(k-mer X dose in the {nref} reference XY males: {min(rx):.3f}-{max(rx):.3f})')
RULES = RULES.replace('Rows are ordered by average-linkage clustering within each group; parents and offspring that are both in the figure are drawn next to each other.', 'Each block is one parent (top) and its offspring (bottom); the label gives the Y state of the unusual member.') if MODE == 'pairs' else RULES
fig.text(0.02, 0.147, '\n'.join(textwrap.fill(l, 285) for l in RULES.split('\n')), ha='left', va='top', fontsize=6.0, color=INK, linespacing=1.3)
fig.savefig(OUT); print(f'written {OUT}')
