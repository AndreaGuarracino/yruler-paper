# Y tree validation: root on A0, haplogroup label monophyly
# usage: python3 ytree_validate.py <tree.nwk> <poznik calls> <haplogrouper calls> <out prefix>
import sys, re, collections
from Bio import Phylo
TR, PZ, HG, OUT = sys.argv[1:5]
t = Phylo.read(TR, 'newick'); tips = {c.name for c in t.get_terminals()}
P = {}   # Poznik long-form ISOGG 2014
for l in open(PZ):
    x = l.split()
    if len(x) >= 6 and x[0] in tips: P[x[0]] = x[-1]
H = {l.split('\t')[0]: l.split('\t')[1] for l in list(open(HG))[1:] if l.split('\t')[0] in tips}   # HaploGrouper ISOGG 2019
a0 = sorted(s for s, h in P.items() if h.startswith('A0')); print('tips', len(tips), 'Poznik', len(P), 'HaploGrouper', len(H), 'A0 outgroup', a0, flush=True)
desc = lambda cl: {c.name for c in cl.get_terminals()}
a0h = sorted(s for s in tips if s in H and H[s].startswith('A0')); a1 = sorted(s for s in tips if s in H and H[s].startswith('A1'))
t.root_with_outgroup(a0h[0]); t.root_with_outgroup(t.common_ancestor(a0h))   # one A0 first, then all A0
print('A0 (HaploGrouper)', a0h, 'monophyletic:', desc(t.common_ancestor(a0h)) == set(a0h), '| A1', len(a1), 'monophyletic:', desc(t.common_ancestor(a1)) == set(a1))
bt = sorted(tips - set(a0h)); print('A1 + BT (all non-A0) form the other root side:', desc(t.common_ancestor(bt)) == set(bt))
tok = lambda h: re.findall(r'[A-Z]+|[a-z]+|\d+', h.split('~')[0].rstrip('*'))
def sup(cl):   # SH-aLRT/UFBoot label of the clade
    m = re.match(r'([\d.]+)/([\d.]+)', str(cl.confidence if cl.confidence is not None else cl.name or ''))
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)
rows = []
for src, L in (('Poznik2016', P), ('HaploGrouper', H)):
    S = collections.defaultdict(set)
    for s, h in L.items():
        tk = tok(h)
        for d in range(1, len(tk) + 1): S[''.join(tk[:d])].add(s)
    for lab, mem in sorted(S.items()):
        if len(mem) < 2 or len(mem) == len(L): continue
        ca = t.common_ancestor(sorted(mem)); dd = desc(ca) & set(L); intr = sorted(dd - mem)
        sh, ub = sup(ca); rows.append((src, lab, len(tok(lab)) if False else len(re.findall(r'[A-Z]+|[a-z]+|\d+', lab)), len(mem), len(intr), sh, ub, ','.join(intr[:10])))
with open(f'{OUT}.labels.tsv', 'w') as o:
    o.write('source\tlabel\tdepth\tcarriers\tintruders\tsh_alrt\tufboot\tintruders_first10\n')
    for r in rows: o.write('\t'.join('' if v is None else str(v) for v in r) + '\n')
for src in ('Poznik2016', 'HaploGrouper'):
    R = [r for r in rows if r[0] == src]; mono = [r for r in R if r[4] == 0]
    print(f'{src}: labels {len(R)}, monophyletic {len(mono)} ({100 * len(mono) / max(1, len(R)):.1f}%)')
    for d in (1, 2, 3, 4, 6, 8):
        Rd = [r for r in R if r[2] == d]
        if Rd: print(f'  depth {d}: {sum(r[4] == 0 for r in Rd)} of {len(Rd)}')
    bad = sorted([r for r in R if r[4] > 0], key=lambda r: -r[4])[:12]
    for r in bad: print(f'  not monophyletic: {r[1]} carriers {r[3]} intruders {r[4]} ({r[7][:80]})')
for s in ('HG03742', 'HG02040'):
    if s in tips:
        cl = t.common_ancestor([s]); par = t.get_path(s)[-2] if len(t.get_path(s)) > 1 else None
        sib = sorted(desc(par) - {s}) if par else []
        print(s, 'Poznik', P.get(s), 'HaploGrouper', H.get(s), 'sister clade size', len(sib), 'sister Poznik labels', collections.Counter(P.get(x, '?')[:3] for x in sib).most_common(4))
j = sorted(s for s, h in H.items() if h.startswith('J2a1a1')); ca = t.common_ancestor(j); extra = sorted(desc(ca) - set(j))
print('J2a1a1 (HaploGrouper):', len(j), 'males; smallest clade has', len(desc(ca)), 'tips; extra', extra)
Phylo.write(t, f'{OUT}.rooted.nwk', 'newick'); print('written', OUT)
