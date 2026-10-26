# Superset: class k-mers in >= 40% of class assemblies, minus CHM13 non-Y; plus old marker sets
import subprocess as sp, collections, math, os, sys
P = '/path/to/project'; O = f'{P}/markers_u'; D = '/path/to/scratch/umk/db'; K = '/path/to/conda/kmc/3.2.4/bin'
BG = f'{P}/ext/kmers/db/bg_noY'; FR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
def run(*a): sp.run(list(a), check=True, stdout=sp.DEVNULL)
A = collections.defaultdict(list)
for l in open(f'{O}/regions.tsv'): a, c = l.split('\t')[:2]; A[c].append(a) if a not in A[c] else None
rc = str.maketrans('ACGT', 'TGCA'); can = lambda k: min(k, k.translate(rc)[::-1])
S = {}
for c in sorted(A):
    acc = None
    for a in A[c]:
        db = f'{D}/{a}.{c}'
        if not os.path.exists(db + '.kmc_pre'): raise SystemExit(f'ERROR missing {db}')
        if acc is None:
            acc = f'{D}/sum.{c}'
            for e in ('pre', 'suf'): sp.run(['cp', f'{db}.kmc_{e}', f'{acc}.kmc_{e}'], check=True)
            continue
        run(f'{K}/kmc_tools', '-t4', 'simple', acc, db, 'union', f'{D}/sum2.{c}', '-ocsum', '-cs65535')
        for e in ('pre', 'suf'): os.replace(f'{D}/sum2.{c}.kmc_{e}', f'{acc}.kmc_{e}')
    need = math.ceil(FR * len(A[c]))
    run(f'{K}/kmc_tools', '-t4', 'simple', acc, f'-ci{need}', BG, 'kmers_subtract', f'{D}/sup.{c}')
    run(f'{K}/kmc_tools', '-t4', 'transform', f'{D}/sup.{c}', 'dump', f'{O}/sup.{c}.txt')
    run(f'{K}/kmc_tools', '-t4', 'transform', acc, 'dump', f'{O}/presence.{c}.txt')   # all k-mers with presence count
    S[c] = [l.split('\t')[0] for l in open(f'{O}/sup.{c}.txt')]
    print(f'{c}: assemblies {len(A[c])}, presence >= {need}, superset {len(S[c])}', flush=True)
old = collections.defaultdict(set)   # old sets: 90% classes, 50% classes, training 4
cl = None
for l in open(f'{P}/yruler_markers/markers_classes.fa'):
    if l[0] == '>': cl = l[1:].strip().rsplit('_', 1)[0]
    else: old[cl].add(can(l.strip()))
for f, c in [('DYZ1', 'dyz1'), ('DYZ2', 'dyz2'), ('DYZ18', 'dyz18'), ('Yqhet_3k1bp', 'yq3k1'), ('Yqhet_2k7bp', 'yq2k7'), ('DYZ19', 'dyz19')]:
    old[c] |= {can(l.split()[0]) for l in open(f'{P}/class_markers_f50/{f}.k31.txt')}
with open(f'{O}/markers_u.fa', 'w') as o:
    for c in ['dyz3', 'yq12', 'ysinglecopy', 'chr2depth', 'dyz1', 'dyz2', 'dyz18', 'yq3k1', 'yq2k7', 'dyz19']:
        ks = sorted(set(S.get(c, [])) | old[c]); n = 0
        for k in ks: n += 1; o.write(f'>{c}_{n}\n{k}\n')
        print(f'fasta {c}: superset {len(S.get(c, []))}, old {len(old[c])}, union {len(ks)}', flush=True)
