# Depth and lineage-out tables for the figure scripts
# usage: [TAG=v9b] python3 um_v9export.py <out dir> [sum dir, default model_v5_v9b_sum] [marker dir, default markers_v9]
# errors: um_v5sum units (Yq12 = DYZ1 + DYZ2 + rest, DYZ2 with its G+C floor); fold, choice, markers, target bp: per-class CV files (DYZ2: the floor run 'ne');
# chr2_depth_test: chr2depth_mean of the test run at that level; level 1 = cv_units.tsv (the CV of cv.tsv) with the full-depth runs
import os, sys, glob, re, collections
P = os.environ['P']; OUT = sys.argv[1]; SD = f'{P}/{sys.argv[2] if len(sys.argv) > 2 else "model_v5_v9b_sum"}'; U = f'{P}/{sys.argv[3] if len(sys.argv) > 3 else "markers_v9"}'; os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get('TAG', 'v9b'); CL = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17']; LV = {'1': 'full', 'd30x': 'd30x', 'd25x': 'd25x', 'd20x': 'd20x', 'd15x': 'd15x', 'd10x': 'd10x', 'd5x': 'd5x', 'd2x': 'd2x', 'd1x': 'd1x'}
def log(*a): print(*a, flush=True)
def units(f):
    R = {}
    for l in list(open(f))[1:]:
        c, r, m, e = l.rstrip('\n').split('\t')
        if c in CL: R[(c, r, m)] = e
    return R
def rows(pat, key):   # (class, key fields, sample) -> row of the per-class CV files
    R = {}
    for c in CL:
        d = pat('ne' if c == 'dyz2' else '')
        fs = [f for f in glob.glob(f'{d}/*_c{c}.tsv')]
        if not fs: sys.exit(f'ERROR no files {d}/*_c{c}.tsv')
        for f in fs:
            for l in open(f):
                x = l.rstrip('\n').split('\t'); R[(c, key(x), x[3] if key(x) is not None else x[2])] = x
    return R
DEP = {}
def dep(lv, s):   # chr2 depth-marker mean count of the test run
    if (lv, s) not in DEP:
        g = f'{U}/depth/{lv}/{s}/run.log'
        if not os.path.exists(g): g = f'{U}/depth/full/{s}/run.log'; log(f'{s}: no {lv} run (full depth below the level); the CV used its full-depth counts (um_eval ALTDIR)')
        m = re.search(r'chr2depth_mean ([0-9.eE+-]+)', open(g).read()) if os.path.exists(g) else None
        if not m: sys.exit(f'ERROR no chr2depth_mean in {g}')
        DEP[(lv, s)] = m.group(1)
    return DEP[(lv, s)]
for f, lv in LV.items():
    E = units(f'{SD}/cv_units.tsv' if f == '1' else f'{SD}/cv_depth_{lv}_units.tsv')
    X = rows(lambda ne: f'{U}/cv_v5{TAG}{ne}' if f == '1' else f'{U}/cv_v5{TAG}{ne}_depth_{lv}', lambda x: x[0])
    n = collections.Counter(); miss = 0
    with open(f'{OUT}/depth_f{f}.tsv', 'w') as o:
        o.write('repeat\tfold\tclass\tsample\tchoice_presence_outside_ruler_filter\tmarkers\terror_pct\treference_bp\testimate_bp\tchr2_depth_test\n')
        for (c, r, m), e in sorted(E.items(), key=lambda z: (CL.index(z[0][0]), int(z[0][1]), z[0][2])):
            x = X.get((c, r, m))
            if x is None: miss += 1; continue
            t = float(x[7]); o.write('\t'.join([r, x[1], c, m, x[4], x[5], e, f'{t:.1f}', f'{t * (1 + float(e) / 100):.0f}', dep(lv, m)]) + '\n'); n[c] += 1
    if miss: sys.exit(f'ERROR level {lv}: {miss} predictions without a CV row')
    log(f'level {f} ({lv}): ' + ', '.join(f'{c} {n[c]}' for c in CL) + f'; men {len({k[2] for k in E})}')
E = units(f'{SD}/lolo_units.tsv'); X = rows(lambda ne: f'{U}/lolo_v5{TAG}{ne}', lambda x: None)
with open(f'{OUT}/lineage_out.tsv', 'w') as o:
    o.write('lineage\tclass\tsample\tchoice_presence_outside_ruler_filter\tmarkers\terror_pct\treference_bp\testimate_bp\n'); k = 0
    for (c, r, m), e in sorted(E.items(), key=lambda z: (CL.index(z[0][0]), z[0][2])):
        x = X.get((c, None, m))
        if x is None: sys.exit(f'ERROR lineage-out {c} {m}: no row')
        t = float(x[6]); o.write('\t'.join([x[0], c, m, x[3], x[4], e, f'{t:.1f}', f'{t * (1 + float(e) / 100):.0f}']) + '\n'); k += 1
log(f'lineage_out: {k} rows, lineages {sorted({x[0] for x in X.values()})}')
