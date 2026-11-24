# Read mixtures summary against whole genomes
# usage: [EST=<estimates file name, default estimates.tsv>] python3 mixtures.py <cohort_1kgp.tsv> <mixreads dir> <out.tsv> <summary out.tsv>; change = 100 x (run length / cohort length of the man - 1), ruler runs only; fraction 0: bp of the female
import sys, os, csv, glob, re, collections as co, numpy as np
COH, D, OUT, SUM = sys.argv[1:5]; C = {r['sample']: r for r in csv.DictReader(open(COH), delimiter='\t')}; CL = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}
num = lambda s: float('0' + s if s[0] == '.' else s); R = []; miss = 0; OLD = {'ruler': 'y_found', 'haploid_depth': 'no_y'}   # old estimate_basis values
for p in sorted(glob.glob(f'{D}/*/f*_d*_s*')):
    m, f = os.path.basename(os.path.dirname(p)).split('_'); fr, dp, s = re.match(r'f(.+)_d(.+)_s(\d+)', os.path.basename(p)).groups()
    lf = next((f'{p}/{n}' for n in (os.environ.get('EST', 'estimates.tsv'), 'lengths.tsv') if os.path.exists(f'{p}/{n}')), None)   # new name, else old
    if lf is None: miss += 1; print(f'WARNING no estimates.tsv or lengths.tsv {p}', flush=True); continue
    L = open(f'{p}/summary.tsv').read().split('\n'); gc = 'gc_curve=0' not in L[0]; nr = next(l.split('\t')[1] for l in L if l and l[0] != '#' and not l.startswith('sample\t'))
    for r in csv.DictReader(open(lf), delimiter='\t'):
        if r['class'] not in CL: continue
        v = r['estimate_bp'] if 'estimate_bp' in r else r['length_bp']; e = float('nan') if v == 'NA' else float(v);   # NA: no_y runs of the NA-rule binary
        b = float(C[m][CL[r['class']] + '_bp']); R.append((m, f, num(fr), num(dp), int(s), CL[r['class']], e, b, 100 * (e / b - 1), float(r['y_dose']), r['y_status'] if 'y_status' in r else OLD[r['estimate_basis']], gc, nr))
print(f'{len(R)} rows, {len({r[:5] for r in R})} runs, missing {miss}', flush=True)
with open(OUT, 'w') as o:
    o.write('male\tfemale\tfraction\tdepth\tseed\tclass\tlength_bp\twhole_bp\tchange_pct\ty_dose\tbasis\tgc_curve\tnreads\n'); o.writelines('\t'.join(map(str, r)) + '\n' for r in R)
G = co.defaultdict(list)
for r in R: G[r[2], r[3], r[5]].append(r)
with open(SUM, 'w') as o:
    o.write('fraction\tdepth\tclass\truns\truns_ruler\tmedian_abs_change_pct\tp90_abs_change_pct\tmax_abs_change_pct\tmedian_signed_pct\tmedian_bp\tmax_bp\ty_dose_median\ty_dose_min\ty_dose_max\tshare_ruler\tshare_gc_curve\n')
    for k in sorted(G, reverse=True):
        v = G[k]; ch = np.array([x[8] for x in v if x[10] == 'y_found']); a = np.abs(ch); y = np.array([x[9] for x in v]); bp = np.array([x[6] for x in v]); bp = bp[~np.isnan(bp)] if (~np.isnan(bp)).any() else np.array([np.nan]); st = f'{np.median(a):.2f}\t{np.percentile(a, 90):.2f}\t{a.max():.2f}\t{np.median(ch):+.2f}' if k[0] > 0 and len(ch) else 'nan\tnan\tnan\tnan'
        o.write(f'{k[0]}\t{k[1]}\t{k[2]}\t{len(v)}\t{len(ch)}\t{st}\t{np.median(bp):.0f}\t{bp.max():.0f}\t{np.median(y):.4f}\t{y.min():.4f}\t{y.max():.4f}\t{np.mean([x[10] == "y_found" for x in v]):.3f}\t{np.mean([x[11] for x in v]):.3f}\n')
print(f'written {OUT} and {SUM}', flush=True)
