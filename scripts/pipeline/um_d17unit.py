# DYZ17 full repeat unit: start-to-start units of DYZ17 consensus hits -> MAFFT -> hmmbuild -> consensus
# usage: python3 um_d17unit.py [min gap 5400] [max gap 5700]; 10 random units per assembly (seed 7), assemblies with a complete DYZ17 region
import sys, os, gzip, random, subprocess as sp
P = os.environ['P']; SAM, MAFFT, HB, HE = 'samtools', 'mafft', 'hmmbuild', 'hmmemit'   # samtools 1.24, MAFFT 7.526, HMMER 3.4
G0, G1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (5400, 5700); O = f'{P}/hmm_d17unit'; os.makedirs(O, exist_ok=True); random.seed(7)
def log(*a): print(*a, flush=True)
A = [l.split('\t')[0] for l in list(open(f'{P}/newclass/truth_all.tsv'))[1:] if l.split('\t')[0] != 'CHM13']; pick = []; nu = 0
for a in A:
    R = sorted((l.split()[0], min(int(l.split()[8]), int(l.split()[9])), max(int(l.split()[8]), int(l.split()[9])), l.split()[11]) for l in gzip.open(f'{P}/hallast_ours/{a}/{a}.ours.DYZ17_ours.hmmer-tblout.txt.gz', 'rt') if not l.startswith('#'))
    U = []
    for x, y in zip(R, R[1:]):
        if x[0] != y[0] or x[3] != y[3]: continue
        if x[3] == '+' and G0 <= y[1] - x[1] <= G1: U.append((x[0], x[1], y[1] - 1, '+'))   # motif start to next motif start
        if x[3] == '-' and G0 <= y[2] - x[2] <= G1: U.append((x[0], x[2] + 1, y[2], '-'))   # minus: motif start is the high end
    nu += len(U); pick += [(a,) + u for u in random.sample(U, min(10, len(U)))]
log(f'assemblies {len(A)}, start-to-start units with gap {G0}-{G1}: {nu}, picked {len(pick)}')
with open(f'{O}/DYZ17.units.fa', 'w') as o:
    for a, c, s, e, st in pick:
        r = sp.run([SAM, 'faidx', f'{P}/ks_truth/fa/{a}.fa.gz', f'{c}:{s}-{e}'] + (['-i'] if st == '-' else []), capture_output=True, text=True, check=True)
        o.write(f'>{a}:{s}-{e}{st}\n' + ''.join(r.stdout.split('\n')[1:]) + '\n')
with open(f'{O}/DYZ17.msa.fa', 'w') as o: sp.run([MAFFT, '--auto', '--thread', os.environ.get('SLURM_CPUS_PER_TASK', '8'), '--quiet', f'{O}/DYZ17.units.fa'], stdout=o, check=True)
r = sp.run([HB, '--dna', '-n', 'DYZ17', f'{O}/DYZ17.hmm', f'{O}/DYZ17.msa.fa'], capture_output=True, text=True, check=True); log([l for l in r.stdout.split('\n') if l and not l.startswith('#')][0])
c = sp.run([HE, '-c', f'{O}/DYZ17.hmm'], capture_output=True, text=True, check=True).stdout.split('\n'); q = ''.join(c[1:])
open(f'{P}/hallast_motifs/DYZ17_unit.fasta', 'w').write(f'>DYZ17_unit\n{q}\n'); log(f'consensus DYZ17_unit {len(q)} bp -> hallast_motifs/DYZ17_unit.fasta')
