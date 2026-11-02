#!/usr/bin/env python3
"""Compare array-reference and CV representation with analysed Y lineages.

Uses the archived ISOGG classification, not the exploratory SNP topology.
Counts describe representation and do not establish coverage of genetic diversity.
"""
import hashlib
import json
from collections import Counter
from pathlib import Path

import frozen_snapshot as f
import validation_plots as v


def main():
    root = f.ROOT
    parents, _ = f.tree_data(root / 'data/figure_metadata/isogg2019_tree.tsv')
    calls = {r['ID']: r['Haplogroup'] for r in f.read_tsv(root / 'data/method/haplogroups_1kgp_isogg2019_decode1.tsv')}
    snapshot = root / 'results/estimates_1kgp.tsv'
    _, population, _ = f.load_snapshot(snapshot)
    for row in population:
        assert calls[row['sample']] == row['haplogroup']
    def group(sample):
        return f.group_of(calls[sample], parents) if sample in calls else 'Unassigned'
    regions = f.read_tsv(root / 'data/markers_u/truth_regions.tsv')
    refs = {k: {r['male'] for r in regions if r['class'] == k} for k in ('dyz3', 'yq12')}
    for k in ('dyz19', 'dyz17'): refs[k] = {r['assembly'] for r in f.read_tsv(root / f'data/method/{k}_clean.tsv') if r['clean'] == 'True'}
    print('assembly references', {k: len(x) for k, x in refs.items()}, flush=True)
    rule = root / 'data/method/truth_qc/x/rule_w4.txt'   # all four assembly QC tracks flag the same 1 kb
    for line in rule.read_text().splitlines():
        target, sample = line.split()[:2]
        if target in refs:
            assert sample in refs[target]
            refs[target].remove(sample)
    print('after the QC rule', {k: len(x) for k, x in refs.items()}, flush=True)
    cv_source = root / 'results/model/depth/depth_f1.tsv'
    predictions = [r for r in f.read_tsv(cv_source) if r['repeat'] == '0']
    files = [cv_source, rule]
    evaluated = {k: {r['sample'] for r in predictions if r['class'] == k} for k in v.KEYS}
    assert evaluated['dyz2'] == evaluated['yq12']
    for key in ('dyz1', 'dyz18s'):   # HG01261 keeps these references
        assert evaluated[key] == evaluated['yq12'] | {'HG01261'}
    for key in refs:
        assert evaluated[key] <= refs[key]
    cohorts = Counter(group(r['sample']) for r in population)
    reference_counts = {k: Counter(map(group, ids)) for k, ids in refs.items()}
    cv_counts = {k: Counter(map(group, evaluated[k])) for k in refs}
    rank = {g: i for i, g in enumerate(f.ORDER)}
    groups = sorted(set(cohorts).union(*(set(c) for c in reference_counts.values())), key=lambda g: (rank.get(g, len(rank)), g))
    rows = []
    for lineage in groups:
        row = {'Y_lineage': lineage, 'cohort_individuals': cohorts[lineage]}
        for key in refs:
            row[key + '_assembly_references'] = reference_counts[key][lineage]
            row[key + '_CV_individuals'] = cv_counts[key][lineage]
        rows.append(row)
    f.write_rows(root / 'results/reference_lineage_coverage.tsv', rows)
    sample_rows = []
    all_ids = {r['sample'] for r in population}.union(*refs.values())
    pop_ids = {r['sample'] for r in population}
    for sample in sorted(all_ids):
        row = {'sample': sample, 'haplogroup': calls.get(sample, ''), 'Y_lineage': group(sample), 'in_analysed_cohort': sample in pop_ids}
        for key in refs:
            row[key + '_assembly_reference'] = sample in refs[key]
            row[key + '_CV_individual'] = sample in evaluated[key]
        sample_rows.append(row)
    f.write_rows(root / 'results/reference_lineage_coverage_samples.tsv', sample_rows)
    tex = [r'\begin{table}[p]', r'\caption{Y-lineage representation of the reference panels and analysed cohort. Each target column gives eligible assembly references / individuals evaluated in cross-validation. The Yq12 column also describes DYZ2; DYZ1 and DYZ18S add HG01261 (lineage R1b), whose Yq12 and DYZ2 lengths are excluded by the assembly-quality rule. Counts use the archived ISOGG 2019 classification hierarchy, independently of the inferred SNP tree. Unassigned denotes missing calls in the archived lookup, not a biological lineage. Reference individuals can also occur in the cohort; columns are not disjoint. A nonzero count does not establish adequate coverage of a lineage, and zero does not establish absence of a related reference in an ancestral or neighbouring clade.}\label{tab:reference-lineages}', r'\centering\small', r'\begin{tabular}{lrrrrr}', r'\toprule', r'Y lineage & Cohort & DYZ3 & Yq12 & DYZ19 & DYZ17 \\', r'\midrule']
    for row in rows:
        cells = [row['Y_lineage'], str(row['cohort_individuals'])]
        cells += [f"{row[k+'_assembly_references']} / {row[k+'_CV_individuals']}" for k in refs]
        tex.append(' & '.join(cells) + r' \\')
    tex += [r'\midrule', 'Total & ' + str(len(population)) + ' & ' + ' & '.join(f'{len(refs[k])} / {len(evaluated[k])}' for k in refs) + r' \\', r'\botrule', r'\end{tabular}', r'\end{table}']
    (root / 'manuscript/reference-lineage-table.tex').write_text('\n'.join(tex) + '\n')
    sources = [snapshot, root / 'data/markers_u/truth_regions.tsv', root / 'data/method/dyz19_clean.tsv', root / 'data/method/dyz17_clean.tsv', root / 'data/method/haplogroups_1kgp_isogg2019_decode1.tsv', root / 'data/figure_metadata/isogg2019_tree.tsv', Path(__file__)] + files
    report = {'population_n': len(population), 'groups_n': len(groups), 'unknown_reference_calls': sorted(s for s in set().union(*refs.values()) if s not in calls), 'cohort_groups_without_evaluated_individuals': {k: [g for g in groups if cohorts[g] and not cv_counts[k][g]] for k in refs}, 'scope': 'Descriptive class-specific lineage representation; no ancestry inference, phylogenetic distances, new model evaluation or new array-QC filtering.', 'inputs': [{'path': str(p.relative_to(root)), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]}
    (root / 'results/reference_lineage_coverage_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: val for k, val in report.items() if k != 'inputs'}, indent=2))


if __name__ == '__main__':
    main()
