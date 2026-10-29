#!/usr/bin/env python3
"""Overlay delivered chrY calls on archived cohort metadata, retaining provenance.

Source metadata are read only. Known SGDP females and samples reported to lack
Y signal are excluded from the paternal-lineage map, with an exclusion audit.
The full metadata join remains available for QC of those samples elsewhere.
"""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from prepare_cohort_pies import (ROOT, META, ALIAS_FIELDS, read_tsv, write_rows,
                                norm, valid_key, major, combine_individuals)

CALLS = ROOT / 'results/unassigned_lineages_calls_2026-10-28.tsv'
OUTPUT = ROOT / 'results/cohort_metadata_2026-10-28'


def build_metadata(calls_path=CALLS):
    base = read_tsv(META / 'cohort_pie_samples.tsv')
    sgdp = read_tsv(META / 'sgdp_metadata.txt', encoding='latin-1')
    calls = read_tsv(calls_path)
    by_call = {(r['cohort'], r['individual_id']): r for r in calls}
    assert len(by_call) == len(calls), 'Duplicate delivered call IDs.'
    available = {(r['cohort'], r['sample']) for r in base}
    assert set(by_call) <= available, 'Delivered IDs absent from the archived map manifest.'
    genders = {}
    for r in sgdp:
        for key in ALIAS_FIELDS:
            if valid_key(r[key]):
                genders.setdefault(norm(r[key]), set()).add(r['Gender'])
    updated, excluded, audit = [], [], []
    for original in base:
        row = dict(original)
        row.update(haplogroup_isogg2019='', genotype_call_status='not supplemented',
                   chrY_sites_genotyped='', net_score='', reported_sex='',
                   map_included=True, map_exclusion_reason='')
        if row['cohort'] == 'SGDP':
            found = genders.get(norm(row['sample'].removeprefix('SGDP_')), set())
            assert len(found) <= 1, ('Conflicting SGDP sex metadata', row['sample'])
            row['reported_sex'] = next(iter(found), 'not reported')
        call = by_call.get((row['cohort'], row['sample']))
        if call:
            row.update(genotype_call_status=call['status'],
                       chrY_sites_genotyped=call['chrY_sites_genotyped'], net_score=call['net_score'])
            if call['status'] == 'called':
                haplogroup = call['haplogroup_isogg2019']
                assert haplogroup not in {'', 'NA'} and major(haplogroup) != 'Unassigned'
                lineage = major(haplogroup)
                assert row['lineage'] in {'Unassigned', lineage}, ('Conflicting existing lineage', row, call)
                assert row['reported_sex'] != 'F', ('Called Y but recorded female; review required', row)
                row.update(lineage=lineage, haplogroup_isogg2019=haplogroup,
                           lineage_source='HaploGrouper ISOGG 2019; delivered chrY genotype calls 2026-10-28',
                           lineage_identifiers=call['individual_id'])
        reasons = []
        if row['reported_sex'] == 'F':
            reasons.append('female in public SGDP metadata')
        if call and call['status'].startswith('no Y signal'):
            reasons.append('no Y signal in delivered chrY audit')
        if reasons:
            row.update(map_included=False, map_exclusion_reason='; '.join(reasons))
            excluded.append(row)
        if call:
            audit.append(dict(call, previous_lineage=original['lineage'], updated_lineage=row['lineage'],
                              reported_sex=row['reported_sex'], map_included=row['map_included'],
                              map_exclusion_reason=row['map_exclusion_reason']))
        updated.append(row)
    included = [r for r in updated if r['map_included']]
    combined, overlaps = combine_individuals(included, sgdp)
    all_combined, _ = combine_individuals(updated, sgdp)
    detail = {(r['cohort'], r['sample']): r for r in updated}
    for row in all_combined:
        records = [detail[tuple(item.split(':', 1))] for item in row['cohort_sample_ids'].split(';')]
        row['map_included'] = any(r['map_included'] for r in records)
        row['map_exclusion_reason'] = '; '.join(sorted({r['map_exclusion_reason'] for r in records if r['map_exclusion_reason']}))
    summary = dict(input_call_records=len(calls),
                   call_status_counts=dict(Counter(r['status'] for r in calls)),
                   assigned_by_cohort=dict(Counter(r['cohort'] for r in calls if r['status'] == 'called')),
                   original_cohort_records=len(base), retained_cohort_records=len(included),
                   retained_by_cohort=dict(Counter(r['cohort'] for r in included)),
                   retained_assigned_by_cohort=dict(Counter(r['cohort'] for r in included if r['lineage'] != 'Unassigned')),
                   excluded_map_records=len(excluded), excluded_samples=[r['sample'] for r in excluded],
                   unique_map_individuals=len(combined), shared_individuals=len(overlaps),
                   reported_locations=len({(round(float(r['latitude']), 6), round(float(r['longitude']), 6)) for r in combined}),
                   lineage_counts=dict(Counter(r['lineage'] for r in combined)),
                   remaining_unassigned=[r['individual_id'] for r in combined if r['lineage'] == 'Unassigned'])
    return dict(updated=updated, included=included, combined=combined, all_combined=all_combined,
                overlaps=overlaps, excluded=excluded, call_audit=audit, summary=summary)


def load_combined_samples(include_excluded=False):
    return build_metadata()['all_combined' if include_excluded else 'combined']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calls', type=Path, default=CALLS)
    p.add_argument('--out', type=Path, default=OUTPUT)
    args = p.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    data = build_metadata(args.calls)
    for key, name in [('updated', 'cohort_pie_samples_all'), ('included', 'cohort_pie_samples'),
                      ('combined', 'cohort_combined_samples'), ('all_combined', 'cohort_combined_samples_all'),
                      ('overlaps', 'cohort_shared_individuals'), ('excluded', 'map_exclusions'),
                      ('call_audit', 'lineage_call_audit')]:
        write_rows(args.out / (name + '.tsv'), data[key])
    locations = defaultdict(list)
    for row in data['combined']:
        locations[round(float(row['latitude']), 6), round(float(row['longitude']), 6)].append(row)
    counts = []
    for (lat, lon), rows in sorted(locations.items()):
        for lineage, n in sorted(Counter(r['lineage'] for r in rows).items()):
            counts.append(dict(latitude=lat, longitude=lon,
                               population_labels=';'.join(sorted({r['population'] for r in rows})),
                               lineage=lineage, n=n, location_n=len(rows)))
    write_rows(args.out / 'cohort_combined_counts.tsv', counts)
    manifest = dict(summary=data['summary'],
                    policy='No lineage inferred from geography. Supplement existing missing calls by exact cohort/sample ID. Exclude documented SGDP females and delivered no-Y-signal samples from the paternal-lineage map only. Preserve complete joins for QC.',
                    inputs=[dict(path=str(path.absolute().relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                            for path in [META/'cohort_pie_samples.tsv', META/'sgdp_metadata.txt', args.calls, Path(__file__).resolve()]],
                    satellite_estimates_changed=False, SNP_tree_changed=False)
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(data['summary'], indent=2))


if __name__ == '__main__':
    main()
