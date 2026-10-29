HGDP/SGDP chrY call integration, 28 October 2026

Input: results/unassigned_lineages_calls_2026-10-28.tsv (21 called, 2 without
Y signal, 8 pending chrY genotypes). Calls are used exactly as delivered;
scores are not interpreted as assignment probabilities. All 21 calls join
existing map records by exact cohort/sample ID; existing assigned broad
lineages are preserved. Full ISOGG calls, scores and site counts are retained
in lineage_call_audit.tsv. Source metadata remain read only.

The paternal-lineage map excludes 8 SGDP records: 7 females in public metadata
(including Peru60, also lacking Y signal), plus NOR111, of unspecified reported
sex and lacking Y signal. The 6 other documented females still have pending
chrY genotypes; they are not represented as unresolved Y lineages in the map.
These exclusions concern the map, not the 1kGP satellite measurements.

Retained cohort records: 1,593 1kGP, 552 HGDP, 175 SGDP; total 2,320.
After the existing exact-ID/documented-alias deduplication: 2,298 individuals
at 132 reported population coordinates; 22 individuals shared between cohorts.
Only HGDP00119 and SGDP_HGDP00725 remain Unassigned in the map.

cohort_pie_samples_all.tsv and cohort_combined_samples_all.tsv retain every
original record/individual for complete QC joins, including excluded samples.
cohort_pie_samples.tsv and cohort_combined_samples.tsv are map-eligible tables.
map_exclusions.tsv records the separate metadata and no-Y-signal reasons.
cohort_combined_counts.tsv gives per-location/lineage sector counts.
manifest.json records inputs, hashes, inclusion policy and all summary counts.

Rebuild: python scripts/cohort_lineages.py
Plot consumers use the same loader in scripts/cohort_lineages.py; the partial
cross-cohort analysis uses the complete join and explicit eligibility flags.
The new calls do not overlap the unassigned assembly-reference individuals.
