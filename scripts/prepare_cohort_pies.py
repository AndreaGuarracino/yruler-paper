#!/usr/bin/env python3
"""Join cohort manifests, reported coordinates and observed Y lineage calls.

No lineage is inferred from a population label. Missing or conflicting calls
are retained as Unassigned. The 1kGP subset matches the dated snapshot or development plots;
HGDP/SGDP are the planned manifests, before their final sample QC.
The combined map counts shared individuals once, using documented identifiers.
"""

import csv
import hashlib
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote, urlparse

import openpyxl


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data/figure_metadata"
ALIAS_FIELDS = ("Sample_ID", "Sample_ID(Aliases)", "SGDP_ID", "Illumina_ID")
COHORT_PRIORITY = {"1kGP": 0, "HGDP": 1, "SGDP": 2}
HALLAST_URL = "https://static-content.springer.com/esm/art%3A10.1007%2Fs00439-020-02204-9/MediaObjects/439_2020_2204_MOESM2_ESM.xlsx"
HALLAST_SHA = "b8ceb96ffbcba51254d9e86735daa178cf7a2397204a2c902167115798d171b3"


def hallast_book():
    # Hallast 2021 Table S2, fetched once
    path = META / "hallast2021_supplement.xlsx"
    if not path.exists():
        print(f"downloading {HALLAST_URL} -> {path}", file=sys.stderr)
        urllib.request.urlretrieve(HALLAST_URL, path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != HALLAST_SHA:
        sys.exit(f"ERROR: {path} sha256 {digest}, expected {HALLAST_SHA}")
    return openpyxl.load_workbook(path, read_only=True, data_only=True)


def read_tsv(path, encoding="utf-8"):
    with Path(path).open(encoding=encoding) as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def norm(value):
    return "".join(c.lower() for c in str(value) if c.isalnum())


def valid_key(value):
    return value is not None and str(value).strip().lower() not in {"", "na", "n/a", "unknown", "none", "-"}


def major(call):
    call = str(call).strip()
    if call.startswith("NO"):
        return "NO"
    if call.startswith(("YRoot", "CT", "BT", "DE", "CF", "P*", "K*")):
        return "Unassigned"
    return call[0] if call and call[0] in "ABCDEFGHIJKLMNOPQRST" else "Unassigned"


def write_rows(path, rows):
    with Path(path).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def combine_individuals(rows, sgdp_meta):
    """Collapse exact sample IDs and documented individual aliases, never populations."""
    parents = {}

    def root(key):
        parents.setdefault(key, key)
        if parents[key] != key:
            parents[key] = root(parents[key])
        return parents[key]

    for metadata in sgdp_meta:
        aliases = [norm(metadata[key]) for key in ALIAS_FIELDS if valid_key(metadata[key])]
        for alias in aliases[1:]:
            parents[root(alias)] = root(aliases[0])
    individuals = defaultdict(list)
    for row in rows:
        individuals[root(norm(row["sample"].removeprefix("SGDP_")))].append(row)
    combined, overlaps = [], []
    for records in individuals.values():
        records.sort(key=lambda row: (COHORT_PRIORITY[row["cohort"]], row["sample"]))
        if len({row["cohort"] for row in records}) != len(records):
            raise ValueError(f"Multiple records within a cohort resolve to one individual: {records}")
        primary = records[0]
        assigned = [row for row in records if row["lineage"] != "Unassigned"]
        lineages = {row["lineage"] for row in assigned}
        lineage = next(iter(lineages)) if len(lineages) == 1 else "Unassigned"
        call_row = assigned[0] if len(lineages) == 1 else primary
        result = {
            "individual_id": primary["sample"],
            "cohorts": ";".join(row["cohort"] for row in records),
            "cohort_sample_ids": ";".join(f"{row['cohort']}:{row['sample']}" for row in records),
            "population": primary["population"],
            "latitude": primary["latitude"], "longitude": primary["longitude"],
            "lineage": lineage,
            "lineage_source": call_row["lineage_source"] if len(lineages) <= 1 else "conflicting individual assignments",
            "lineage_identifiers": call_row["lineage_identifiers"],
            "coordinate_source": primary["coordinate_source"],
            "coordinate_sample": f"{primary['cohort']}:{primary['sample']}",
            "coordinate_rule": "first available cohort: 1kGP, HGDP, SGDP; reported coordinates retained",
        }
        combined.append(result)
        if len(records) > 1:
            overlaps.append({**result, "all_population_coordinates": ";".join(
                f"{row['cohort']}:{row['population']}:{row['latitude']},{row['longitude']}" for row in records)})
    assert len(combined) + sum(len(row["cohorts"].split(";")) - 1 for row in overlaps) == len(rows)
    return combined, overlaps


def prepare(snapshot_path=None):
    coords = {(row["cohort"], norm(row["population"])): row
              for row in read_tsv(META / "cohort_population_locations.tsv")}
    sgdp_meta = read_tsv(META / "sgdp_metadata.txt", encoding="latin-1")
    sgdp_aliases = defaultdict(list)
    for row in sgdp_meta:
        for key in ALIAS_FIELDS:
            if valid_key(row[key]):
                sgdp_aliases[norm(row[key])].append(row)
    workbook = hallast_book()
    table = list(workbook["Suppl2"].values)
    published = defaultdict(set)
    for fields in table[2:]:
        row = dict(zip(table[1], fields))
        call = major(row["hg_ISOGG_v03.10.19"])
        for key in ("Sample_ID", "vcf_ID"):
            if valid_key(row[key]) and call != "Unassigned":
                published[norm(row[key])].add(call)

    def published_call(keys):
        found = set()
        for key in keys:
            found.update(published.get(norm(key), set()))
        identifiers = ";".join(sorted({str(key) for key in keys if published.get(norm(key))}))
        if len(found) == 1:
            return next(iter(found)), "Hallast 2021 Supplementary Table 2", identifiers
        reason = "no matched published assignment" if not found else "conflicting published assignments"
        return "Unassigned", reason, identifiers

    rows = []
    calls = {row["ID"]: row["Haplogroup"] for row in read_tsv(ROOT / "data/method/haplogroups_1kgp_isogg2019_decode1.tsv")}
    if snapshot_path is None:
        source_rows = [row for row in read_tsv(ROOT / "data/catalog_d376/out_rule/catalog.tsv")
                       if row["sample"] not in {"HG02300", "NA19226"} and float(row["y_dosage"]) >= 0.9]
        stage = "development catalog after dosage filter"
    else:
        source_rows = [row for row in read_tsv(snapshot_path)
                       if row["flag"] in {"ok", "y_dose_below_0.9"} and row["sex_ped"] == "1"]
        calls = {row["sample"]: row["haplogroup"] for row in source_rows}
        stage = f"analysed frozen-model snapshot: {Path(snapshot_path).name}"
    for row in source_rows:
        coordinate = next(value for (cohort, _), value in coords.items()
                          if cohort == "1kGP" and value["population_id"] == row["population"])
        rows.append({"cohort": "1kGP", "sample": row["sample"], "population": row["population"],
                     "latitude": coordinate["latitude"], "longitude": coordinate["longitude"],
                     "lineage": major(calls[row["sample"]]), "lineage_source": "HaploGrouper ISOGG 2019 deCODE",
                     "lineage_identifiers": row["sample"],
                     "coordinate_source": "IGSR population metadata", "sample_stage": stage})

    missing_coordinates = []
    for line in (META / "planned_hgdp_sgdp_samples.tsv").read_text().splitlines():
        sample, url, _ = line.split("\t")
        cohort = "SGDP" if sample.startswith("SGDP_") else "HGDP"
        population = unquote(urlparse(url).path.split("/data/", 1)[1].split("/", 1)[0])
        identifier = sample.removeprefix("SGDP_")
        keys = [identifier]
        for row in sgdp_aliases.get(norm(identifier), []) if cohort == "SGDP" else []:
            # A matching individual supplies additional exact identifiers.
            keys.extend(row[key] for key in ALIAS_FIELDS if valid_key(row[key]))
        lineage, source, identifiers = published_call(keys)
        coordinate = coords.get((cohort, norm(population)))
        coordinate_source = "IGSR population metadata"
        if coordinate is None and cohort == "SGDP":
            # A reported public population location fills a missing IGSR row.
            candidates = [row for row in sgdp_meta if norm(row["Population_ID"]) == norm(population)
                          and row["Embargo"] == "FullyPublic"]
            values = {(row["Latitude"], row["Longitude"]) for row in candidates}
            if len(values) == 1:
                lat, lon = next(iter(values))
                coordinate = {"latitude": lat, "longitude": lon}
                coordinate_source = "SGDP published public-sample metadata"
        if coordinate is None:
            missing_coordinates.append(sample)
            continue
        rows.append({"cohort": cohort, "sample": sample, "population": population,
                     "latitude": coordinate["latitude"], "longitude": coordinate["longitude"],
                     "lineage": lineage, "lineage_source": source,
                     "lineage_identifiers": identifiers,
                     "coordinate_source": coordinate_source, "sample_stage": "planned manifest before final QC"})
    if missing_coordinates:
        raise ValueError(f"Missing reported coordinates for {missing_coordinates}")
    assert len({(row["cohort"], row["sample"]) for row in rows}) == len(rows)
    write_rows(META / "cohort_pie_samples.tsv", rows)
    combined, overlaps = combine_individuals(rows, sgdp_meta)
    write_rows(META / "cohort_combined_samples.tsv", combined)
    write_rows(META / "cohort_shared_individuals.tsv", overlaps)
    combined_groups = defaultdict(list)
    for row in combined:
        combined_groups[(round(float(row["latitude"]), 6), round(float(row["longitude"]), 6))].append(row)
    combined_counts = []
    for (latitude, longitude), individuals in sorted(combined_groups.items()):
        counts = Counter(row["lineage"] for row in individuals)
        populations = ";".join(sorted({row["population"] for row in individuals}))
        for lineage, n in sorted(counts.items()):
            combined_counts.append({"latitude": latitude, "longitude": longitude,
                                    "population_labels": populations, "lineage": lineage,
                                    "n": n, "location_n": len(individuals)})
    write_rows(META / "cohort_combined_counts.tsv", combined_counts)
    groups = defaultdict(Counter)
    for row in rows:
        groups[(row["cohort"], row["population"], row["latitude"], row["longitude"])][row["lineage"]] += 1
    with (META / "cohort_pie_counts.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["cohort", "population", "latitude", "longitude", "lineage", "n", "population_n"])
        for (cohort, population, lat, lon), counts in sorted(groups.items()):
            for lineage, n in sorted(counts.items()):
                writer.writerow([cohort, population, lat, lon, lineage, n, sum(counts.values())])
    for cohort in ("1kGP", "HGDP", "SGDP"):
        selected = [row for row in rows if row["cohort"] == cohort]
        print(cohort, "samples", len(selected), "population records", len({row["population"] for row in selected}),
              "with lineage", sum(row["lineage"] != "Unassigned" for row in selected),
              "unassigned", sum(row["lineage"] == "Unassigned" for row in selected))
    print("Combined unique individuals", len(combined), "shared individuals", len(overlaps),
          "reported locations", len(combined_groups),
          "unassigned", sum(row["lineage"] == "Unassigned" for row in combined))


if __name__ == "__main__":
    prepare()
