from pathlib import Path
from collections import Counter
import argparse
import csv
import hashlib
import json
import math


SETTINGS = ["Hospital", "WWTP", "Community", "Wet.market"]


def read_unique(path, required):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty input: {path}")
    keys = [row["arg_subtype"].strip() for row in rows]
    assert len(keys) == len(set(keys)), f"Duplicate subtype: {path}"
    assert all(keys), f"Empty subtype: {path}"
    return dict(zip(keys, rows))


def present(value):
    return str(value).strip().upper() not in {"", "NA", "NAN", "FALSE", "0"}


def integer(value):
    number = float(value)
    assert math.isfinite(number) and number == int(number) and number >= 0
    return int(number)


def selected(mobile, total, breadth, mobility_tenths, host_threshold):
    assert total > 0
    return 10 * mobile >= mobility_tenths * total and breadth >= host_threshold


def save_csv(root, name, rows):
    with (root / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Tier I threshold sensitivity; Supplementary Data S5-2 and S5-3.")
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--assessed", type=Path, required=True)
    parser.add_argument("--original-tier1", type=Path, required=True)
    parser.add_argument("--hosts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if not __debug__:
        raise RuntimeError("Run Python without -O; input checks must remain enabled.")
    inputs = {"core": args.core, "assessed": args.assessed,
              "original_tier1": args.original_tier1}
    root = args.output
    root.mkdir(parents=True, exist_ok=True)
    hosts = {s.strip() for s in args.hosts.read_text(encoding="utf-8-sig").splitlines() if s.strip()}
    if not hosts:
        raise ValueError("The host eligibility file is empty.")
    setting_columns = [f"Site.{s}" for s in SETTINGS]
    core = read_unique(args.core, ["arg_subtype", "core_flag", *setting_columns])
    required = ["arg_subtype", "arg_type", "risk_level", "species_list", "host_breadth",
                "chromosome_contigs", "plasmid_contigs", "virus_contigs", "mobile_contigs",
                "total_contigs", "mobility_ratio", *setting_columns]
    assessed = read_unique(args.assessed, required)
    original = read_unique(args.original_tier1, ["arg_subtype"])
    assert set(assessed) <= set(core) and set(original) <= set(assessed)
    assert all(present(row["core_flag"]) for row in core.values())
    levels = Counter(row["risk_level"] for row in assessed.values())
    scenario_definitions = [
        {"scenario": f"M{m * 10}_H{h}", "mobility_tenths": m,
         "mobility_threshold": m / 10, "host_threshold": h}
        for m in [6, 7, 8] for h in [2, 3, 4]
    ]
    detail = []
    for name, core_row in core.items():
        flags = [int(present(core_row[f"Site.{s}"])) for s in SETTINGS]
        assert any(flags)
        row = {
            "arg_subtype": name,
            "arg_type": name.split("__", 1)[0],
            "original_assessed_set": name in assessed,
            "original_tier": None,
            "original_TierI": int(name in original),
            "hospital_core": flags[0], "wwtp_core": flags[1],
            "community_core": flags[2], "wet_market_core": flags[3],
            "nonhospital_additional_core": int(not flags[0] and any(flags[1:])),
            "host_breadth": None, "chromosome_contigs": None,
            "plasmid_contigs": None, "phage_contigs": None,
            "mobile_contigs": None, "classified_contigs": None,
            "mobility_fraction": None, "species_list": None,
        }
        if name in assessed:
            source = assessed[name]
            assert [int(present(source[f"Site.{s}"])) for s in SETTINGS] == flags
            host_list = [x.strip() for x in source["species_list"].split(";") if x.strip()]
            breadth = integer(source["host_breadth"])
            assert len(host_list) == len(set(host_list)) == breadth
            assert set(host_list) <= hosts and 1 <= breadth <= len(hosts)
            chrom = integer(source["chromosome_contigs"])
            plasmid = integer(source["plasmid_contigs"])
            phage = integer(source["virus_contigs"])
            mobile = integer(source["mobile_contigs"])
            total = integer(source["total_contigs"])
            assert total > 0 and mobile == plasmid + phage
            assert total == chrom + plasmid + phage
            ratio = mobile / total
            assert abs(ratio - float(source["mobility_ratio"])) < 1e-12
            assert selected(mobile, total, breadth, 7, 3) == (name in original)
            assert (source["risk_level"] == "Level I") == (name in original)
            if name in original:
                for key in set(original[name]) & set(source):
                    assert original[name][key] == source[key], f"Roster/trait conflict: {name}, {key}"
            row.update(
                arg_type=source["arg_type"], original_tier=source["risk_level"],
                host_breadth=breadth, chromosome_contigs=chrom,
                plasmid_contigs=plasmid, phage_contigs=phage,
                mobile_contigs=mobile, classified_contigs=total,
                mobility_fraction=ratio, species_list=source["species_list"],
            )
            for definition in scenario_definitions:
                row[definition["scenario"]] = int(selected(
                    mobile, total, breadth,
                    definition["mobility_tenths"], definition["host_threshold"],
                ))
            row["selected_scenario_count"] = sum(row[d["scenario"]] for d in scenario_definitions)
        else:
            for definition in scenario_definitions:
                row[definition["scenario"]] = None
            row["selected_scenario_count"] = None
        detail.append(row)
    detail.sort(key=lambda row: (
        not row["original_TierI"], not row["original_assessed_set"], row["arg_subtype"],
    ))
    nonhospital = {row["arg_subtype"] for row in detail if row["nonhospital_additional_core"]}
    original_set = set(original)
    original_nonhospital = nonhospital & original_set
    scenario_sets = {
        d["scenario"]: {r["arg_subtype"] for r in detail if r[d["scenario"]] == 1}
        for d in scenario_definitions
    }
    summaries = []
    for definition in scenario_definitions:
        selected_set = scenario_sets[definition["scenario"]]
        retained = len(selected_set & original_set)
        summary = dict(definition)
        summary.pop("mobility_tenths")
        summary.update(
            selected_TierI_count=len(selected_set), original_102_retained=retained,
            original_102_retention=retained / len(original_set),
            original_102_not_selected=len(original_set - selected_set),
            alternative_only_count=len(selected_set - original_set),
            jaccard_with_original_102=retained / len(selected_set | original_set),
            nonhospital_additional_selected=len(selected_set & nonhospital),
            original_nonhospital_10_retained=len(selected_set & original_nonhospital),
        )
        summaries.append(summary)
    for a in scenario_definitions:
        for b in scenario_definitions:
            if a["mobility_tenths"] <= b["mobility_tenths"] and a["host_threshold"] <= b["host_threshold"]:
                assert scenario_sets[b["scenario"]] <= scenario_sets[a["scenario"]]
    assert selected(7, 10, 3, 7, 3)
    assert not selected(699, 1000, 3, 7, 3)
    assert not selected(7, 10, 2, 7, 3)
    assert scenario_sets["M70_H3"] == original_set
    common = set.intersection(*scenario_sets.values())
    checks = {
        "core_subtypes": len(core), "original_assessed_subtypes": len(assessed),
        "not_in_original_assessed_set": len(core) - len(assessed),
        "original_tiers": dict(levels), "baseline_exact_match": True,
        "nonhospital_additional_core": len(nonhospital),
        "nonhospital_additional_assessed": len(nonhospital & set(assessed)),
        "original_nonhospital_TierI": len(original_nonhospital),
        "selected_in_all_nine": len(common),
        "original_102_selected_in_all_nine": len(common & original_set),
        "original_nonhospital_10_selected_in_all_nine": len(common & original_nonhospital),
        "monotonicity_and_threshold_boundary_checks": "passed",
        "inputs": {
            key: {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for key, path in {**inputs, "hosts": args.hosts}.items()
        },
    }
    save_csv(root, "TierI_threshold_summary.csv", summaries)
    save_csv(root, "TierI_subtype_membership.csv", detail)
    (root / "sensitivity_results.json").write_text(
        json.dumps({"summary": summaries, "detail": detail, "checks": checks}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps({"summary": summaries, "checks": checks}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
