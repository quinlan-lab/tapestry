"""Compare parent-specific WhatsHap 2.8 switches with nearby calls and BED masks."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

HEADER = (
    "#child_id chromosome position1 position2 transmitted_hap_father1 "
    "transmitted_hap_father2 transmitted_hap_mother1 transmitted_hap_mother2 recombination_cost"
).split()
DISTANCE_BP = 10_000


def switches(path):
    with path.open() as handle:
        if handle.readline().split() != HEADER:
            raise ValueError(f"{path}: expected the WhatsHap 2.8 recombination header")
        for line in handle:
            if not line.strip():
                continue
            child, chrom, left, right, father1, father2, mother1, mother2, cost = line.split()
            left, right, cost = int(left), int(right), int(cost)
            states = (father1, father2, mother1, mother2)
            if not (0 < left < right and cost >= 0 and set(states) <= {"0", "1"}):
                raise ValueError(f"{path}: invalid recombination record: {line.strip()}")
            if father1 == father2 and mother1 == mother2:
                raise ValueError(f"{path}: no parental switch: {line.strip()}")
            for parent, before, after in (("father", father1, father2), ("mother", mother1, mother2)):
                if before != after:
                    yield {
                        "event_id": f"{child}:{chrom}:{left}-{right}:{parent}",
                        "#child_id": child, "chromosome": chrom, "parent": parent,
                        "position1": left, "position2": right, "state_before": before,
                        "state_after": after, "recombination_cost": cost,
                    }


def nearest_gaps(calls):
    """Bracket edge distances, pooling both parental roles within each child/contig."""
    groups = defaultdict(list)
    for call in calls:
        groups[call["#child_id"], call["chromosome"]].append(call)
    gaps = {}
    for group in groups.values():
        ordered = sorted(group, key=lambda c: c["position1"])
        max_end = 0
        for i, call in enumerate(ordered):
            neighbors = []
            if i > 0:
                neighbors.append(max(0, call["position1"] - max_end))
            if i + 1 < len(ordered):
                neighbors.append(max(0, ordered[i + 1]["position1"] - call["position2"]))
            gaps[call["event_id"]] = min(neighbors, default=None)
            max_end = max(max_end, call["position2"])
    return gaps


def bed_intervals(path):
    intervals = defaultdict(set)
    with path.open() as handle:
        for line in handle:
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            if not 0 <= start < end:
                raise ValueError(f"{path}: invalid BED interval: {line.strip()}")
            intervals[chrom].add((start, end))
    return intervals


def overlaps(call, intervals):
    """Intersect one-based inclusive marker brackets with zero-based BED intervals."""
    left, right = call["position1"] - 1, call["position2"]
    return any(start < right and end > left for start, end in intervals.get(call["chromosome"], ()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("segdups", type=Path)
    parser.add_argument("blacklist", type=Path)
    args = parser.parse_args()
    calls = list(switches(args.run_dir / "recombination.txt"))
    gaps = nearest_gaps(calls)
    if len(calls) != len(gaps):
        raise ValueError("Duplicate crossover bracket for the same child and parent")
    segdups = bed_intervals(args.segdups)
    blacklist = bed_intervals(args.blacklist)
    columns = ["event_id", "#child_id", "chromosome", "parent", "position1", "position2",
               "state_before", "state_after", "recombination_cost", "nearest_gap_bp",
               "clustered", "segdup_overlap", "blacklist_overlap"]
    with (args.run_dir / "crossovers.tsv").open("w") as table, (args.run_dir / "crossovers.bed").open("w") as bed:
        writer = csv.DictWriter(table, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        counts = Counter()
        for call in calls:
            gap = gaps[call["event_id"]]
            clustered = gap is not None and gap <= DISTANCE_BP
            segdup = overlaps(call, segdups)
            blacklisted = overlaps(call, blacklist)
            writer.writerow({**call, "nearest_gap_bp": gap, "clustered": int(clustered),
                             "segdup_overlap": int(segdup), "blacklist_overlap": int(blacklisted)})
            print(call["chromosome"], call["position1"] - 1, call["position2"],
                  call["parent"] + ("_clustered" if clustered else ""), sep="\t", file=bed)
            counts.update({"clustered": int(clustered), "segdup": int(segdup), "blacklist": int(blacklisted)})
    parents = Counter(call["parent"] for call in calls)
    print(f"{args.run_dir}: {len(calls)} switches ({parents['father']} paternal, "
          f"{parents['mother']} maternal), {counts['clustered']} clustered, "
          f"{counts['segdup']} segdup overlaps, {counts['blacklist']} blacklist overlaps")


if __name__ == "__main__":
    main()
