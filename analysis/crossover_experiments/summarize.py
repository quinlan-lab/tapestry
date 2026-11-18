"""Compare parent-specific WhatsHap 2.8 switches with nearby calls and BED masks."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from report_whatshap_crossovers import COLUMNS, DISTANCE_BP, nearest_gaps, switches


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
    segdups = bed_intervals(args.segdups)
    blacklist = bed_intervals(args.blacklist)
    columns = [*COLUMNS, "segdup_overlap", "blacklist_overlap"]
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
