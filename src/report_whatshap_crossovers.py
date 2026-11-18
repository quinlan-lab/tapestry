"""Flag WhatsHap crossover candidates whose marker brackets are close together."""

import argparse
import csv
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import TypedDict

HEADER = (
    "#child_id chromosome position1 position2 transmitted_hap_father1 "
    "transmitted_hap_father2 transmitted_hap_mother1 transmitted_hap_mother2 recombination_cost"
).split()
DISTANCE_BP = 10_000
COLUMNS = ["event_id", "#child_id", "chromosome", "parent", "position1", "position2",
           "state_before", "state_after", "recombination_cost", "nearest_gap_bp", "clustered"]
Switch = TypedDict("Switch", {
    "event_id": str, "#child_id": str, "chromosome": str, "parent": str,
    "position1": int, "position2": int, "state_before": str, "state_after": str,
    "recombination_cost": int,
})


def switches(path: Path) -> Iterator[Switch]:
    with path.open() as handle:
        if handle.readline().split() != HEADER:
            raise ValueError(f"{path}: expected the WhatsHap 2.8 recombination header")
        for number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            try:
                child, chrom, left_text, right_text, father1, father2, mother1, mother2, cost_text = line.split()
                left, right, cost = int(left_text), int(right_text), int(cost_text)
            except ValueError as exc:
                raise ValueError(f"{path}:{number}: invalid recombination record: {line.strip()}") from exc
            states = (father1, father2, mother1, mother2)
            if not (0 < left < right and cost >= 0 and set(states) <= {"0", "1"}):
                raise ValueError(f"{path}:{number}: invalid recombination record: {line.strip()}")
            if father1 == father2 and mother1 == mother2:
                raise ValueError(f"{path}:{number}: no parental switch: {line.strip()}")
            for parent, before, after in (("father", father1, father2), ("mother", mother1, mother2)):
                if before != after:
                    yield {
                        "event_id": f"{child}:{chrom}:{left}-{right}:{parent}",
                        "#child_id": child, "chromosome": chrom, "parent": parent,
                        "position1": left, "position2": right, "state_before": before,
                        "state_after": after, "recombination_cost": cost,
                    }


def nearest_gaps(calls: list[Switch]) -> dict[str, int | None]:
    """Measure bracket edge distances across both parents for each child and contig."""
    groups: dict[tuple[str, str], list[Switch]] = defaultdict(list)
    for call in calls:
        groups[call["#child_id"], call["chromosome"]].append(call)
    gaps: dict[str, int | None] = {}
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
    if len(calls) != len(gaps):
        raise ValueError("Duplicate crossover bracket for the same child and parent")
    return gaps


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="WhatsHap recombination report")
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--distance-bp", type=int, default=DISTANCE_BP,
                        help="Maximum nearest bracket gap to flag (default: %(default)s)")
    parser.add_argument("--tsv", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    args = parser.parse_args()
    if args.distance_bp < 0:
        parser.error("--distance-bp must be a non-negative integer")

    print(f"Screening {args.input} on {args.chromosome} at gaps <= {args.distance_bp} bp", flush=True)
    calls = sorted(switches(args.input), key=lambda c: (c["position1"], c["position2"], c["event_id"]))
    if any(call["chromosome"] != args.chromosome for call in calls):
        raise ValueError(f"{args.input}: report contains a crossover outside {args.chromosome}")
    gaps = nearest_gaps(calls)
    flagged = 0
    with args.tsv.open("w") as table, args.bed.open("w") as bed:
        writer = csv.DictWriter(table, COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for call in calls:
            gap = gaps[call["event_id"]]
            clustered = int(gap is not None and gap <= args.distance_bp)
            flagged += clustered
            writer.writerow({**call, "nearest_gap_bp": gap, "clustered": clustered})
            start, end = call["position1"] - 1, call["position2"]
            print(call["chromosome"], start, end, call["event_id"], 0, ".", start, end, 0,
                  gap if gap is not None else ".", clustered, sep="\t", file=bed)
    print(f"Wrote {len(calls)} candidates ({flagged} flagged) to {args.tsv} and {args.bed}")


if __name__ == "__main__":
    main()
