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
UnresolvedSpan = tuple[str, int, int]  # chromosome, zero-based start, exclusive end


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


def read_unresolved_spans(
    path: Path, child_id: str, distance_bp: int = DISTANCE_BP,
) -> dict[str, list[UnresolvedSpan]]:
    """Read the screened TSV and return excluded spans for each parent.

    A distance cluster is a connected component of brackets across both parents.
    It excludes the entire interval from the first to last marker, including
    gaps between calls, for every parent with a switch in that component.
    """
    if distance_bp < 0:
        raise ValueError("Crossover distance must be non-negative")
    brackets: dict[str, list[tuple[int, int, str, int]]] = defaultdict(list)
    with path.open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if rows.fieldnames != COLUMNS:
            raise ValueError(f"{path}: expected the screened WhatsHap crossover TSV header")
        for row in rows:
            try:
                left, right = int(row["position1"]), int(row["position2"])
                clustered = int(row["clustered"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{rows.line_num}: invalid crossover coordinates or flag") from exc
            chrom, parent = row["chromosome"], row["parent"]
            if not (row["#child_id"] == child_id and chrom and parent in ("father", "mother")
                    and 0 < left < right and clustered in (0, 1)):
                raise ValueError(f"{path}:{rows.line_num}: invalid crossover for child {child_id}")
            brackets[chrom].append((left, right, parent, clustered))

    spans: dict[str, list[UnresolvedSpan]] = {"father": [], "mother": []}
    for chrom, calls in brackets.items():
        ordered = sorted(calls)
        first = 0
        while first < len(ordered):
            left, right, _, _ = ordered[first]
            stop = first + 1
            while stop < len(ordered) and ordered[stop][0] - right <= distance_bp:
                right = max(right, ordered[stop][1])
                stop += 1
            component = ordered[first:stop]
            if any(flag != int(len(component) > 1) for _, _, _, flag in component):
                raise ValueError(f"{path}: clustered flags disagree with the {distance_bp} bp cutoff "
                                 f"at {chrom}:{left}-{right}")
            for parent in {parent for _, _, parent, _ in component}:
                spans[parent].append((chrom, left - 1, right))
            first = stop
    return spans


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
