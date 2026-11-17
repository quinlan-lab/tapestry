#!/usr/bin/env python3
"""Screen crossover marker pairs using indexed long-read BAMs and pedigree phase.

Thresholds are provisional. RETAIN_SUPPORTED is not biological truth, and
insufficient spanning evidence is not a false call. HP tags are never evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
from dataclasses import asdict, dataclass
import hashlib
import json
import logging
import math
from pathlib import Path
from typing import Any

import polars as pl
import pysam

from util.phase_sets import read_ps

LOG = logging.getLogger(__name__)
REQUIRED = ['event_id', '#child_id', 'parent', 'chromosome', 'position1', 'position2']
CALL_COLUMNS = [*REQUIRED, 'parent_id', 'status', 'reasons', 'max_anchor_conflict_fraction',
                'parent_haplotypes', 'parent_haplotype_read_counts', 'parent_phase_concordance',
                'child_transmitted_pattern', 'child_recombinant_reads']
SITE_COLUMNS = ['event_id', '#child_id', 'parent', 'chromosome', 'sample', 'position',
                'unique_reads', 'conflicting_reads', 'conflict_fraction']
OUTPUTS = ['calls.tsv', 'retained.tsv', 'flagged.tsv', 'unresolved.tsv',
           'anchor-evidence.tsv', 'filter.qc.json']


@dataclass(frozen=True)
class Settings:
    mapq: int = 20
    baseq: int = 20
    min_anchor_unique_reads: int = 10
    min_conflicting_reads: int = 3
    min_conflict_fraction: float = 0.20
    min_parent_haplotype_reads: int = 3
    min_parent_concordance: float = 0.90
    min_child_recombinant_reads: int = 3
    child_paternal_haplotype: int = 1

    def validate(self) -> None:
        for name in ('mapq', 'baseq'):
            if not 0 <= getattr(self, name) <= 255:
                raise ValueError(f'{name} must be between 0 and 255')
        for name in ('min_anchor_unique_reads', 'min_conflicting_reads',
                     'min_parent_haplotype_reads', 'min_child_recombinant_reads'):
            if getattr(self, name) < 1:
                raise ValueError(f'{name} must be positive')
        for name in ('min_conflict_fraction', 'min_parent_concordance'):
            value = getattr(self, name)
            if not math.isfinite(value) or not 0 < value <= 1:
                raise ValueError(f'{name} must be finite and in (0, 1]')
        if self.child_paternal_haplotype not in (1, 2):
            raise ValueError('child_paternal_haplotype must be 1 or 2')


def families(path: Path) -> dict[str, tuple[str, str]]:
    result = {}
    seen = set()
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        if len(fields) != 6:
            raise ValueError('PED must have six columns')
        _, child, father, mother, _, _ = fields
        if child in seen:
            raise ValueError(f'Duplicate PED individual: {child}')
        seen.add(child)
        if father != '0' and mother != '0':
            if len({child, father, mother}) != 3:
                raise ValueError(f'Invalid trio for {child}')
            result[child] = (father, mother)
    if not result:
        raise ValueError('PED contains no complete trios')
    return result


def read_calls(path: Path, pedigree: dict[str, tuple[str, str]]) -> list[dict[str, Any]]:
    frame = pl.read_csv(path, separator='\t', infer_schema=False)
    if not set(REQUIRED).issubset(frame.columns):
        raise ValueError(f'Crossover TSV requires columns: {REQUIRED}')
    rows = frame.to_dicts()
    seen = set()
    for row in rows:
        if any(row[k] is None or not str(row[k]).strip() for k in REQUIRED):
            raise ValueError('Crossover TSV contains missing required values')
        if row['event_id'] in seen:
            raise ValueError(f'Duplicate event_id: {row["event_id"]}')
        seen.add(row['event_id'])
        if row['#child_id'] not in pedigree or row['parent'] not in ('father', 'mother'):
            raise ValueError(f'Event {row["event_id"]}: unknown child or parent role')
        row['position1'], row['position2'] = int(row['position1']), int(row['position2'])
        if not 0 < row['position1'] < row['position2']:
            raise ValueError('Marker positions must be positive, one-based and increasing')
        parent = pedigree[row['#child_id']][0 if row['parent'] == 'father' else 1]
        if row.get('parent_id') and row['parent_id'] != parent:
            raise ValueError(f'Event {row["event_id"]}: parent_id disagrees with PED')
        row['parent_id'] = parent
    return rows


def evidence(bam: pysam.AlignmentFile, sample: str, chrom: str, left: int, right: int,
             settings: Settings) -> tuple[list[dict[str, Any]], dict[str, int]]:
    observations: dict[int, dict[str, set[str]]] = {p: defaultdict(set) for p in (left, right)}
    primary_pairs: dict[str, set[str]] = defaultdict(set)
    # Fetch two marker windows only, even when the crossover bracket is huge.
    for index, pos in enumerate((left, right)):
        for read in bam.fetch(chrom, pos - 1, pos):
            if (read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail
                    or read.mapping_quality < settings.mapq):
                continue
            end = read.reference_end
            if end is None:
                continue
            name = read.query_name
            if name is None:
                raise ValueError(f'BAM for {sample} contains an alignment without a read name')
            if index and read.reference_start <= left - 1 < end:
                continue  # Already processed from the left marker window.
            sequence, qualities = read.query_sequence, read.query_qualities
            if sequence is None or qualities is None:
                continue
            bases = {r + 1: sequence[q] for q, r in read.get_aligned_pairs(matches_only=True)
                     if r is not None and q is not None and r + 1 in observations
                     and qualities[q] >= settings.baseq and sequence[q] in 'ACGT'}
            for p, base in bases.items():
                observations[p][name].add(base)
            if not read.is_supplementary and left in bases and right in bases:
                primary_pairs[name].add(bases[left] + '-' + bases[right])
    bad = {name for reads in observations.values() for name, bases in reads.items() if len(bases) > 1}
    sites = []
    for pos, reads in observations.items():
        count = len(reads)
        conflicts = sum(len(bases) > 1 for bases in reads.values())
        sites.append({'sample': sample, 'position': pos, 'unique_reads': count,
                      'conflicting_reads': conflicts,
                      'conflict_fraction': conflicts / count if count else None})
    pairs = Counter(next(iter(bases)) for name, bases in primary_pairs.items()
                    if len(bases) == 1 and name not in bad)
    return sites, dict(pairs)


def alignment_status(sites: list[dict[str, Any]], settings: Settings) -> str:
    if any(s['unique_reads'] >= settings.min_anchor_unique_reads
           and s['conflicting_reads'] >= settings.min_conflicting_reads
           and s['conflict_fraction'] >= settings.min_conflict_fraction for s in sites):
        return 'FLAG_ALIGNMENT_AMBIGUITY'
    if any(s['unique_reads'] < settings.min_anchor_unique_reads for s in sites):
        return 'REVIEW_LOW_COVERAGE'
    return 'CLEAN_ANCHORS'


def predicted(vcf: pysam.VariantFile, sample: str, chrom: str,
              left: int, right: int) -> list[str] | None:
    alleles = []
    phase_sets = set()
    for pos in (left, right):
        matches = [r for r in vcf.fetch(chrom, pos - 1, pos)
                   if r.pos == pos and r.ref is not None and len(r.ref) == 1
                   and r.alts and len(r.alts) == 1 and len(r.alts[0]) == 1]
        if len(matches) != 1:
            return None
        record = matches[0]
        call = record.samples[sample]
        gt = call.get('GT')
        if gt is None or len(gt) != 2 or any(a not in (0, 1) for a in gt):
            return None
        if gt[0] != gt[1]:
            ps = read_ps([call.get('PS')], 0)
            if not call.phased or ps is None:
                return None
            phase_sets.add(ps)
        values = record.alleles
        if values is None:
            return None
        alleles.append([values[a] for a in gt])
    if len(phase_sets) > 1:
        return None
    return ['-'.join(pair[h] for pair in alleles) for h in (0, 1)]


def classify(call: dict[str, Any], sites: list[dict[str, Any]], pairs: dict[str, dict[str, int]],
             parent_haps: list[str] | None, child_haps: list[str] | None,
             settings: Settings) -> dict[str, Any]:
    parent, child = call['parent_id'], call['#child_id']
    paternal_index = settings.child_paternal_haplotype - 1
    index = paternal_index if call['parent'] == 'father' else 1 - paternal_index
    expected_child = child_haps[index] if child_haps else None
    support = [pairs[parent].get(h, 0) for h in parent_haps] if parent_haps else None
    total = sum(pairs[parent].values())
    concordance = (sum(pairs[parent].get(h, 0) for h in set(parent_haps)) / total
                   if parent_haps and total else None)
    child_support = pairs[child].get(expected_child, 0) if expected_child else None
    reasons = []
    alignment = alignment_status(sites, settings)
    if alignment != 'CLEAN_ANCHORS':
        reasons.append(alignment)
    phase_resolved = call.get('phase_status', 'not_requested') in ('resolved', 'not_requested')
    if not phase_resolved or not parent_haps or len(set(parent_haps)) != 2 or not expected_child:
        reasons.append('REVIEW_UNRESOLVED_PHASE')
    elif any(a == b for a, b in zip(parent_haps[0].split('-'), parent_haps[1].split('-'))):
        reasons.append('REVIEW_UNINFORMATIVE_PARENT_MARKERS')
    elif expected_child in parent_haps:
        reasons.append('REVIEW_NO_RECOMBINANT_PATTERN')
    if support is None or min(support) < settings.min_parent_haplotype_reads:
        reasons.append('REVIEW_INSUFFICIENT_PARENT_SPANNING_READS')
    if concordance is not None and concordance < settings.min_parent_concordance:
        reasons.append('FLAG_PARENT_PHASE_DISAGREEMENT')
    if child_support is None or child_support < settings.min_child_recombinant_reads:
        reasons.append('REVIEW_INSUFFICIENT_CHILD_SPANNING_READS')
    status = ('FLAG_READ_CONFLICT' if any(r.startswith('FLAG_') for r in reasons)
              else 'REVIEW_INSUFFICIENT_EVIDENCE' if reasons else 'RETAIN_SUPPORTED')
    return {**{k: call[k] for k in REQUIRED}, 'parent_id': parent, 'status': status,
            'reasons': ';'.join(reasons),
            'max_anchor_conflict_fraction': max((s['conflict_fraction'] or 0 for s in sites), default=0),
            'parent_haplotypes': parent_haps, 'parent_haplotype_read_counts': support,
            'parent_phase_concordance': concordance, 'child_transmitted_pattern': expected_child,
            'child_recombinant_reads': child_support, 'clean_primary_pair_counts': pairs, 'anchors': sites}


def write_table(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    values = [{k: json.dumps(r.get(k)) if isinstance(r.get(k), (list, dict)) else r.get(k)
               for k in columns} for r in rows]
    frame = pl.DataFrame(values, strict=False) if values else pl.DataFrame(schema={k: pl.String for k in columns})
    frame.select(columns).write_csv(path, separator='\t')


def sha256(path: Path) -> str:
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def run(crossovers: Path, phased_vcf: Path, ped: Path, bams: dict[str, Path],
        output_dir: Path, settings: Settings = Settings()) -> dict[str, Any]:
    settings.validate()
    pedigree = families(ped)
    calls = read_calls(crossovers, pedigree)
    samples = sorted({s for c in calls for s in (c['#child_id'], *pedigree[c['#child_id']])})
    if not set(samples).issubset(bams):
        raise ValueError(f'Missing BAM mappings: {sorted(set(samples) - set(bams))}')
    inputs = [crossovers, phased_vcf, ped, *(bams[s] for s in samples)]
    if {p.resolve() for p in inputs} & {(output_dir / name).resolve() for name in OUTPUTS}:
        raise ValueError('Output path collides with an input file')
    LOG.info('Inputs: crossovers=%s VCF=%s PED=%s BAMs=%s', crossovers, phased_vcf, ped, bams)
    LOG.info('Settings: %s; %d candidate events', asdict(settings), len(calls))
    results, site_rows = [], []
    with ExitStack() as stack:
        vcf = stack.enter_context(pysam.VariantFile(str(phased_vcf)))
        if vcf.index is None:
            raise ValueError('Phased VCF must have a readable index')
        if not set(samples).issubset(vcf.header.samples):
            raise ValueError('Phased VCF is missing a selected pedigree member')
        readers = {}
        for sample in samples:
            bam = stack.enter_context(pysam.AlignmentFile(str(bams[sample]), 'rb'))
            if not bam.has_index():
                raise ValueError(f'BAM for {sample} must have a readable index')
            names = {r['SM'] for r in bam.header.to_dict().get('RG', []) if 'SM' in r}
            if names and names != {sample}:
                raise ValueError(f'BAM read-group samples {names} disagree with mapping {sample}')
            readers[sample] = bam
        for call in calls:
            chrom, left, right = call['chromosome'], call['position1'], call['position2']
            child = call['#child_id']
            trio = (child, *pedigree[child])
            if chrom not in vcf.header.contigs or any(chrom not in readers[s].references for s in trio):
                raise ValueError(f'Event {call["event_id"]}: contig {chrom} missing from VCF or BAM')
            if any(right > readers[s].get_reference_length(chrom) for s in trio):
                raise ValueError(f'Event {call["event_id"]}: marker exceeds BAM contig length')
            sites, pairs = [], {}
            for sample in trio:
                ss, pp = evidence(readers[sample], sample, chrom, left, right, settings)
                sites.extend(ss)
                pairs[sample] = pp
                site_rows.extend({**{k: call[k] for k in ('event_id', '#child_id', 'parent', 'chromosome')},
                                  **s} for s in ss)
            results.append(classify(call, sites, pairs,
                                    predicted(vcf, call['parent_id'], chrom, left, right),
                                    predicted(vcf, child, chrom, left, right), settings))
    report = {'settings': asdict(settings), 'vcf_sha256': sha256(phased_vcf), 'sample_ids': samples,
              'inputs_sha256': {str(p): sha256(p) for p in (crossovers, phased_vcf, ped)},
              'bam_metadata': {s: {'path': str(bams[s]), 'bytes': bams[s].stat().st_size,
                                    'mtime_ns': bams[s].stat().st_mtime_ns} for s in samples},
              'pysam_version': pysam.__version__, 'calls': results,
              'status_counts': dict(Counter(r['status'] for r in results)),
              'coordinates': 'Input and output marker positions are one-based inclusive.',
              'limitations': [
                  'Provisional thresholds developed on one chr22 family fixture, not independent calibration.',
                  'Assumes independent long-read molecules have unique read names within each sample BAM.',
                  'Child haplotypes must be oriented to parental origin as specified; HP tags are not support.',
                  'Both markers must lie in a single primary alignment for positive spanning support.',
                  'Insufficient or unresolved evidence is not a false call.',
                  'No document overlap, pedigree site whitelist, event spacing or coordinate mask is used.']}
    output_dir.mkdir(parents=True, exist_ok=True)
    write_table(output_dir / 'calls.tsv', results, CALL_COLUMNS)
    for name, status in [('retained', 'RETAIN_SUPPORTED'), ('flagged', 'FLAG_READ_CONFLICT'),
                         ('unresolved', 'REVIEW_INSUFFICIENT_EVIDENCE')]:
        write_table(output_dir / f'{name}.tsv', [r for r in results if r['status'] == status], CALL_COLUMNS)
    write_table(output_dir / 'anchor-evidence.tsv', site_rows, SITE_COLUMNS)
    (output_dir / 'filter.qc.json').write_text(json.dumps(report, indent=2) + '\n')
    LOG.info('Completed %d events: %s; outputs: %s', len(results), report['status_counts'], output_dir)
    return report


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--crossovers', type=Path, required=True, help='Parent-specific crossover TSV from report_whatshap_crossovers.py')
    p.add_argument('--phased-vcf', type=Path, required=True, help='Indexed VCF with oriented child haplotypes')
    p.add_argument('--ped', type=Path, required=True, help='Six-column pedigree; may contain multiple families')
    p.add_argument('--bam', nargs=2, action='append', required=True, metavar=('SAMPLE', 'PATH'), help='Repeat for each sample used by the crossover TSV')
    p.add_argument('--output-dir', type=Path, required=True, help='Directory for evidence tables, classified subsets and QC')
    for flag, dest, default in [('--min-mapq', 'mapq', 20), ('--min-baseq', 'baseq', 20),
                                ('--min-anchor-reads', 'min_anchor_unique_reads', 10),
                                ('--min-conflicting-reads', 'min_conflicting_reads', 3),
                                ('--min-parent-haplotype-reads', 'min_parent_haplotype_reads', 3),
                                ('--min-child-recombinant-reads', 'min_child_recombinant_reads', 3)]:
        help_text = {
            'mapq': 'Minimum alignment mapping quality',
            'baseq': 'Minimum base quality at each marker',
            'min_anchor_unique_reads': 'Minimum unique molecules per marker in each trio member',
            'min_conflicting_reads': 'Minimum within-molecule conflicts to flag a marker',
            'min_parent_haplotype_reads': 'Minimum clean spanning molecules per parental haplotype',
            'min_child_recombinant_reads': 'Minimum child molecules supporting the recombinant pair',
        }[dest]
        p.add_argument(flag, dest=dest, type=int, default=default, help=help_text)
    p.add_argument('--min-conflict-fraction', type=float, default=0.20,
                   help='Minimum conflicting-molecule fraction to flag an adequately covered marker')
    p.add_argument('--min-parent-concordance', type=float, default=0.90,
                   help='Minimum fraction of clean parent pairs matching the predicted haplotypes')
    p.add_argument('--child-paternal-haplotype', type=int, choices=(1, 2), default=1,
                   help='One-based GT allele index of the paternal child haplotype; the other is maternal')
    return p


def main() -> None:
    p = parser()
    args = p.parse_args()
    bams = {}
    for sample, path in args.bam:
        if sample in bams:
            p.error(f'Duplicate BAM mapping: {sample}')
        bams[sample] = Path(path)
    settings = Settings(**{k: getattr(args, k) for k in asdict(Settings())})
    try:
        run(args.crossovers, args.phased_vcf, args.ped, bams, args.output_dir, settings)
    except (ValueError, OSError, KeyError, pl.exceptions.PolarsError) as exc:
        p.exit(2, f'Error: {exc}\n')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    main()
