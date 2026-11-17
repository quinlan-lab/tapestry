"""Read-counting, phase and family invariants for crossover evidence screening."""
from dataclasses import replace
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from filter_whatshap_crossovers_by_reads import (
    REQUIRED, OUTPUTS, Settings, alignment_status, classify, evidence, families, predicted, run,
)

CHROM = 'contig_custom_7'
LEFT, RIGHT = 101, 201


def alignment(name, pair='A-C', start=80, length=200, flag=0, mapq=60, baseq=40):
    read = pysam.AlignedSegment()
    read.query_name = name
    bases = list('A' * length)
    for pos, base in zip((LEFT, RIGHT), pair.split('-')):
        if start <= pos - 1 < start + length:
            bases[pos - 1 - start] = base
    read.query_sequence = ''.join(bases)
    read.query_qualities = [baseq] * length
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigar = [(0, length)]
    return read


def write_bam(path, sample, reads):
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': CHROM, 'LN': 1000}], 'RG': [{'ID': sample, 'SM': sample}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as writer:
        for read in sorted(reads, key=lambda r: r.reference_start):
            read.set_tag('RG', sample)
            writer.write(read)
    pysam.index(str(path))


def write_vcf(path, samples, overrides=None):
    overrides = overrides or {}
    header = ('##fileformat=VCFv4.2\n'
              f'##contig=<ID={CHROM},length=1000>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
              '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
    rows = []
    for pos, ref, alt in [(LEFT, 'A', 'G'), (RIGHT, 'C', 'T')]:
        calls = []
        for sample in samples:
            role = sample.split('_')[0]
            gt = ('0|1' if pos == LEFT else '1|0') if role == 'kid' else (
                '0|1' if role == 'dad' else ('1|0' if pos == LEFT else '0|1'))
            # Different samples can legitimately have different PS identifiers.
            ps = {'kid': 11, 'dad': 22, 'mom': 33}[role]
            calls.append(overrides.get((sample, pos), f'{gt}:{ps}'))
        rows.append(f'{CHROM}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT:PS\t' + '\t'.join(calls) + '\n')
    plain = path.with_suffix('')
    plain.write_text(header + ''.join(rows))
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset='vcf', force=True)


def write_calls(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=REQUIRED)
        writer.writeheader()
        writer.writerows(rows)


def fixture(root, nfamilies=1):
    ped, calls, vcf = root / 'pedigree.ped', root / 'events.tsv', root / 'phased.vcf.gz'
    samples, mappings, rows, pedigree = [], {}, [], []
    for i in range(nfamilies):
        child, dad, mom = [f'{role}_{i}' for role in ('kid', 'dad', 'mom')]
        pedigree.append(f'F{i} {child} {dad} {mom} 1 -9\n')
        rows.append(dict(zip(REQUIRED, [f'event{i}', child, 'father', CHROM, LEFT, RIGHT])))
        for sample, haps in [(child, ['A-T', 'G-C']), (dad, ['A-C', 'G-T']), (mom, ['G-C', 'A-T'])]:
            samples.append(sample)
            path = root / f'{sample}.bam'
            reads = [alignment(f'{sample}-{h}-{j}', h) for h in haps for j in range(3)]
            write_bam(path, sample, reads)
            mappings[sample] = path
    ped.write_text(''.join(pedigree))
    write_calls(calls, rows)
    write_vcf(vcf, samples)
    return calls, vcf, ped, mappings


class MoleculeTests(unittest.TestCase):
    def test_unique_molecules_overlap_conflicts_quality_and_no_stitching(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'reads.bam'
            reads = [alignment('ok'), alignment('ok'), alignment('ok', flag=2048),
                     alignment('conflict'), alignment('conflict', 'G-C', flag=2048),
                     alignment('split', 'A-T', length=40),
                     alignment('split', 'A-T', start=180, length=40, flag=2048),
                     alignment('supp-only', 'A-T', flag=2048),
                     alignment('lowmap', mapq=19), alignment('lowbase', baseq=19),
                     alignment('secondary', flag=256), alignment('duplicate', flag=1024),
                     alignment('qcfail', flag=512)]
            write_bam(path, 'dad', reads)
            with pysam.AlignmentFile(str(path)) as bam:
                sites, pairs = evidence(bam, 'dad', CHROM, LEFT, RIGHT, Settings())
            self.assertEqual([s['unique_reads'] for s in sites], [4, 4])
            self.assertEqual([s['conflicting_reads'] for s in sites], [1, 0])
            self.assertEqual(pairs, {'A-C': 1})
            self.assertEqual(alignment_status(sites, replace(Settings(), min_anchor_unique_reads=3,
                min_conflicting_reads=1, min_conflict_fraction=.25)), 'FLAG_ALIGNMENT_AMBIGUITY')
            self.assertEqual(alignment_status(sites, Settings()), 'REVIEW_LOW_COVERAGE')

    def test_fetches_only_two_marker_windows_for_wide_bracket(self):
        class Reader:
            def __init__(self):
                self.windows = []
            def fetch(self, *window):
                self.windows.append(window)
                return []
        reader = Reader()
        sites, pairs = evidence(reader, 's', CHROM, 1, 100000000, Settings())
        self.assertEqual(reader.windows, [(CHROM, 0, 1), (CHROM, 99999999, 100000000)])
        self.assertEqual(pairs, {})
        self.assertTrue(all(s['conflict_fraction'] is None for s in sites))
        self.assertEqual(alignment_status(sites, Settings()), 'REVIEW_LOW_COVERAGE')


class PhaseAndDecisionTests(unittest.TestCase):
    def test_actual_sample_phase_sets_and_homozygous_missing_phase(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path = root / 'v.vcf.gz'
            samples = ['kid_0', 'dad_0', 'mom_0']
            write_vcf(path, samples)
            with pysam.VariantFile(str(path)) as v:
                self.assertEqual(predicted(v, 'kid_0', CHROM, LEFT, RIGHT), ['A-T', 'G-C'])
                self.assertEqual(predicted(v, 'dad_0', CHROM, LEFT, RIGHT), ['A-C', 'G-T'])
            for changed in ['1|0:999', '0/1:.', '1|0:-1', '.|0:11']:
                write_vcf(path, samples, {('kid_0', RIGHT): changed})
                with pysam.VariantFile(str(path)) as v:
                    self.assertIsNone(predicted(v, 'kid_0', CHROM, LEFT, RIGHT))
            write_vcf(path, samples, {('kid_0', LEFT): '0/0:.'})
            with pysam.VariantFile(str(path)) as v:
                self.assertEqual(predicted(v, 'kid_0', CHROM, LEFT, RIGHT), ['A-T', 'A-C'])

    def test_thresholds_orientation_parent_role_and_unresolved_are_distinct(self):
        call = dict(zip(REQUIRED, ['id', 'c', 'father', CHROM, LEFT, RIGHT]))
        call['parent_id'] = 'p'
        sites = [{'unique_reads': 10, 'conflicting_reads': 0, 'conflict_fraction': 0}]
        pairs = {'p': {'A-C': 3, 'G-T': 3}, 'c': {'A-T': 3, 'G-C': 4}}
        settings = Settings()
        def check(**kwargs):
            return classify(kwargs.get('call', call), sites, kwargs.get('pairs', pairs),
                kwargs.get('parent_haps', ['A-C', 'G-T']), kwargs.get('child_haps', ['A-T', 'G-C']),
                kwargs.get('settings', settings))
        self.assertEqual(check()['status'], 'RETAIN_SUPPORTED')
        self.assertEqual(check(settings=replace(settings, min_child_recombinant_reads=4))['status'], 'REVIEW_INSUFFICIENT_EVIDENCE')
        self.assertEqual(check(settings=replace(settings, child_paternal_haplotype=2))['child_transmitted_pattern'], 'G-C')
        self.assertEqual(check(call={**call, 'parent': 'mother'})['child_transmitted_pattern'], 'G-C')
        self.assertEqual(check(call={**call, 'phase_status': 'conflicting_flanks'})['status'], 'REVIEW_INSUFFICIENT_EVIDENCE')
        self.assertEqual(check(parent_haps=None)['status'], 'REVIEW_INSUFFICIENT_EVIDENCE')
        self.assertIn('REVIEW_UNINFORMATIVE_PARENT_MARKERS', check(parent_haps=['A-C', 'A-T'])['reasons'])
        self.assertIn('REVIEW_NO_RECOMBINANT_PATTERN', check(child_haps=['A-C', 'G-C'])['reasons'])
        bad = {'p': {**pairs['p'], 'A-T': 2}, 'c': pairs['c']}
        self.assertIn('FLAG_PARENT_PHASE_DISAGREEMENT', check(pairs=bad)['reasons'])
        self.assertEqual(check(pairs=bad)['status'], 'FLAG_READ_CONFLICT')

    def test_invalid_settings(self):
        for change in [dict(mapq=-1), dict(baseq=256), dict(min_anchor_unique_reads=0),
                       dict(min_conflicting_reads=0), dict(min_conflict_fraction=0),
                       dict(min_conflict_fraction=float('nan')), dict(min_parent_concordance=1.1),
                       dict(min_child_recombinant_reads=0), dict(child_paternal_haplotype=3)]:
            with self.assertRaises(ValueError):
                replace(Settings(), **change).validate()


class InputAndCLITests(unittest.TestCase):
    def test_multifamily_cli_nonstandard_contig_and_parameter_effect(self):
        with tempfile.TemporaryDirectory(prefix='read evidence ') as temp:
            root = Path(temp)
            calls, vcf, ped, mappings = fixture(root, 2)
            output = root / 'results'
            command = [sys.executable, str(ROOT / 'src/filter_whatshap_crossovers_by_reads.py'),
                       '--crossovers', str(calls), '--phased-vcf', str(vcf), '--ped', str(ped),
                       '--output-dir', str(output), '--min-anchor-reads', '3']
            for sample, path in mappings.items():
                command.extend(['--bam', sample, str(path)])
            completed = subprocess.run(command, cwd='/tmp', capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            qc = json.loads((output / 'filter.qc.json').read_text())
            self.assertEqual(qc['status_counts'], {'RETAIN_SUPPORTED': 2})
            self.assertEqual(qc['settings']['min_anchor_unique_reads'], 3)
            self.assertEqual(len(qc['sample_ids']), 6)
            self.assertEqual({r['#child_id'] for r in qc['calls']}, {'kid_0', 'kid_1'})
            for name in OUTPUTS:
                self.assertTrue((output / name).is_file())
            with (output / 'unresolved.tsv').open() as f:
                self.assertEqual(len(list(csv.DictReader(f, delimiter='\t'))), 0)
            qc = run(calls, vcf, ped, mappings, output, replace(Settings(), min_anchor_unique_reads=3,
                     min_parent_haplotype_reads=4))
            self.assertEqual(qc['status_counts'], {'REVIEW_INSUFFICIENT_EVIDENCE': 2})
            duplicate = subprocess.run(command + ['--bam', 'kid_0', str(mappings['kid_0'])],
                                       capture_output=True, text=True)
            self.assertNotEqual(duplicate.returncode, 0)
            self.assertIn('Duplicate BAM mapping', duplicate.stderr)

    def test_empty_report_has_headers_and_no_retained_events(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            calls, vcf, ped, mappings = fixture(root)
            write_calls(calls, [])
            report = run(calls, vcf, ped, {}, root / 'out')
            self.assertEqual(report['calls'], [])
            self.assertEqual(report['status_counts'], {})
            for filename in OUTPUTS[:-1]:
                with (root / 'out' / filename).open() as handle:
                    reader = csv.DictReader(handle, delimiter='\t')
                    self.assertTrue(reader.fieldnames)
                    self.assertEqual(list(reader), [])

    def test_invalid_mapping_ids_phase_inputs_and_collisions(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            calls, vcf, ped, mappings = fixture(root)
            with self.assertRaisesRegex(ValueError, 'Missing BAM'):
                run(calls, vcf, ped, {}, root / 'out')
            with self.assertRaisesRegex(ValueError, 'collides'):
                copy = root / 'calls.tsv'
                copy.write_text(calls.read_text())
                run(copy, vcf, ped, mappings, root)
            original = calls.read_text()
            calls.write_text(original.replace('kid_0', 'unknown'))
            with self.assertRaisesRegex(ValueError, 'unknown child'):
                run(calls, vcf, ped, mappings, root / 'out')
            calls.write_text(original + original.splitlines()[1] + '\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate event_id'):
                run(calls, vcf, ped, mappings, root / 'out')
            calls.write_text(original)
            wrong = root / 'wrong.bam'
            write_bam(wrong, 'wrong_sample', [alignment('a')])
            with self.assertRaisesRegex(ValueError, 'read-group samples'):
                run(calls, vcf, ped, {**mappings, 'dad_0': wrong}, root / 'out')
            Path(str(mappings['kid_0']) + '.bai').unlink()
            with self.assertRaisesRegex(ValueError, 'readable index'):
                run(calls, vcf, ped, mappings, root / 'out')
            with self.assertRaisesRegex(ValueError, 'Phased VCF'):
                run(calls, vcf.with_suffix(''), ped, mappings, root / 'out')
            ped.write_text(ped.read_text() * 2)
            with self.assertRaisesRegex(ValueError, 'Duplicate PED'):
                families(ped)
