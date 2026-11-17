"""Trio filtering must retain homozygotes and test every selected sample."""
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from prefilter_whatshap_vcf import prefilter


class PrefilterTests(unittest.TestCase):
    def test_trio_quality_missing_calls_filter_status_and_homozygotes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / 'input.vcf'
            header = ('##fileformat=VCFv4.2\n##contig=<ID=chr22,length=1000>\n'
                      '##FILTER=<ID=FAIL,Description="Fail">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
                      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkid\tdad\tmom\tother\n')
            good = ['0/1:20:10', '0/0:30:40', '1/1:99:190']
            cases = [
                (20, '.', good), (20, 'PASS', good),
                (19, '.', good), (20, 'FAIL', good),
                (20, '.', ['0/1:19:10', *good[1:]]),
                (20, '.', [good[0], '0/0:30:9', good[2]]),
                (20, '.', [*good[:2], './.:99:30']),
                (20, '.', [*good[:2], '1/1:.:30']),
                (20, '.', [*good[:2], '1/1:99:.']),
            ]
            source.write_text(header + ''.join(
                f'chr22\t{i}\t.\tA\tG\t{qual}\t{status}\t.\tGT:GQ:DP\t'
                + '\t'.join([*calls, './.:.:.']) + '\n'
                for i, (qual, status, calls) in enumerate(cases, 1)))
            target = root / 'filtered.vcf.gz'
            report = prefilter(source, target, ['kid', 'dad', 'mom'])
            self.assertEqual(report['retained_records'], 2)
            with pysam.VariantFile(str(target)) as reader:
                self.assertEqual(list(reader.header.samples), ['kid', 'dad', 'mom'])
                rows = list(reader.fetch('chr22'))
                self.assertEqual([r.pos for r in rows], [1, 2])
                self.assertEqual(rows[0].samples['dad']['GT'], (0, 0))
                self.assertEqual(rows[0].samples['mom']['GT'], (1, 1))
            with self.assertRaisesRegex(ValueError, 'distinct samples'):
                prefilter(source, target, ['kid', 'kid', 'mom'])


if __name__ == '__main__':
    unittest.main()
