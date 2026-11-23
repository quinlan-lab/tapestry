"""Model-only publication invariants, including fully unresolved chromosomes."""
import argparse
import tempfile
import unittest
from pathlib import Path

import polars as pl
import pysam
import pyBigWig

from trio_model_methylation import run


class TrioModelTests(unittest.TestCase):
    def test_model_labels_nulls_original_variants_and_indexes(self):
        for resolved in (True, False):
            with self.subTest(resolved=resolved), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                ref = root / 'ref.fa'
                ref.write_text('>chr22\nCGCGCGCGCG\n')
                pysam.faidx(str(ref))
                model = root / 'model'
                model.mkdir()
                for role in ('kid', 'dad', 'mom'):
                    for hap, value in [('hap1', 20), ('hap2', 80), ('combined', 50)]:
                        bed = model / f'{role}.{hap}.bed'
                        bed.write_text('##pileup-mode=model\n#chrom\tbegin\tend\tmod_score\ttype\tcov\n' +
                                       ''.join(f'chr22\t{i}\t{i+1}\t{value}\tmod\t12\n' for i in (0, 2, 4, 6)))
                        pysam.tabix_compress(str(bed), str(bed) + '.gz')
                        pysam.tabix_index(str(bed) + '.gz', preset='bed')
                for side, hap in [('paternal', 'A'), ('maternal', 'D')]:
                    pl.DataFrame([('chr22', 0, 2, hap, 1., 5)] if resolved else [], schema={
                        'chrom': pl.String, 'start': pl.Int64, 'end': pl.Int64,
                        f'{side}_haplotype': pl.String, 'haplotype_concordance': pl.Float64,
                        'num_het_SNVs_in_parent': pl.Int64}, orient='row').write_csv(
                            root / f'{side}.hap-map.tsv', separator='\t')
                    (root / f'K.bit-vector-sites-mismatches.{side}.bed').touch()
                vcf = root / 'original.vcf'
                # This low-quality variant would have been excluded from WhatsHap.
                vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tK\tF\tM\n'
                               'chr22\t2\t.\tG\tA\t1\tLowQual\t.\tGT\t0/1\t0/0\t0/1\n')
                out = root / 'out'
                run(argparse.Namespace(output_dir=out, reference=ref, target='chr22:1-10',
                                       model_dir=model, map_dir=root, annotation_vcf=vcf,
                                       kid_id='K', dad_id='F', mom_id='M'))
                bed = out / 'trio.dna-methylation.all_cpgs.sorted.bed.gz'
                frame = pl.read_csv(bed, separator='\t', comment_prefix='##')
                self.assertEqual(frame.height, 5)  # includes uncovered reference CpG
                self.assertFalse(any(c.endswith('_count') for c in frame.columns))
                first = frame.row(0, named=True)
                self.assertTrue(first['cpg_is_allele_specific_kid'])
                self.assertFalse(first['cpg_is_allele_specific_dad'])
                for label, expected in [('kid_pat', .2), ('kid_mat', .8), ('dad_A', .2),
                                        ('dad_B', .8), ('mom_C', .2), ('mom_D', .8)]:
                    self.assertEqual(first[f'methylation_level_{label}_model'], expected if resolved else None)
                    self.assertEqual(frame.tail(4)[f'methylation_level_{label}_model'].null_count(), 4)
                self.assertEqual(first['methylation_level_kid_model'], .5)
                self.assertIsNone(frame.row(-1, named=True)['methylation_level_kid_model'])
                with pysam.TabixFile(str(bed)) as tabix:
                    self.assertEqual(len(list(tabix.fetch('chr22'))), 5)
                tracks = list(out.glob('*.bw'))
                self.assertEqual(len(tracks), 9)
                for track in tracks:
                    with pyBigWig.open(str(track)) as bw:
                        self.assertEqual(bw.chroms(), {'chr22': 10})


if __name__ == '__main__':
    unittest.main()
