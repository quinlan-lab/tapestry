"""Check scientific labels, unresolved spans, and indexes after the chr22 smoke run."""
import argparse
from pathlib import Path
import polars as pl
import pyBigWig
import pysam

from get_meth_hap1_hap2 import read_meth_level


def check(results, docker_runtime=False):
    output = results / 'per_chrom_methylation/chr22'
    bed = output / 'trio.dna-methylation.all_cpgs.sorted.bed.gz'
    frame = pl.read_csv(bed, separator='\t', comment_prefix='##', infer_schema_length=None)
    assert not any(c.endswith('_count') for c in frame.columns)
    assert frame['start_cpg'].min() >= 20499999
    assert frame['start_cpg'].max() < 21500000
    assert frame['start_cpg'].is_sorted()
    assert set(frame['paternal_haplotype'].drop_nulls()) == {'A', 'B'}
    assert set(frame['maternal_haplotype'].drop_nulls()) <= {'C', 'D'}
    assert frame.filter(pl.col('start_cpg') < 20959351)['paternal_haplotype'].drop_nulls().unique().to_list() == ['A']
    assert frame.filter(pl.col('start_cpg') >= 20961970)['paternal_haplotype'].drop_nulls().unique().to_list() == ['B']
    for side, labels in [('paternal', ('kid_pat', 'dad_A', 'dad_B')),
                         ('maternal', ('kid_mat', 'mom_C', 'mom_D'))]:
        unresolved = frame.filter(pl.col(f'{side}_haplotype').is_null())
        assert unresolved.height > 0
        for label in labels:
            assert unresolved[f'methylation_level_{label}_model'].null_count() == unresolved.height
        assert frame.filter(pl.col(f'{side}_haplotype').is_not_null()).height > 0
    early_span = ((20502762, 20503731, ('maternal',)) if docker_runtime else
                  (20502762, 20510653, ('paternal', 'maternal')))
    for start, end, sides in [early_span, (20959351, 20961970, ('paternal',))]:
        span = frame.filter(pl.col('start_cpg').is_between(start, end - 1))
        assert span.height > 0
        for side in sides:
            assert span[f'{side}_haplotype'].null_count() == span.height
    # Compare every resolved haplotype measurement with the corresponding caller BED.
    for role, labels in [('kid', ('kid_pat', 'kid_mat')), ('dad', ('dad_A', 'dad_B')),
                         ('mom', ('mom_C', 'mom_D'))]:
        for hap, label in enumerate(labels, 1):
            raw = read_meth_level(results / f'per_chrom_model/chr22/model/{role}.hap{hap}.bed.gz', 'model')
            observed = frame.filter(pl.col(f'methylation_level_{label}_model').is_not_null())
            assert observed.height > 0, label
            joined = observed.join(raw, left_on=['#chrom', 'start_cpg'], right_on=['chromosome', 'start'])
            assert joined.height == observed.height
            assert joined[f'methylation_level_{label}_model'].equals(joined['methylation_level'])
    with pysam.TabixFile(str(bed)) as tabix:
        assert list(tabix.fetch('chr22', 20959351, 20961970))
    for model_bed in (results / 'per_chrom_model/chr22/model').glob('*.bed.gz'):
        with pysam.TabixFile(str(model_bed)) as tabix:
            assert list(tabix.fetch('chr22'))
    tracks = list(output.glob('*.bw'))
    assert len(tracks) == 9
    for track in tracks:
        with pyBigWig.open(str(track)) as bw:
            assert 'chr22' in bw.chroms()
            assert bw.intervals('chr22')
    for bam in (results / 'per_chrom_model/chr22/haplotagged').glob('*.bam'):
        with pysam.AlignmentFile(str(bam)) as reads:
            assert reads.has_index()
            tagged = [r for r in reads.fetch('chr22', 20500000, 21500000) if r.has_tag('HP')]
            assert tagged
    print('PASS: model labels, unresolved nulls, caller parity, nine BigWigs and BED/BAM indexes')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--docker-runtime', action='store_true',
                        help='Check the pinned Linux/Python 3.11 fixture baseline')
    args = parser.parse_args()
    check(args.results, args.docker_runtime)
