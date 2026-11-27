#!/usr/bin/env python3
"""Assign chromosome-local model calls using screened trio maps and expand CpGs."""
import argparse
import json
import logging
import re
from pathlib import Path

import polars as pl
import pysam
import pyBigWig

from get_meth_hap1_hap2 import read_meth_hap1_hap2, read_meth_level
from phase_meth_to_parent_haps import phase_meth_to_parent_haplotypes
from expand_to_all_cpgs_trio import (
    expand_meth_to_all_cpgs, compute_proximity_to_mismatched_heterozygous_sites,
    get_joint_called_variants, label_with_variants, label_cpgs_as_allele_specific,
)
from write_data import add_bigwig_entries, write_methylation


def clip(frame, chrom, start, end):
    column = 'chrom' if 'chrom' in frame.columns else 'chromosome'
    return frame.filter((pl.col(column) == chrom) &
                        (pl.col('start') >= start) & (pl.col('start') < end))


def run(args):
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    with pysam.FastaFile(str(args.reference)) as ref:
        chrom = args.target.split(':')[0]
        length = ref.get_reference_length(chrom)
        start, end = 0, length
        if ':' in args.target:
            left, right = args.target.split(':')[1].split('-')
            start, end = int(left) - 1, int(right)
        if end > length:
            raise ValueError('Target exceeds reference chromosome length')
        # Fetch the next base to recognize a CpG whose C is the last target base.
        sequence = ref.fetch(chrom, start, min(end + 1, length)).upper()
    coordinates = [(chrom, start + m.start(), start + m.start() + 1)
                   for m in re.finditer('CG', sequence) if start + m.start() < end]
    reference = pl.DataFrame(coordinates, schema={'chrom': pl.String, 'start': pl.Int64,
                                                 'end': pl.Int64}, orient='row')
    models = []
    combined = []
    for role in ('kid', 'dad', 'mom'):
        models.append(clip(read_meth_hap1_hap2('model', args.model_dir / f'{role}.hap1.bed.gz',
                                              args.model_dir / f'{role}.hap2.bed.gz'), chrom, start, end))
        combined.append(clip(read_meth_level(args.model_dir / f'{role}.combined.bed.gz', 'model'),
                             chrom, start, end).rename({
                                 'chromosome': 'chrom', 'total_read_count': f'total_read_count_{role}',
                                 'methylation_level': f'methylation_level_{role}_model'}))
    maps = [pl.read_csv(args.map_dir / f'{side}.hap-map.tsv', separator='\t', schema_overrides={
                'chrom': pl.String, 'start': pl.Int64, 'end': pl.Int64,
                f'{side}_haplotype': pl.String, 'haplotype_concordance': pl.Float64,
                'num_het_SNVs_in_parent': pl.Int64})
            for side in ('paternal', 'maternal')]
    phased = phase_meth_to_parent_haplotypes(None, models[0], None, models[1],
                                           None, models[2], *maps)
    unphased = combined[0]
    for frame in combined[1:]:
        unphased = unphased.join(frame, on=['chrom', 'start', 'end'], how='full', coalesce=True)
    frame = expand_meth_to_all_cpgs(reference, unphased, phased)
    frame = compute_proximity_to_mismatched_heterozygous_sites(
        frame, args.map_dir / f'{args.kid_id}.bit-vector-sites-mismatches.paternal.bed',
        args.map_dir / f'{args.kid_id}.bit-vector-sites-mismatches.maternal.bed')
    variants = get_joint_called_variants(args.kid_id, args.dad_id, args.mom_id, args.annotation_vcf)
    frame = label_cpgs_as_allele_specific(label_with_variants(frame, variants))
    write_methylation(frame, str(out / 'trio.dna-methylation.all_cpgs.bed'),
                      source=f'{__file__} with args {vars(args)}')
    # BigWigs retain cytosine-base coordinates; the annotated table spans CG.
    specs = [(args.kid_id, 'pat', 'kid_pat'), (args.kid_id, 'mat', 'kid_mat'),
             (args.dad_id, 'A', 'dad_A'), (args.dad_id, 'B', 'dad_B'),
             (args.mom_id, 'C', 'mom_C'), (args.mom_id, 'D', 'mom_D')]
    tracks = [(uid, label, person, phased) for uid, label, person in specs]
    tracks += [(uid, 'combined', role, data) for uid, role, data in
               zip((args.kid_id, args.dad_id, args.mom_id), ('kid', 'dad', 'mom'), combined)]
    for uid, label, person, data in tracks:
        with pyBigWig.open(str(out / f'{uid}.dna-methylation.{label}.model.hg38.bw'), 'w') as bw:
            bw.addHeader([(chrom, length)])
            add_bigwig_entries(bw, data, person, 'model')
    qc = {'target': args.target, 'measurement': 'model', 'rows': frame.height,
          'paternal_unresolved_rows': frame['paternal_haplotype'].null_count(),
          'maternal_unresolved_rows': frame['maternal_haplotype'].null_count()}
    (out / 'trio.model.qc.json').write_text(json.dumps(qc, indent=2) + '\n')
    logging.info('Wrote model methylation, nine BigWigs, index and QC: %s', qc)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('model-dir', 'map-dir', 'reference', 'annotation-vcf', 'output-dir'):
        parser.add_argument('--' + name, required=True, type=Path)
    for name in ('target', 'kid-id', 'dad-id', 'mom-id'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    logging.info('Inputs: %s', vars(args))
    run(args)


if __name__ == '__main__':
    main()
