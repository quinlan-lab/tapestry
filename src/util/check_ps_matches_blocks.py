"""Verify that the pedmec-phased VCF carries PS tags that match the blocks TSV.

The phase-set assignment fix in phasing_trio.py assigns each SNV to a phase
block by the per-sample PS (phase-set) FORMAT tag rather than by genomic-range
overlap. That is only sound if, for the production VCF:

  1. PS is actually present on phased genotypes, and
  2. every PS value equals a `phase_set` listed for that sample in the WhatsHap
     blocks TSV (get_phase_blocks renames `phase_set` -> phase_block_id).

This script checks both for one sample. Run it on CHPC against the real
pedmec-phasing output BEFORE trusting the fix / before the regression run:

    PYTHONPATH=src:src/util .venv/bin/python src/util/check_ps_matches_blocks.py \\
        --vcf  <pedmec-phasing>/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz \\
        --blocks_tsv <pedmec-phasing>/...phased.NA12878.blocks.tsv \\
        --sample NA12878

Exit code is non-zero if PS is missing on phased het sites or any PS does not
match a block (i.e. the precondition for the fix does not hold).
"""

import argparse
import logging
import sys

from cyvcf2 import VCF  # type: ignore

from phasing import is_snv_het
from phasing_trio import get_phase_blocks, _read_ps


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", required=True, help="WhatsHap pedmec-phased trio VCF")
    parser.add_argument("--blocks_tsv", required=True, help="WhatsHap blocks TSV for the sample")
    parser.add_argument("--sample", required=True, help="Sample id to check (e.g. NA12878)")
    parser.add_argument("--max_report", type=int, default=20,
                        help="Max number of offending sites to print (default: 20)")
    parser.add_argument("--progress_every", type=int, default=1_000_000,
                        help="Log progress every N variants scanned (default: 1,000,000)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    logger = logging.getLogger(__name__)

    block_ids = set(get_phase_blocks(args.blocks_tsv, args.sample)["phase_block_id"].to_list())
    logger.info(f"{len(block_ids)} phase blocks in {args.blocks_tsv} for {args.sample}")

    vcf = VCF(args.vcf, strict_gt=True)
    si = vcf.samples.index(args.sample)

    try:
        vcf.get_header_type("PS")
    except KeyError:
        logger.error("VCF does not declare a 'PS' FORMAT tag. The PS-based fix cannot be used.")
        vcf.close()
        sys.exit(2)

    n_phased_het = 0
    n_no_ps = 0
    n_ps_not_in_blocks = 0
    examples_no_ps = []
    examples_unmatched = []

    n_total = 0
    last_chrom = None
    for variant in vcf:
        # Progress: cyvcf2 streams the whole VCF, so this loop is the slow part.
        # Log every args.progress_every variants and on each new contig.
        n_total += 1
        if variant.CHROM != last_chrom:
            last_chrom = variant.CHROM
            logger.info(f"... {variant.CHROM} (variants seen: {n_total:,}; "
                        f"phased het: {n_phased_het:,})")
        elif n_total % args.progress_every == 0:
            logger.info(f"... {variant.CHROM}:{variant.POS} (variants seen: "
                        f"{n_total:,}; phased het: {n_phased_het:,})")

        if not variant.is_snp or len(variant.ALT) != 1:
            continue
        if not is_snv_het(variant, si):
            continue
        gt = variant.genotypes[si]
        if not gt[2]:  # require phased
            continue
        n_phased_het += 1
        ps = _read_ps(variant.format("PS"), si)
        if ps is None:
            n_no_ps += 1
            if len(examples_no_ps) < args.max_report:
                examples_no_ps.append(f"{variant.CHROM}:{variant.POS}")
        elif ps not in block_ids:
            n_ps_not_in_blocks += 1
            if len(examples_unmatched) < args.max_report:
                examples_unmatched.append(f"{variant.CHROM}:{variant.POS} PS={ps}")
    vcf.close()

    logger.info(f"phased het sites for {args.sample}: {n_phased_het}")
    logger.info(f"  with no PS tag:           {n_no_ps}")
    logger.info(f"  PS not matching a block:  {n_ps_not_in_blocks}")
    if examples_no_ps:
        logger.info("  examples missing PS: " + ", ".join(examples_no_ps))
    if examples_unmatched:
        logger.info("  examples unmatched PS: " + ", ".join(examples_unmatched))

    if n_phased_het == 0:
        logger.error("No phased het sites found for this sample; cannot verify.")
        sys.exit(2)
    if n_no_ps or n_ps_not_in_blocks:
        logger.error("FAIL: PS precondition does not hold for the fix.")
        sys.exit(1)
    logger.info("PASS: every phased het site has a PS that matches a block. The PS-based fix is sound.")


if __name__ == "__main__":
    main()
