"""Write a VCF of markers at which a kid has ≥1 deduced haplotype, for IGV.

Input is a haplotype-map "markers" table such as
    /scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/
        CEPH1463.GRCh38/CEPH1463.GRCh38.markers.sorted.txt
whose first (header) line names the columns

    #chom pos <sampleID> <sampleID> ...

and whose data lines carry, per sample, a diploid haplotype call whose two
entries are either a haplotype LABEL (a letter A, B, C, ... naming the
grandparental-of-origin haplotype) or '?' (that haplotype could not be
deduced at this marker), e.g.

    chr1 16297 A/B C/D E/F G/H I/J K/L A|? A|? B|? ...

The defining/founder samples appear as e.g. "A/B" (their two haplotypes are the
labels by definition); the samples whose haplotypes are being INFERRED appear as
e.g. "?|?", "?|C" or "A|?" (a '|' separates the two, and a letter means that
haplotype was deduced at this marker).

A haplotype is "deduced" for a sample at a marker when its entry is a letter,
not '?'. This tool emits, for the requested kid sample(s), the markers at which
AT LEAST ONE of the kid's two haplotypes is deduced (entry != "?|?"), so they
can be shown as an IGV variant track. Pass --sample more than once to pool
several kids: a marker is kept when at least one named kid has ≥1 deduced
haplotype there. Each kept record carries the kid call(s) in its INFO/HAP field
so a hover in IGV shows which haplotype(s) were deduced.

The output has no REF/ALT nucleotides (the markers table does not carry them),
so REF is written as 'N' and ALT as '.'; the record is a position-only tick.

The output filename is DERIVED from the markers filename and the sample(s):
    <dataset>.<sample[+sample...]>.deduced-haplotype-markers.vcf
e.g. "CEPH1463.GRCh38.markers.sorted.txt" + --sample NA12879 gives
"CEPH1463.GRCh38.NA12879.deduced-haplotype-markers.vcf". It is written next to
the markers file unless --out-dir is given.

Usage:
    PYTHONPATH=src:src/util .venv/bin/python \\
        src/util/write_deduced_haplotype_markers_vcf.py \\
        --markers CEPH1463.GRCh38.markers.sorted.txt \\
        --sample  NA12879

By default the VCF is then bgzip-compressed and tabix-indexed (needs bgzip and
tabix on PATH); pass --no-index to leave a plain .vcf.
"""

import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess


def derive_out_path(markers_path, samples, out_dir):
    """Derive the output VCF path from the markers filename and the sample(s).

    Strips the …[.markers][.sorted].txt[.gz] tail off the markers basename to
    recover the dataset prefix, then appends the '+'-joined sample IDs, e.g.
    (CEPH1463.GRCh38.markers.sorted.txt, [NA12879]) ->
        <out_dir>/CEPH1463.GRCh38.NA12879.deduced-haplotype-markers.vcf
    """
    dataset = os.path.basename(markers_path)
    for suffix in (r"\.gz$", r"\.txt$", r"\.sorted$", r"\.markers$"):
        dataset = re.sub(suffix, "", dataset)
    if out_dir is None:
        out_dir = os.path.dirname(markers_path) or "."
    sample_tag = "+".join(samples)
    return os.path.join(out_dir, f"{dataset}.{sample_tag}.deduced-haplotype-markers.vcf")


def _open(path):
    """Open a possibly gzipped text file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _deduced_alleles(entry):
    """The deduced (non-'?') haplotype labels in a diploid entry like '?|C'.

    Entries use '|' between the inferred sample's two haplotypes; the founder
    columns use '/'. Both are accepted so the same test works on any column.
    """
    return [a for a in entry.replace("/", "|").split("|") if a != "?"]


def select_markers(markers_path, samples, logger):
    """Yield (chrom, pos, info) for markers where a kid has ≥1 deduced haplotype.

    info is the VCF INFO string, e.g. "HAP=NA12879:?|C" (or a comma-joined list
    across several deduced kids), naming each kept kid and its raw call.
    """
    with _open(markers_path) as f:
        header = f.readline().rstrip("\n").split()
        # header[0] = '#chom', header[1] = 'pos', header[2:] = sample IDs, in the
        # same order as the genotype fields on each data line.
        col_of = {sample_id: i for i, sample_id in enumerate(header[2:])}
        missing = [s for s in samples if s not in col_of]
        if missing:
            raise SystemExit(
                f"sample(s) not found in markers header: {', '.join(missing)}\n"
                f"available: {', '.join(header[2:])}"
            )

        n_kept = 0
        for line in f:
            fields = line.split()
            if not fields:
                continue
            calls = fields[2:]
            deduced = [
                f"{s}:{calls[col_of[s]]}"
                for s in samples
                if _deduced_alleles(calls[col_of[s]])
            ]
            if not deduced:
                continue
            n_kept += 1
            yield fields[0], fields[1], "HAP=" + ",".join(deduced)
        logger.info(f"Selected {n_kept} marker(s) with ≥1 deduced haplotype "
                    f"for: {', '.join(samples)}")


def write_vcf(records, out_path, samples):
    """Write selected (chrom, pos, info) records to a minimal VCFv4.2."""
    with open(out_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(
            '##INFO=<ID=HAP,Number=.,Type=String,Description='
            '"Deduced haplotype call(s) as sampleID:call, for kid sample(s) '
            f'{",".join(samples)}; a letter is a deduced grandparental-of-origin '
            'haplotype, ? is undeduced">\n'
        )
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, info in records:
            f.write(f"{chrom}\t{pos}\t.\tN\t.\t.\t.\t{info}\n")


def bgzip_and_index(vcf_path, logger):
    """bgzip-compress and tabix-index the VCF if the tools are available."""
    if not (shutil.which("bgzip") and shutil.which("tabix")):
        logger.warning("bgzip/tabix not found on PATH; leaving plain VCF "
                       f"'{vcf_path}' (pass --no-index to silence this)")
        return
    subprocess.run(["bgzip", "--force", vcf_path], check=True)
    subprocess.run(["tabix", "--force", "--preset", "vcf", f"{vcf_path}.gz"],
                   check=True)
    logger.info(f"Wrote, compressed and indexed: '{vcf_path}.gz' (+ .tbi)")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--markers", required=True,
                        help="Path to the markers table (…markers.sorted.txt[.gz])")
    parser.add_argument("--sample", required=True, action="append", dest="samples",
                        metavar="SAMPLE_ID",
                        help="Kid sample ID (e.g. NA12879); repeat to pool kids")
    parser.add_argument("--out-dir", default=None,
                        help="Directory for the derived-name output VCF "
                             "(default: alongside the markers file)")
    parser.add_argument("--no-index", action="store_true",
                        help="Leave a plain .vcf; do not bgzip/tabix it")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    logger = logging.getLogger(__name__)

    out_path = derive_out_path(args.markers, args.samples, args.out_dir)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    records = select_markers(args.markers, args.samples, logger)
    write_vcf(records, out_path, args.samples)
    if args.no_index:
        logger.info(f"Wrote: '{out_path}'")
    else:
        bgzip_and_index(out_path, logger)


if __name__ == "__main__":
    main()
