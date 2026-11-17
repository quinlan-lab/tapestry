"""Apply the provisional trio-wide quality filter before WhatsHap phasing."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import pysam

EXPRESSION = 'QUAL>=20 && N_MISSING=0 && N_PASS(FMT/GQ>=20 & FMT/DP>=10)=3'


def prefilter(input_vcf: Path, output_vcf: Path, samples: list[str]) -> dict:
    if len(samples) != 3 or len(set(samples)) != 3:
        raise ValueError('Expected three distinct samples: child, father, mother')
    with pysam.VariantFile(str(input_vcf)) as reader:
        if not set(samples).issubset(reader.header.samples):
            raise ValueError('Trio samples are missing from the input VCF')
        for tag in ('GT', 'GQ', 'DP'):
            if tag not in reader.header.formats:
                raise ValueError(f'Missing FORMAT/{tag}')
        before = sum(1 for _ in reader)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    # Subset first: samples outside the trio must not affect N_PASS/N_MISSING.
    subset = ['bcftools', 'view', '-s', ','.join(samples), '-Ou', str(input_vcf)]
    select = ['bcftools', 'view', '-f', 'PASS,.', '-i', EXPRESSION,
              '-Oz', '-o', str(output_vcf)]
    print(f'Prefilter input: {input_vcf}; samples: {samples}; expression: {EXPRESSION}', flush=True)
    with subprocess.Popen(subset, stdout=subprocess.PIPE) as upstream:
        try:
            subprocess.run(select, stdin=upstream.stdout, check=True)
        finally:
            if upstream.stdout is not None:
                upstream.stdout.close()
            status = upstream.wait()
        if status:
            raise subprocess.CalledProcessError(status, subset)
    subprocess.run(['bcftools', 'index', '-f', '-t', str(output_vcf)], check=True)
    with pysam.VariantFile(str(output_vcf)) as reader:
        after = sum(1 for _ in reader)
    with input_vcf.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    report = {'enabled': True, 'input': str(input_vcf), 'input_sha256': digest,
              'output': str(output_vcf), 'samples': samples, 'filter_status': ['PASS', '.'],
              'expression': EXPRESSION, 'input_records': before, 'retained_records': after,
              'excluded_records': before - after, 'commands': [subset, select],
              'bcftools_version': subprocess.check_output(['bcftools', '--version'], text=True).splitlines()[0],
              'interpretation': 'Provisional quality filter; no calibrated crossover threshold or truth set.'}
    Path(str(output_vcf) + '.qc.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'Prefilter retained {after}/{before} records: {output_vcf}', flush=True)
    if not after:
        raise ValueError('Prefilter retained zero records; inspect its QC report before phasing')
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--samples', nargs=3, required=True, metavar='SAMPLE')
    args = parser.parse_args()
    prefilter(args.input, args.output, args.samples)


if __name__ == '__main__':
    main()
