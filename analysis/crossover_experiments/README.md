# WhatsHap crossover experiments on chr22

This is the analysis for the
[crossover-filter slide deck](https://docs.google.com/presentation/d/1DNbL8Fv1w71YkX4RbH7z1vQOdtUXww-2XkzPDmkIMVo/edit).
The subsequent [crossover-aware trio methylation analysis](../nextflow-trio-model/README.md)
records the phase blocks, resolved bit-vector matches, and methylation figures
in the second deck. See the [analysis index](../README.md) for both output branches.

The checked-in GRCh38 fixture covers chr22:20,500,000–21,500,000 for NA12879
(child), NA12877 (father) and NA12878 (mother). These are candidate switches from
WhatsHap 2.8 with bcftools 1.23. A raw WhatsHap row can report a switch for both
parents. Each parent contributes one switch to the table below. See the
[slide notes](slides.md) for code references and upstream thresholds.

| Phasing input | VCF records | Paternal | Maternal | Within 10 kb of another call | Overlap with either public track |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline | 5,995 | 17 | 22 | 35 | 38 |
| Trio quality filter | 2,303 | 3 | 1 | 3 | 0 |
| Quality + public mask | 1,507 | 2 | 1 | 2 | 0 |
| Quality + four-SNP mask | 2,299 | 1 | 0 | 0 | 0 |
| Quality + DP ≥18 | 1,658 | 0 | 0 | 0 | 0 |

- **Quality filtering** requires FILTER PASS or `.`, QUAL ≥20, and called
  genotypes with GQ ≥20 and DP ≥10 in all three samples. The retained paternal
  bracket at 20,959,352–20,961,970 was identical in the baseline, quality, public
  mask and four-SNP mask runs. Raising the depth floor to 18 removed its right
  marker, where the child had DP 17, and yielded zero events in the interval.
- **Public annotations and spacing** answer different questions. All 35
  clustered baseline switches overlapped both public tracks. Three isolated
  baseline switches also overlapped both. After quality filtering, the three
  clustered switches were near the left edge of the fixture and outside both
  tracks. The two distinct brackets were 24 bp apart, and one bracket switched
  both parental states. Masking all UCSC segmental duplications and ENCODE
  blacklist intervals before rerunning WhatsHap left two of these switches.
  Spacing flags can catch clusters in regions that these annotations miss, and
  can leave isolated calls in annotated regions. The overlap measurement uses
  the full marker bracket. The mask excludes VCF records by POS before phasing.
- **Read evidence** flagged the three clustered quality-filtered switches.
  At chr22:20,503,731, overlapping alignments assigned conflicting bases to
  15/23 child, 40/74 father and 47/86 mother read names with usable evidence.
  Removing the four inspected SNPs at 20,503,731, 20,503,741, 20,503,743 and
  20,503,755 before phasing removed the cluster. The retained paternal
  bracket had 17 G–A and 17 A–G primary spanning reads in the father and
  eight recombinant G–G primary spanning reads in the child. The unordered
  paternal allele pairs agree with the checked-in HiPhase VCF.
- **Pedigree comparison** supports a paternal change in the same neighborhood.
  The broad candidate in the inheritance review notes (20,955,482–20,983,035,
  assembly unverified in those notes) overlaps the WhatsHap bracket. The
  checked-in GRCh38 inheritance map changes the child's paternal label from
  E to G between 20,982,976 and 20,983,035. That bracket starts 21,006 bp after
  the right WhatsHap marker. The maternal label remains I. Treat the map and
  review notes as comparison evidence. [Inheritance review notes][review] ·
  [Checked-in inheritance map][iht]

The spacing rule is exploratory. It flags the nearest bracket gap ≤10,000 bp
within one child and chromosome, pooling both parents and using zero for
overlapping brackets. A single bracket with no neighbor gets a blank gap. The
read audit uses MAPQ and baseQ ≥20 and counts conflicting alignments by query
name. It flags an anchor with ≥10 covering molecules, ≥3 conflicts and ≥20%
conflict, and requires ≥3 clean spanning reads per parental haplotype, ≥90%
parental concordance and ≥3 recombinant child reads for read support. These
cutoffs and the four-SNP mask were examined in the same small window. The
window lacks a complete, independently adjudicated event set, so precision,
sensitivity and genome-wide cutoffs remain unmeasured.

## Reproduce the five branches

From the repository root, use the checked-in fixture, the repository `.venv`,
WhatsHap 2.8, bcftools, bgzip and samtools. Use a fresh output directory for
each rerun. The experiment phases the trio and audits the indexed input BAMs
directly, using aligned bases for evidence. It does not require methylation
calling.

```bash
./trio_dev_data_get_ref.sh trio_dev_data
mkdir -p .test-work/crossover-annotations
curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz \
  -o .test-work/crossover-annotations/genomicSuperDups.txt.gz
curl -fL https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz \
  -o .test-work/crossover-annotations/hg38-blacklist.v2.bed.gz
./analysis/crossover_experiments/reproduce.sh \
  .test-work/crossover-annotations/genomicSuperDups.txt.gz \
  .test-work/crossover-annotations/hg38-blacklist.v2.bed.gz
```

SHA-256 for the inputs used in the table:

```text
ce3ee1ca39356238f7aee438a40a88b4f1b9d80b316b263e16fb12402212d10f  chr22.fa
3f4388ed3eb7462034e89c1f45314160bdb14558502092326063105c095061df  genomicSuperDups.txt.gz
c92e763af17271446194991e71917ac220593a5a3d40a06667be24178ef08cf2  hg38-blacklist.v2.bed.gz
```

## Find the results

By default, the runner writes all five branches to
`trio_dev_data/output/crossover-experiments/`. Pass a new directory as the
optional third argument to write the same layout elsewhere. Start with
`quality/crossovers.tsv` for the four candidates and
`quality/read-evidence/calls.tsv` for the evidence behind each classification.

```text
trio_dev_data/output/crossover-experiments/
├── baseline/                 Original VCF
├── quality/                  Trio quality filter
│   └── read-evidence/        Classifications and per-anchor counts
├── public-mask/              Quality filter plus public region masks
├── four-snp-mask/            Quality filter with four suspect SNPs removed
└── depth18/                  Quality filter with a higher depth floor
```

Each branch contains `crossovers.tsv`, `crossovers.bed`, the raw
`recombination.txt` report and the indexed `phased.vcf.gz`. The other files in
each branch record the phasing command, logs and inputs for reproduction.
Shared annotation BEDs and filtered input VCFs are at the experiment root.
TSV marker positions are one-based inclusive. BED is zero-based half-open.

The trio methylation workflow writes to sibling directories such as
`trio_dev_data/output/pedmec-phasing/` when run.

For IGV, select hg38 and load the local `crossovers.bed`, `phased.vcf.gz` and
the indexed fixture BAMs. Load `segdup.bed` and `blacklist.bed` from the
experiment output root to compare annotations. Start at
`chr22:20,502,763-20,510,653` for the cluster or
`chr22:20,959,352-20,961,970` for the paternal candidate. IGV can save a session
from the loaded tracks.

[review]: https://docs.google.com/document/d/1HGBhjqk4cqZaWIkh0GeAaBC_fcZYWRRnYO5cjiifKEY/edit
[iht]: ../../trio_dev_data/input/pedigree/CEPH1463.GRCh38.iht.sorted.txt
