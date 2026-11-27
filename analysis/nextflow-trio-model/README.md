# Crossover-aware trio model methylation on chr22

This analysis records the experiments behind
[the phase-block and methylation deck](https://docs.google.com/presentation/d/1fMsmfMTXjtrsqQmHAlmu7VouI0DEg1I7fLOuYRM2h8c/edit).
It follows the [crossover-filter experiments](../crossover_experiments/README.md)
on the same GRCh38 chr22:20,500,000–21,500,000 fixture: child NA12879, father
NA12877, mother NA12878. See [slide notes](slides.md) for the figure-to-code map.

The first analysis asks which candidate crossovers merit follow-up. This one
uses the provisional trio quality filter and 10 kb distance screen to split
bit-vector matching into resolved intervals, then assigns model methylation to
parental haplotypes. An isolated call is not automatically a verified crossover.
The retained paternal candidate has the read support described in the first
analysis; the distance screen itself does not perform that read audit.

## Recorded results

These values were recovered from the existing `nextflow-trio-model` run records
and regenerated with `summarize.py`. They describe this fixture, not whole chr22.

| Stage | Observation |
| --- | --- |
| Variant prefilter | 2,303 records; PASS or `.`, QUAL ≥20, called trio genotypes, GQ ≥20 and DP ≥10 in each member |
| WhatsHap phase blocks | One block per sample: child 20,500,019–21,349,175 (447 variants); father 20,500,019–21,485,412 (949); mother 20,500,927–21,485,412 (1,065) |
| Crossovers | Four parent-specific switches: three paternal, one maternal; three flagged at 10 kb |
| Cluster exclusion | Both parents unresolved across BED `[20,502,762, 20,510,653)`, including the gap between flagged brackets |
| Isolated paternal bracket | One-based 20,959,352–20,961,970; paternal assignment unresolved across BED `[20,959,351, 20,961,970)` |
| All-CpG publication | 19,206 rows; nine model BigWigs (six haplotype and three combined tracks) |

The resolved maps have these zero-based, half-open intervals:

| Parent | Start | End | Haplotype | Informative SNVs | Concordance | Mismatches |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| Paternal | 20,500,018 | 20,502,762 | A | 3 | 1.0 | 0 |
| Paternal | 20,510,653 | 20,959,351 | A | 620 | 1.0 | 0 |
| Paternal | 20,961,970 | 21,349,175 | B | 189 | 1.0 | 0 |
| Maternal | 20,500,926 | 20,502,762 | C | 1 | 1.0 | 0 |
| Maternal | 20,510,653 | 21,349,175 | D | 918 | 1.0 | 0 |

A/B label the father's fixed haplotypes; C/D label the mother's. The paternal
assignment changes A→B across the isolated bracket. The maternal D interval
continues across it. Perfect concordance in the very short left-edge intervals
rests on only three paternal SNVs and one maternal SNV. It is not independent
evidence of reliable transmission there. Each informative SNV belongs to its
sample's own `FORMAT/PS`, including when phase-set spans nest.

| Non-null model measurements | Child paternal | Child maternal | Father A | Father B | Mother C | Mother D | Child combined | Father combined | Mother combined |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CpG rows | 5,931 | 5,073 | 12,205 | 12,195 | 12,410 | 12,527 | 16,759 | 18,524 | 18,826 |

The workflow removes old HP/PS BAM tags while preserving MM/ML, haplotags against
pedMEC, then runs pb-CpG-tools v3.0.0 in **model** mode with minimum coverage 10
and MAPQ 1. Haplotype measurements remain null where that parent's assignment is
unresolved; missing measurements are never zero-filled. Combined calls include
untagged reads and remain available outside resolved maps. The all-CpG table
spans the CG dinucleotide; BigWigs use the cytosine base. Stored methylation
levels are fractions; the figure displays percentages.

## Reproduce

Run from the repository root. Use the repository `.venv` dependencies, WhatsHap
2.8, Nextflow, bcftools (the recorded experiment used 1.23), bgzip, tabix,
samtools, and pb-CpG-tools v3.0.0 (`aligned_bam_to_cpg_scores`) on `PATH`.
The latter is not in the pedigree Docker image; use a trio-capable Linux/CHPC
environment. The checked-in indexed BAMs and joint VCF are inputs. Prepare the
reference once:

```bash
./trio_dev_data_get_ref.sh trio_dev_data
./analysis/nextflow-trio-model/reproduce.sh
```

The runner invokes the default trio workflow on the 1 Mb fixture with
explicit scientific settings, records its command, tool versions, log, and
Nextflow trace, generates per-sample block lists, runs the existing scientific
fixture checks, and generates the figures. It refuses an existing destination.
For a repeat, supply a fresh path:

```bash
./analysis/nextflow-trio-model/reproduce.sh trio_dev_data/output/nextflow-trio-model-rerun
```

To regenerate only the summary and figures from an existing published run:

```bash
.venv/bin/python analysis/nextflow-trio-model/summarize.py \
  --run-dir trio_dev_data/output/nextflow-trio-model \
  --output-dir .test-work/trio-deck-figures
```

This reads `per_chrom_phased/chr22.crossovers.tsv`, both
`per_chrom_hap_maps/chr22/*.hap-map.tsv` tables and mismatch BEDs, the final
all-CpG BED, and `phase-blocks/NA1287{7,8,9}.blocks.tsv`. If block lists are
missing, regenerate each with `whatshap stats --sample SAMPLE --block-list
DESTINATION per_chrom_phased/chr22.phased.vcf.gz`. Omitting `--output-dir`
replaces the derived figures and summaries in the run directory. Regenerated
figures use the deck's data and filenames; they do not preserve the original
PNG styling byte for byte.

## Find the results

```text
trio_dev_data/output/nextflow-trio-model/
├── per_chrom_phased/           Indexed pedMEC VCF, raw report, screened TSV/BED
├── per_chrom_hap_maps/chr22/   Parental maps and separate mismatch tracks
├── per_chrom_model/chr22/      Re-haplotagged BAMs and raw model BEDs on rerun
├── per_chrom_methylation/chr22/  All-CpG BED/index, nine BigWigs, QC
├── phase-blocks/              Sample block lists, stats, phase-blocks.png
├── crossover-brackets.png     Cluster and isolated candidate, separate zooms
├── resolved-intervals.png     Full maps and short left-edge interval zoom
├── parent_phased_methylation.png  Assignments and child model measurements
├── resolved-intervals.tsv     Coordinates, labels, concordance, SNVs, mismatches
└── summary.json               Row and non-null measurement counts
```

Existing local runs may omit intermediate BAMs/model BEDs and historical
command/version records. The runner creates those on a fresh run; summary
generation alone cannot recover missing execution provenance.

The existing [`tests/run_nextflow_trio_smoke.sh`](../../tests/run_nextflow_trio_smoke.sh)
checks the fixture's candidate brackets, A→B assignment, unresolved cluster and
bracket, model-call parity, null preservation, and indexes. It also exercises
distance cutoffs of 24 and 23 bp. Those sensitivity runs are test cases, not
additional experiments shown in this deck.

The deck's next steps—whole-chr22/genome-wide evaluation, manuscript integration,
and scans for de novo methylation or compound heterozygosity—remain future work.
The fixture does not establish genome-wide precision, sensitivity, or optimal
filter thresholds.
