# Crossover slide notes

Slide numbers refer to the [November crossover deck][november]. The
[run instructions and results](README.md) give the complete chr22 comparison.

| Slide | Code or evidence to cite |
| --- | --- |
| 6–7 | [Runner](reproduce.sh) and [trio VCF quality filter][quality]. Cite 39 → 4 switches and 5,995 → 2,303 variants. |
| 9–10 | [Read audit][reads] and the four-SNP mask in the [runner](reproduce.sh). Cite the conflicting bases and the unchanged paternal bracket after rephasing. |
| 11 | [Spacing and BED overlap code](summarize.py) and the public-mask branch in the [runner](reproduce.sh). Cite the three quality-filtered calls flagged by spacing outside both public tracks. |
| 12 | [Read audit][reads], [original paternal HiPhase VCF][hiphase] and [GRCh38 inheritance map][iht]. Use the two pedigree comparisons with their respective brackets. |
| 13–14 | The pinned source references below for genotype and marker filters. The DP ≥18 run tests one Recombinator cutoff in WhatsHap. |

### Slide 13 · Recombinator

The pinned marker selector uses site call rate ≥0.95, accepts records for
which cyvcf2 reports no failing FILTER, and screens non-SNP records by REF
length ≤3, one ALT and ALT length ≤3. It supports an optional exclusion BED.
Within a family it requires called genotypes with GQ ≥20 and DP ≥18 by
default. For the unphased path, the transmitting heterozygous parent's ALT
fraction must be strictly between 0.31 and 0.69, and the children's fractions
must be >0.31 and ≤0.69 where evaluated. The fully phased path branches before
that allele-balance screen. [CLI defaults][recomb-cli] · [Marker selector][recomb-variants]

The downstream caller defaults to at least 20 informative sites per block and
also suppresses adjacent runs of rapid state flipping. These are marker-count
rules. [Caller configuration][recomb-block-cli] · [Block filtering][recomb-blocks]

### Slide 14 · Platinum Pedigree Inheritance

At Tapestry's pinned `gtg-ped-map` revision, defaults are QUAL ≥20, family
member depth ≥10 and `--run 10`. The code skips length-changing alleles, rejects
missing family-member alleles, and also requires each sample's depth to be
within **one** chromosome-specific standard deviation of its mean, as enforced
by the implementation. It tracks informative alleles through the pedigree.
The marker-masking rule examines same-state runs in marker order and masks a
position when the counts in both directions, each including that position, are
below `--run`. [CLI defaults][gtg-cli] · [Record filter][gtg-variants] ·
[Depth filter][gtg-depth] · [Marker selection][gtg-markers] · [Run masking][gtg-run]

### July deck · Slide 6

The July slide used concordance 0.57 to predict a minority segment under a single
crossover and observed about 82% of maternal mismatches in that segment. The
archived scripts ranked mismatch concentration under a roughly uniform
heterozygous-marker density assumption. They had no workflow callers beyond
their standalone wrapper and can be recovered from Git history for that example.
[July slide 6][july] · [Historical wrapper][heuristic-shell] · [Historical scoring code][heuristic-py]

[november]: https://docs.google.com/presentation/d/1DNbL8Fv1w71YkX4RbH7z1vQOdtUXww-2XkzPDmkIMVo/edit
[july]: https://docs.google.com/presentation/d/1WNP8aAWJM2DS_gKk1CvkkQO4UIj4-K_WlG8FwQ9KZNU/edit#slide=id.g3eeef4e7a8e_0_10
[quality]: ../../src/prefilter_whatshap_vcf.py
[reads]: ../../src/filter_whatshap_crossovers_by_reads.py
[iht]: ../../trio_dev_data/input/pedigree/CEPH1463.GRCh38.iht.sorted.txt
[hiphase]: ../../trio_dev_data/input/pedigree/NA12877.GRCh38.deepvariant.glnexus.phased.vcf.gz
[heuristic-shell]: https://github.com/quinlan-lab/tapestry/blob/19df1034b45cf15c957e2218e97da2fd7cee125f/find_clean_single_crossover.sh
[heuristic-py]: https://github.com/quinlan-lab/tapestry/blob/19df1034b45cf15c957e2218e97da2fd7cee125f/src/util/find_clean_single_crossover.py
[recomb-cli]: https://github.com/quinlan-lab/recombinator/blob/a164c2ea5e91debacbe658e85fa38e89ebafad05/recombinator/recombinator.py#L76-L85
[recomb-variants]: https://github.com/quinlan-lab/recombinator/blob/a164c2ea5e91debacbe658e85fa38e89ebafad05/recombinator/recombinator.py#L425-L493
[recomb-block-cli]: https://github.com/quinlan-lab/recombinator/blob/a164c2ea5e91debacbe658e85fa38e89ebafad05/recombinator/recombinator.py#L32-L49
[recomb-blocks]: https://github.com/quinlan-lab/recombinator/blob/a164c2ea5e91debacbe658e85fa38e89ebafad05/recombinator/recombinator.py#L124-L194
[gtg-cli]: https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L48-L58
[gtg-variants]: https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L156-L183
[gtg-depth]: https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L217-L236
[gtg-markers]: https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L1082-L1119
[gtg-run]: https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L935-L967
