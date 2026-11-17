# Crossover speaker notes

Copy the text inside each block into the corresponding slide’s speaker notes.
The links use full URLs so they travel with the notes. The Tapestry links fixed
to commit `3d5894a` will resolve on GitHub once that commit is pushed.

## November deck

[Discovering crossovers for methylation phasing in trios](https://docs.google.com/presentation/d/1DNbL8Fv1w71YkX4RbH7z1vQOdtUXww-2XkzPDmkIMVo/edit)

### Slide 6 · Baseline calls

```text
We ran pedigree-aware WhatsHap 2.8 on GRCh38 chr22:20,500,000–21,500,000 for child NA12879, father NA12877 and mother NA12878. The baseline used 5,995 VCF records and yielded 39 candidate transmission switches, 17 paternal and 22 maternal. A report row that switches both parents counts as two events.

Phasing command and reproduction results
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/analysis/crossover_experiments/reproduce.sh#L61-L89
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/analysis/crossover_experiments/README.md
```

### Slide 7 · Trio quality filtering

```text
Before phasing, we accepted VCF FILTER=PASS or “.”, required site QUAL ≥20, and required called genotypes with GQ ≥20 and DP ≥10 in all three trio members. The filter retained 2,303 of 5,995 records. Rerunning WhatsHap yielded four switches: three paternal and one maternal.

Trio VCF quality filter
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/src/prefilter_whatshap_vcf.py
```

### Slide 8 · Suspect variants

```text
The baseline maternal candidate at chr22:21,125,468–21,125,744 disappeared after trio quality filtering. This comparison shows that the call depends on the input variant set. Read evidence and phase continuity are needed to adjudicate an individual switch.

Reproduce the baseline and quality-filtered reports
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/analysis/crossover_experiments/README.md
```

### Slide 9 · Ambiguous alignment cluster

```text
Three of the four quality-filtered switches fell near chr22:20.50 Mb. At position 20,503,731, overlapping alignments assigned conflicting bases to 15 of 23 child read names, 40 of 74 father read names and 47 of 86 mother read names with usable evidence. The audit counts evidence by read name and requires mapping and base qualities of at least 20. These conflicts make the marker unreliable for locating a crossover.

Read audit
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/src/filter_whatshap_crossovers_by_reads.py#L110-L159
```

### Slide 10 · Four-SNP masking experiment

```text
We removed four inspected SNPs at chr22 positions 20,503,731, 20,503,741, 20,503,743 and 20,503,755 from the quality-filtered VCF, then reran WhatsHap. All three boundary switches disappeared. The paternal bracket at chr22:20,959,352–20,961,970 was unchanged. The four positions were selected after inspecting reads in this window.

Four-SNP mask and rephasing command
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/analysis/crossover_experiments/reproduce.sh#L53-L89
```

### Slide 11 · Spacing and public masks

```text
As a review flag, we marked a switch if its bracket lay within 10 kb of another bracket for the same child and chromosome. We pooled both parents and counted overlapping brackets as zero distance. All 35 flagged baseline switches overlapped both the UCSC segmental duplication track and the ENCODE hg38 blacklist. Three isolated baseline switches overlapped both tracks too.

After quality filtering, the spacing rule flagged all three boundary switches. Each was outside both public tracks. Excluding variants in the two tracks before rerunning WhatsHap still left two boundary switches. Spacing can flag problematic clusters missed by these annotations and can leave isolated calls in annotated regions. The 10 kb rule was explored in one window.

Spacing and overlap calculation
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/analysis/crossover_experiments/summarize.py

Public mask before phasing
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/analysis/crossover_experiments/reproduce.sh#L41-L57
```

### Slide 12 · Read-supported paternal candidate

```text
For chr22:20,959,352–20,961,970, the father had 17 G–A and 17 A–G primary spanning reads. The child had eight recombinant G–G primary spanning reads. The paternal allele pairs agree with the original HiPhase VCF.

The broad candidate in the inheritance review notes spans 20,955,482–20,983,035 and overlaps the WhatsHap bracket. The assembly for that review interval is unverified. The checked-in GRCh38 inheritance map places the paternal E-to-G transition at 20,982,976–20,983,035, with its left marker 21,006 bp beyond the right WhatsHap marker. The maternal label stays I. The pedigree comparison supports a paternal change in the same neighborhood.

Read audit and pedigree comparison
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/src/filter_whatshap_crossovers_by_reads.py
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/trio_dev_data/input/pedigree/CEPH1463.GRCh38.iht.sorted.txt
https://docs.google.com/document/d/1HGBhjqk4cqZaWIkh0GeAaBC_fcZYWRRnYO5cjiifKEY/edit
```

### Slide 13 · Recombinator variant filtering

```text
Recombinator requires site call rate ≥0.95 and accepts records with no failing VCF FILTER. For non-SNP records it screens for REF length ≤3, one ALT and ALT length ≤3. It can also take an exclusion BED. The defaults within each family are called genotypes with GQ ≥20 and DP ≥18.

For unphased records, the transmitting heterozygous parent’s ALT fraction must be >0.31 and <0.69, and the children’s fractions must be >0.31 and ≤0.69 where evaluated. Fully phased records take a separate path before that allele-balance check. Downstream, the caller defaults to 20 informative sites per block and suppresses adjacent short runs of rapid state flipping.

Applying just the DP ≥18 floor to all three members of our quality-filtered WhatsHap input removed the right marker of the supported paternal candidate, where the child had DP 17. That rerun yielded zero switches in the window.

Recombinator defaults, marker selection and block filtering
https://github.com/quinlan-lab/recombinator/blob/a164c2ea5e91debacbe658e85fa38e89ebafad05/recombinator/recombinator.py#L76-L85
https://github.com/quinlan-lab/recombinator/blob/a164c2ea5e91debacbe658e85fa38e89ebafad05/recombinator/recombinator.py#L425-L493
https://github.com/quinlan-lab/recombinator/blob/a164c2ea5e91debacbe658e85fa38e89ebafad05/recombinator/recombinator.py#L32-L49
https://github.com/quinlan-lab/recombinator/blob/a164c2ea5e91debacbe658e85fa38e89ebafad05/recombinator/recombinator.py#L124-L194

WhatsHap depth sensitivity run
https://github.com/quinlan-lab/tapestry/blob/3d5894a68d5a7c2eb406065bcf19392997ec03b9/analysis/crossover_experiments/reproduce.sh#L58-L59
```

### Slide 14 · Platinum Pedigree Inheritance variant filtering

```text
At Tapestry’s pinned revision, gtg-ped-map defaults to QUAL ≥20, DP ≥10 for every family member and a marker-run setting of 10. It skips length-changing alleles and sites with missing family-member alleles. For each sample it also requires depth within one standard deviation of that sample’s chromosome-wide mean, as implemented in the code.

The algorithm tracks informative alleles through the pedigree. It masks a marker if its same-state counts on both sides, each including the marker itself, are below the run setting. This threshold is in markers.

gtg-ped-map defaults, record filter, depth filter and run masking
https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L48-L58
https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L156-L183
https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L217-L236
https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L935-L967
```

## July deck

[WhatsHap phase blocks and recombination](https://docs.google.com/presentation/d/1WNP8aAWJM2DS_gKk1CvkkQO4UIj4-K_WlG8FwQ9KZNU/edit#slide=id.g3eeef4e7a8e_0_10)

### Slide 6 · Finding a crossover from mismatches

```text
The single-crossover heuristic used maternal concordance of 0.57 to predict a terminal minority segment under roughly uniform heterozygous-marker density. About 82% of maternal mismatches fell in the predicted segment. The score ranked blocks by mismatch concentration for follow-up review. The standalone scripts for this example are archived and remain available at the historical revision below.

Historical shell wrapper and Python scoring code
https://github.com/quinlan-lab/tapestry/blob/19df1034b45cf15c957e2218e97da2fd7cee125f/find_clean_single_crossover.sh
https://github.com/quinlan-lab/tapestry/blob/19df1034b45cf15c957e2218e97da2fd7cee125f/src/util/find_clean_single_crossover.py
```
