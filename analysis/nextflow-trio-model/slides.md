# Phase-block and trio methylation slide notes

Source: [the second slide deck](https://docs.google.com/presentation/d/1fMsmfMTXjtrsqQmHAlmu7VouI0DEg1I7fLOuYRM2h8c/edit).
These notes map its six slides to reproducible evidence. Paths below are relative
to the repository root. See [README.md](README.md) for commands and coordinate
conventions. The linked deck is not modified by this analysis.

## Slide 1 · TODOs from last time (261103)

The preceding [crossover-filter analysis](../crossover_experiments/README.md)
examined variant quality, public masks, targeted SNP removal, depth, and read
evidence. The current experiment applies crossover-aware bit-vector matching
and model methylation in the Nextflow trio path. The 1 Mb fixture is the scope
of the evidence; manuscript integration is a subsequent task.

## Slide 2 · WhatsHap phase blocks

Figure: `trio_dev_data/output/nextflow-trio-model/phase-blocks/phase-blocks.png`.

WhatsHap stats reports one phase block for each member. The child's block spans
20,500,019–21,349,175 with 447 variants; the father's spans
20,500,019–21,485,412 with 949; the mother's spans 20,500,927–21,485,412
with 1,065. Coordinates in the block lists are one-based inclusive. Phase blocks
can span transmission changes, so matching child and parent across an entire
block can obscure a crossover. Block-list variant counts and informative SNV
counts used by the matcher measure different things.

Evidence: `phase-blocks/NA12879.blocks.tsv`, `NA12877.blocks.tsv`, and
`NA12878.blocks.tsv` within the run. Commands: [reproduce.sh](reproduce.sh).
Figure generation: [summarize.py](summarize.py).

## Slide 3 · Likely false and true crossovers

Figure: `trio_dev_data/output/nextflow-trio-model/crossover-brackets.png`.

The quality-filtered input yields four parent-specific switches. The bracket
20,502,763–20,503,731 switches both parents; a second paternal bracket is
20,503,755–20,510,653. Their nearest gaps are zero, zero, and 24 bp, so all
three events are flagged at 10 kb. The isolated paternal bracket
20,959,352–20,961,970 is 448,699 bp from its nearest neighbor. The first
analysis supplies the read-audit evidence for suspect cluster anchors and the
supported isolated candidate. Distance alone does not establish truth.

Evidence: `per_chrom_phased/chr22.crossovers.tsv` and the raw recombination
report. Caller/screen: [report_whatshap_crossovers.py](../../src/report_whatshap_crossovers.py).

## Slide 4 · Matching bitvectors either side of “true” crossover

Figure: `trio_dev_data/output/nextflow-trio-model/resolved-intervals.png`.

Matching is restricted to resolved portions of each child/parent phase-set pair.
The paternal A interval before the isolated bracket has 620 informative SNVs;
the B interval after it has 189. Both have concordance 1.0 and zero mismatches.
The maternal D interval has 918 informative SNVs and crosses the paternal
bracket without interruption. The bracket remains unresolved for paternal
assignment. The flagged left-edge cluster is unresolved for both parents,
including the gap between brackets. The very short blocks before the cluster
have only three paternal SNVs and one maternal SNV, despite perfect concordance.

Evidence: `resolved-intervals.tsv`, the two parental map tables, and separate
mismatch BEDs. Matching: [hap_map_trio.py](../../src/hap_map_trio.py).
Exclusion spans: [report_whatshap_crossovers.py](../../src/report_whatshap_crossovers.py).

## Slide 5 · Phased methylation in the child

Figure: `trio_dev_data/output/nextflow-trio-model/parent_phased_methylation.png`.

The all-CpG output has 19,206 rows. Child paternal, maternal, and combined model
measurements are non-null on 5,931, 5,073, and 16,759 rows respectively.
Stored fractions are plotted as percentages. Missing calls remain null.
Combined measurements include untagged reads; they do not require a resolved
parental assignment. Gray shading in each child haplotype lane marks unresolved
assignment for that parent; the combined lane is not masked this way.

Evidence: the final all-CpG BED, `trio.model.qc.json`, and `summary.json`.
Assignment/publication: [trio_model_methylation.py](../../src/trio_model_methylation.py)
and [trio_phasing.nf](../../workflows/trio_phasing.nf).
Scientific checks: [check_trio_model_fixture.py](../../tests/check_trio_model_fixture.py).

## Slide 6 · Next steps

Whole-chr22 and genome-wide evaluation, manuscript integration, and downstream
de novo methylation/compound-heterozygosity tools are proposals in the deck.
They are not demonstrated by this fixture. No additional experiment or
performance measurement is implied by this slide.
