# tapestry

Trace PacBio HiFi CpG methylation to parental and founder haplotypes.
Both Nextflow workflows use one published container image.

| Workflow | Inputs | Command |
|---|---|---|
| [Trio](#trio-workflow) | Trio PED, joint VCF, three BAMs, reference | `nextflow run . -entry trio -profile docker ...` |
| [Pedigree](#pedigree-workflow) | Completed family-WDL `outputs.json`, PED, GRCh38 reference | `nextflow run . -entry pedigree -profile docker ...` |

**Migrating existing commands:** WDL-output runs now require `-entry pedigree`.
Trio runs now require `-entry trio`; `--container` is the image override for both modes.
Changes to entry-point names or workflow module layout may invalidate older
task caches; new runs resume normally with `-resume`.

## Setup: one image for both workflows

Run from a repository checkout with Nextflow **24.04.2 or later** and Docker.
Nextflow runs on the host; Docker automatically pulls
`ghcr.io/quinlan-lab/tapestry:dev` and supplies the scientific tools, including
WhatsHap **2.8**, samtools **1.21**, and pb-CpG-tools **3.0.0**. No local Python
setup or image build is required. Images target **Linux amd64**; Apple Silicon
requires amd64 emulation. Docker tasks use the launching user's UID/GID.
Nextflow execution reports go to `.nextflow-reports/` in the launch directory;
`-with-*` options can redirect them.

For **CHPC/Slurm**, replace `-profile docker` with `-profile slurm,apptainer`.
Configure your site's account, queue, and resources in a local Nextflow config;
see [nextflow.config](nextflow.config). To use a SIF already on shared storage:

```bash
apptainer pull tapestry.sif docker://ghcr.io/quinlan-lab/tapestry:dev
# Add --container /absolute/path/to/tapestry.sif to either workflow command.
```

The [GHCR package](https://github.com/orgs/quinlan-lab/packages/container/package/tapestry)
provides the moving `dev` tag and `sha-<full-git-commit>` tags. CI publishes only
after checks pass on main-branch pushes or manual runs, recording image digests
in the publishing job summary. Pin an exact image for reproducible runs with
`--container ghcr.io/quinlan-lab/tapestry@sha256:<digest>`. `dev` is a development
image, not a stable release. Cluster verification and release gates are tracked
in [ROADMAP.md](ROADMAP.md).

<details>
<summary><h3>Build locally or use installed tools</h3></summary>

```bash
docker build --platform linux/amd64 -t tapestry:dev .
# Add --container tapestry:dev to either Docker workflow command.
```

The [Dockerfile](Dockerfile) and [requirements.txt](requirements.txt) pin the
shared runtime dependencies; samtools and pb-CpG-tools archives are SHA-256
checked. Without a container profile, provide these tools and Python dependencies
on each worker's `PATH`. For local notebooks, figures, IGV serving, or legacy
scripts, see the [dependency guide](docs/dependencies.md).

</details>

## Trio workflow

Supply a six-column PED defining one trio, an indexed joint VCF, three indexed
BAMs with read-group `SM` tags matching the PED IDs, and an indexed reference.
BAMs must retain the `MM`/`ML` methylation tags. Keep `.bai`/`.csi` BAM indexes,
`.tbi`/`.csi` VCF indexes, and the reference `.fai` beside their files; a BGZF
FASTA also needs its `.gzi`.

```bash
nextflow run . -entry trio -profile docker \
  --ped /path/to/trio.ped \
  --joint-vcf /path/to/joint.vcf.gz \
  --reference-fasta /path/to/GRCh38.fa \
  --child-bam /path/to/child.bam \
  --father-bam /path/to/father.bam \
  --mother-bam /path/to/mother.bam \
  --chromosomes chr22 \
  --outdir /path/to/trio-results \
  -resume
```

The workflow runs WhatsHap pedMEC, screens crossovers, matches parental
haplotypes, haplotags all three BAMs, calls model methylation with pb-CpG-tools,
and expands to reference and sample CpGs. Trio labels are child hap1=paternal,
hap2=maternal; father hap1/hap2=A/B; mother hap1/hap2=C/D.

| Option | Default / purpose |
|---|---|
| `--chromosomes chr1,chr2,...` | Required chromosome list; processed in parallel |
| `--region chr22:20500000-21500000` | Optional one-based inclusive interval; requires a single matching chromosome |
| `--min-coverage` | `10` |
| `--trio-min-mapq` | `1` |
| `--crossover-distance-bp` | `10000`; nearby brackets form unresolved clusters |

### Trio results

Outputs are chromosome-local under `--outdir`; chromosomes are not concatenated.

| Directory | Contents |
|---|---|
| `per_chrom_phased/` | Indexed phased VCF, raw recombination report, screened crossover TSV/BED |
| `per_chrom_hap_maps/CHROM/` | Paternal/maternal maps, mismatch BEDs/VCFs, headers/indexes, full-precision `*.hap-map.tsv` |
| `per_chrom_model/CHROM/` | Indexed haplotagged BAMs and pb-CpG-tools caller outputs |
| `per_chrom_methylation/CHROM/` | `trio.dna-methylation.all_cpgs.sorted.bed.gz` and `.tbi`, six haplotype and three combined model BigWigs, `trio.model.qc.json` |

Unresolved sides retain null haplotype measurements and coverage; combined
measurements remain available. Uncovered reference CpGs are retained. Trio
Nextflow outputs omit count-methylation columns. Final BED intervals span the
CG dinucleotide; BigWigs span the cytosine base and use a 0–1 scale.

<details>
<summary><h3>Trio output column definitions</h3></summary>

BED coordinates are zero-based, half-open. Read each file's header for its
emitted columns and metadata. This reference also includes legacy shell
count-methylation fields, which the Nextflow trio table omits.

### Output format of trio workflow

Legacy shell BEDs include a header containing run metadata and column headings. The `##source=` line records the script path followed by `with args` and a dictionary of the command-line arguments passed to `src/expand_to_all_cpgs_trio.py`. The final header line lists the column names in the form "#col1 col2 col3...".

Column | Definition
:--- | :---
chrom | chromosome on which CpG site is found
start_cpg | start coordinate of CpG dinucleotide
end_cpg | end coordinate of CpG dinucleotide
total_read_count_kid | total number of HiFi reads overlapping CpG site in the kid
methylation_level_kid_count | count-based unphased methylation level in the kid
methylation_level_kid_model | model-based unphased methylation level in the kid
total_read_count_dad | total number of HiFi reads overlapping CpG site in dad
methylation_level_dad_count | count-based unphased methylation level in dad
methylation_level_dad_model | model-based unphased methylation level in dad
total_read_count_mom | total number of HiFi reads overlapping CpG site in mom
methylation_level_mom_count | count-based unphased methylation level in mom
methylation_level_mom_model | model-based unphased methylation level in mom
total_read_count_kid_pat | number of reads from the kid's paternal haplotype
methylation_level_kid_pat_count | count-based methylation level on the kid's paternal haplotype
total_read_count_kid_mat | number of reads from the kid's maternal haplotype
methylation_level_kid_mat_count | count-based methylation level on the kid's maternal haplotype
methylation_level_kid_pat_model | model-based methylation level on the kid's paternal haplotype
methylation_level_kid_mat_model | model-based methylation level on the kid's maternal haplotype
total_read_count_dad_A | number of reads from dad's haplotype A (hap1)
methylation_level_dad_A_count | count-based methylation level on dad's haplotype A
total_read_count_dad_B | number of reads from dad's haplotype B (hap2)
methylation_level_dad_B_count | count-based methylation level on dad's haplotype B
methylation_level_dad_A_model | model-based methylation level on dad's haplotype A
methylation_level_dad_B_model | model-based methylation level on dad's haplotype B
total_read_count_mom_C | number of reads from mom's haplotype C (hap1)
methylation_level_mom_C_count | count-based methylation level on mom's haplotype C
total_read_count_mom_D | number of reads from mom's haplotype D (hap2)
methylation_level_mom_D_count | count-based methylation level on mom's haplotype D
methylation_level_mom_C_model | model-based methylation level on mom's haplotype C
methylation_level_mom_D_model | model-based methylation level on mom's haplotype D
start_hap_map_block_pat | start of the paternal hap-map block
end_hap_map_block_pat | end of the paternal hap-map block
paternal_haplotype | which of dad's haplotypes (A or B) corresponds to the kid's paternal haplotype in this block
paternal_concordance | bit-vector concordance between the kid's paternal allele sequence and dad's assigned haplotype (A or B) over the het SNVs in the paternal hap-map block. Both vectors come from the same `whatshap phase --ped` run. The assigned haplotype is defined as whichever of dad's two haplotypes gives the higher concordance, so the value lies in `[0.5, 1.0]` by construction
num_het_SNVs_in_dad | number of heterozygous sites in dad within the paternal hap-map block
start_hap_map_block_mat | start of the maternal hap-map block
end_hap_map_block_mat | end of the maternal hap-map block
maternal_haplotype | which of mom's haplotypes (C or D) corresponds to the kid's maternal haplotype in this block
maternal_concordance | bit-vector concordance between the kid's maternal allele sequence and mom's assigned haplotype (C or D) over the het SNVs in the maternal hap-map block. Both vectors come from the same `whatshap phase --ped` run. The assigned haplotype is defined as whichever of mom's two haplotypes gives the higher concordance, so the value lies in `[0.5, 1.0]` by construction
num_het_SNVs_in_mom | number of heterozygous sites in mom within the maternal hap-map block
cpg_is_within_50bp_of_mismatch_site_pat | is the CpG site within 50bp of a het SNV in dad at which the kid's paternal allele disagrees with dad's assigned haplotype (a within-pedMEC bit-vector mismatch on the paternal side)
cpg_is_within_50bp_of_mismatch_site_mat | is the CpG site within 50bp of a het SNV in mom at which the kid's maternal allele disagrees with mom's assigned haplotype (a within-pedMEC bit-vector mismatch on the maternal side)
cpg_overlaps_at_least_one_snv | does the CpG dinucleotide overlap an SNV?
snv_genotypes_kid | are the SNVs (if any) that overlap the CpG dinucleotide `hom` or `het` in the kid?
cpg_is_allele_specific_kid | does the CpG dinucleotide appear in the reads of only one haplotype in the kid?
snv_genotypes_dad | are the SNVs (if any) that overlap the CpG dinucleotide `hom` or `het` in dad?
cpg_is_allele_specific_dad | does the CpG dinucleotide appear in the reads of only one haplotype in dad?
snv_genotypes_mom | are the SNVs (if any) that overlap the CpG dinucleotide `hom` or `het` in mom?
cpg_is_allele_specific_mom | does the CpG dinucleotide appear in the reads of only one haplotype in mom?

</details>

<details>
<summary><h3>Trio phasing, crossover rules, and output details</h3></summary>

For each chromosome, `PHASE_TRIO_CHROMOSOME` subsets the joint VCF, applies the
provisional filter (PASS or `.`, QUAL ≥20, called genotypes with GQ ≥20 and
DP ≥10 in all three samples), clears existing phase, and runs pedMEC. It
publishes `<chrom>.phased.vcf.gz`, `<chrom>.phased.vcf.gz.tbi`, and
WhatsHap’s raw `<chrom>.recombination.txt` report under
`<outdir>/per_chrom_phased/`.

`SCREEN_TRIO_CROSSOVERS` reads the raw report and emits one record per parent
whose transmitted haplotype switches between the two variant markers reported
by WhatsHap. Each record spans both flanking markers. If both parents switch
in one raw row, each parent gets its own record for that interval.

For each record, the process finds the smallest gap to another inferred
crossover for the same child and chromosome, considering both parents. For
disjoint intervals, the gap is the later interval’s left marker position minus
the earlier interval’s right marker position. Intervals that overlap or share
a marker have a gap of zero. The process sets `clustered` to 1 when the smallest
gap is at most `--crossover-distance-bp` and to 0 otherwise. The cutoff is a
non-negative integer that defaults to 10,000 bp. See the
[chr22 experiment](analysis/crossover_experiments/README.md).

The process publishes all inferred crossovers with their `clustered` flags in
`<chrom>.crossovers.tsv` and `<chrom>.crossovers.bed` alongside the raw report.
The TSV records the two one-based marker positions, the parent, the transmitted
haplotype at each marker, the WhatsHap recombination cost, and `nearest_gap_bp`.
The BED spans the same markers in zero-based half-open coordinates and puts
the event ID in column 4, the gap in column 10, and the flag in column 11.
When a crossover has no neighbor, the gap is blank in TSV and `.` in BED.
Configure process resources with the `trio_crossover_screening` label.

After screening, `MATCH_TRIO_BIT_VECTORS` pairs each chromosome's phased VCF
with its screened TSV, generates the trio's WhatsHap block lists, and runs
the crossover-aware matcher with the same distance cutoff. It publishes the
following files under `<outdir>/per_chrom_hap_maps/<chrom>/` for each parent.

- `<child>.hap-map-blocks.<paternal|maternal>.bed.gz` with its tabix index and
  BED column header file. The fourth BED field contains the parental haplotype,
  concordance, and informative SNV count.
- `<child>.bit-vector-sites-mismatches.<paternal|maternal>.bed` with its column
  header file, plus an indexed VCF for visualization. A parent with zero
  mismatches has an empty BED and a header-only VCF.

`CALL_TRIO_MODEL` strips pre-existing HP/PS tags, retains MM/ML, haplotags
all three BAMs against the chromosome's pedMEC VCF, indexes them, and runs
pb-CpG-tools with `--pileup-mode model`, `--min-coverage 10` and `--min-mapq 1`.
Override the latter two with `--min-coverage` and `--trio-min-mapq`.
Resources come from the `trio_model_calling` and `trio_model_expansion` labels.

`EXPAND_TRIO_MODEL` consumes the screened paternal/maternal maps, retaining
child hap1=paternal, hap2=maternal; father hap1/hap2=A/B; mother hap1/hap2=C/D.
Each unresolved side has null haplotype measurements and coverage; combined
model measurements remain available. Original joint-call SNVs, including those
excluded from phasing, annotate allele-specific CpGs. No count measurement is
created. Reference CpGs and sample-only combined CpGs are retained, including
uncovered sites. A region restricts final CpGs by cytosine coordinate; annotation
also includes the following G base. As in the legacy trio table, final intervals
span the CG dinucleotide; BigWigs span the cytosine base and use a 0–1 scale.

</details>

<details>
<summary><h3>Run the checked-in chr22 example</h3></summary>

The fixture covers GRCh38 chr22:20,500,000–21,500,000. Download its reference
once, then run the trio workflow:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" -w /work ghcr.io/quinlan-lab/tapestry:dev \
  bash ./trio_dev_data_get_ref.sh trio_dev_data
fixture=trio_dev_data/input
nextflow run . -entry trio -profile docker \
  --ped "$fixture/trio.ped" \
  --joint-vcf "$fixture/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz" \
  --reference-fasta "$fixture/chr22.fa" \
  --child-bam "$fixture/NA12879.GRCh38.haplotagged.bam" \
  --father-bam "$fixture/NA12877.GRCh38.haplotagged.bam" \
  --mother-bam "$fixture/NA12878.GRCh38.haplotagged.bam" \
  --chromosomes chr22 --region chr22:20500000-21500000 \
  --outdir trio_dev_data/output/nextflow-trio \
  -resume
```

Use the same pinned runtime when comparing phasing results. WhatsHap selects
different reads in the Linux/Python 3.11 image and the original macOS/Python
3.13 host fixture: two versus three raw crossover rows. Both retain the isolated
paternal switch at 20,959,352–20,961,970. Tests keep strict expectations for each
runtime. See the [crossover experiments](analysis/crossover_experiments/README.md)
and [trio methylation analysis](analysis/nextflow-trio-model/README.md).

</details>

## Pedigree workflow

Use completed miniwdl outputs from the PacBio HiFi human WGS `family` WDL
**v3.3.0 or v3.3.1**, a six-column PED, and an indexed **GRCh38** FASTA.
Methylation calling has already run upstream. This workflow consumes the WDL's
HiPhase and model methylation outputs,
reconciles them with pedigree inheritance, and publishes founder-phased
methylation for `chr1`–`chr22`. PED members with both parents present are
selected by default.

```bash
nextflow run . -entry pedigree -profile docker \
  --outputs-json /path/to/miniwdl-run/outputs.json \
  --ped /path/to/family.ped \
  --reference-fasta /path/to/GRCh38.fa \
  --outdir /path/to/pedigree-results \
  -resume
```

### Pedigree options and validation

| Option | Default / purpose |
|---|---|
| `--samples CHILD1,CHILD2` | All PED members with both parents present |
| `--regions chr1,chr2` | All GRCh38 autosomes |
| `--reference-index PATH` | `<reference-fasta>.fai` |
| `--min-coverage` | `10` |
| `--mismatch-window-bp` | `50` |
| `--bigwig false` | Disable BigWigs; enabled by default |
| `--map-min-qual`, `--map-min-depth`, `--min-run-markers` | `20`, `10`, `10` |
| `--concordance-min-qual`, `--concordance-min-depth` | `20`, `5` |
| `--qc-discordance-threshold`, `--qc-min-paired-cpgs` | `0.4`, `100` |

Every pedigree run validates the PED, sample sets, WDL release, reference,
artifact headers/indexes, and output collisions first. Settings are printed and saved
in `pipeline_info/resolved-run.json`; canonical inputs and the validation report
are saved alongside it. Users do not prepare an intermediate manifest.

### Pedigree results

| Path under `--outdir` | Contents |
|---|---|
| `samples/` | Per-sample all-CpG BEDs, indexes, BigWigs, and QC |
| `reference/`, `inheritance/` | Reference CpGs and inheritance-phasing outputs |
| `visualizations/` | Offline haplotype ancestry and transmission QC views |
| `pipeline_info/` | Validated inputs, settings, and provenance |
| `results-manifest.json` | Published relative paths, indexes, and sample statuses |

The pedigree workflow uses model methylation only; count mode is disabled.
Completion prints the manifest path and sample status counts. Processing is
chromosome-wise to bound memory. BigWigs collapse repeated intervals only when
their non-null values agree; conflicts fail.

<details>
<summary><h3>Pedigree output column definitions</h3></summary>

BED coordinates are zero-based, half-open. Read each file's header for its
emitted columns and metadata. This reference also covers legacy shell outputs.
The model-only Nextflow pedigree table retains count fields as nulls and names
its configurable mismatch-proximity flag `cpg_is_within_mismatch_window`.

### Output format of pedigree workflow

The BED header records the script and arguments in `##source=` and lists column
names on the final `#chrom` line. BED coordinates are zero-based, half-open.

Definitions of column headings:

Column | Definition
:--- | :---
chrom | chromosome on which CpG site is found
start_cpg | start coordinate of CpG dinucleotide
end_cpg | end coordinate of CpG dinucleotide
total_read_count | total number of HiFi reads overlapping CpG site
methylation_level_count | methylation level: fraction of reads that are methylated (inferred from sequencing kinetics) at given CpG site, c.f., https://github.com/PacificBiosciences/pb-CpG-tools?tab=readme-ov-file#output-modes-and-option-details
methylation_level_model | methylation level: neural network maps sequencing kinetics at multiple CpG sites to methylation probability at a central CpG site, c.f., "Methods" section of https://www.pacb.com/wp-content/uploads/poster_saunders.pdf
start_hap_map_block | start of hap-map block (a genomic interval in which parental phasing of reads can be performed)
end_hap_map_block | end of hap-map block
haplotype_concordance_in_hap_map_block | the degree of similarity (on a scale of 0 to 1, with 1 representing perfect similarity) of the bit vectors: sequences of alleles observed at heterozygous sites on haplotypes inferred by read-backed vs inheritance-backed phasing
num_het_SNVs_in_hap_map_block | number of heterozygous sites in the given hap-map block
total_read_count_pat | number of reads from the (deduced) paternal haplotype
total_read_count_mat | number of reads from the (deduced) maternal haplotype
founder_haplotype_pat | a label indicating which of the haplotypes in a founder of the pedigree the paternal haplotype in the given individual was descended from
founder_haplotype_mat | a similar label corresponding to the maternal haplotype of the sample
methylation_level_pat_count | count-based methylation level on the paternal haplotype (fraction of paternal reads that are methylated)
methylation_level_mat_count | count-based methylation level on the maternal haplotype (fraction of maternal reads that are methylated)
methylation_level_pat_model | model-based methylation level on the paternal haplotype (deep-learning estimated methylation using only the paternal reads)
methylation_level_mat_model | model-based methylation level on the maternal haplotype (deep-learning estimated methylation using only the maternal reads)
cpg_is_within_50bp_of_mismatch_site | is the given CpG site within 50bp of a heterozygous site at which the read-based and inheritance-based bit vectors (allele sequences corresponding to each haplotype) are mismatched
cpg_overlaps_at_least_one_snv | does the CpG dinucleotide overlap an SNV?
snv_genotypes | are the SNVs (if any) that overlap CpG dinucleotide `hom` or `het`?
cpg_is_allele_specific | does the CpG dinucleotide appear in the reads corresponding to one haplotype, but not the reads of the other haplotype?

</details>

<details>
<summary><h3>Explore ancestry and transmission QC</h3></summary>

Open `<outdir>/visualizations/haplotype-ancestry/index.html` locally. Select a
sample and click an inheritance block to trace its founder haplotype through the
pedigree and compare mean methylation on the two haplotypes. Hover for counts
and QC details. Blank spans indicate missing inheritance blocks;
`inheritance map only` indicates no downstream methylation output. Neither
label describes sequencing coverage. Copy the entire `haplotype-ancestry/`
directory to move the offline visualization.

The **Transmission QC** link compares parent–child methylation, separating
completeness from concordance. Select both members for downstream processing to
obtain a concordance estimate. Click a pair for chromosome details and hover for
metrics and excluded-CpG counts. Colors show cohort-relative outlier severity;
scores stay blank when fewer than four pairs are available.

CpGs are paired by exact coordinates and matching transmitted founder labels,
excluding mismatch windows and ambiguous labels or values. Callable fraction
uses label-matched CpGs after exclusions as its denominator; inherited
specificity additionally requires a measurement on the parent's non-transmitted
haplotype. Ambiguous coordinates remain in shared-CpG counts but are excluded
from eligible and paired counts. These comparisons alone do not establish
inheritance of methylation state.

</details>

## Further reading

- [Supporting documentation](docs/README.md): legacy shell workflows, the
  pedigree diagram, WhatsHap evidence, and human review records.
- [Wiki](wiki/): scientific walkthroughs.
- [Analysis index](analysis/README.md): reproducible experiments and figures.
- [Roadmap](ROADMAP.md) and [completed work](completed.md): release gates,
  decisions, and verification evidence.
- [Developer tests](tests/README.md): container and workflow regression checks.
