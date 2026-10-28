# Human review ledger

This ledger records scoped human reviews of specific repository revisions.
Pending entries are a review checklist, not evidence that a review occurred.
Automated checks and AI inspection do not constitute human review.

## How to record a review

Each entry has five fields: Status, Reviewer, Review date, Scope, and Review
questions. An unchecked Status box means pending; replace it with a green ✅
marker when the entire stated scope has been manually reviewed. Replace
individual Scope boxes with ✅ markers as you finish those parts, and supply
the reviewer's name and actual review date (`YYYY-MM-DD`). Use only the green
marker for completed items, without a checkbox. These markers record review
completion, not approval or proof that the code is correct.

Only record completion on the reviewer's explicit statement. If a review finds
issues, record them alongside the relevant review questions. Reading code does
not establish that tests passed or scientific parity was demonstrated. Commit
completed entries to preserve when the review was recorded. For each future
HR-XXX completion requested by Peter McHale, replace the Status and all Scope
boxes with ✅ markers, record the reviewer and actual review date, then automatically commit
and push the update.

A review applies only to the stated revision and scope. To extend it to a newer
revision, inspect the intervening diff and add a dated follow-up entry naming
both revisions and the scope checked. Preserve the original entry. Line numbers
below refer to the pinned revision, not current `main`.

## Brent's generic-workflow changes

- Base revision, before Brent's first commit:
  `700e30c9779c2bef85679e10c3dbdcc1c05b96ec`.
- Target revision, Brent's final commit:
  `e2502e2c78cfd50009d10c301043c5be9e32c4a1`.
- Commit dates: August 14–September 15, 2026; all 18 commits in this range are
  authored by Brent Pedersen.
- Net diff size: 50 files changed, 12,134 additions, and 228 deletions.
- [Net diff](https://github.com/quinlan-lab/tapestry/compare/700e30c9779c2bef85679e10c3dbdcc1c05b96ec..e2502e2c78cfd50009d10c301043c5be9e32c4a1).

Every entry below uses that base and target revision. Later fixes and cleanup
on `main` are outside this review scope. The entries cover selected high-priority
areas, not the entire net diff. Review questions are prompts to investigate,
not human findings or approvals.

### Overall change summary

The net change adds a generic Nextflow workflow downstream of completed PacBio
HiFi human WGS family-WDL runs. It consumes miniwdl `outputs.json`, a six-column
PED, and an indexed GRCh38 reference directly, supporting WDL v3.3.0/v3.3.1 and
autosomal pedigree/model-only processing. Validation generates the internal
manifest and resolved settings before scientific stages run.

The workflow normalizes the joint VCF, runs pinned `gtg` inheritance mapping and
concordance, filters model methylation BEDs, reconciles read-backed and
inheritance-backed phase, and publishes all-CpG tables, BigWigs, QC, provenance,
and a results manifest. Founder phasing processes variants and methylation by
chromosome to reduce peak memory, preserves null-founder outputs when phase
cannot be assigned, and explicitly handles duplicate BigWig intervals.

New offline ancestry and transmission-QC views display inherited founder
haplotypes and methylation summaries. Their analysis introduces parent–child
comparison metrics, mismatch and ambiguity exclusions, and cohort-relative
outlier scores. Shared BED/VCF writers also change, creating review obligations
for legacy callers as well as the generic path.

A pinned container environment, Docker/Apptainer/Slurm profiles, unit and
regression tests, and GitHub CI support execution and reproducibility. The
integration fixture checks workflow completion and complete resume caching.
Real WDL/CEPH scientific parity and target-cluster runtime checks remain release
gates at the target revision. This summary describes the implementation scope;
it does not assert human approval or completion of those gates.

### HR-001: All-CpG construction and biological annotations — Highest scientific priority

- ✅ Review complete
- Reviewer: Peter McHale
- Review date: 2026-10-26
- Scope:
  - ✅ [expand_model_to_all_cpgs.py:208–245 — joins and null handling](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/expand_model_to_all_cpgs.py#L208-L245)
  - ✅ [expand_model_to_all_cpgs.py:118–139 — mismatch proximity](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/expand_model_to_all_cpgs.py#L118-L139)
  - ✅ [expand_model_to_all_cpgs.py:143–191 — variant annotations](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/expand_model_to_all_cpgs.py#L143-L191)
- Review questions: Are reference-only/sample-only CpGs and null measurements
  preserved? Is overlap multiplicity intentional? Are coordinate conversions,
  the strict mismatch-distance boundary, multiallelic handling, and the
  allele-specific flag scientifically appropriate?

### HR-002: VCF normalization and inheritance inference — Highest scientific priority

- ✅ Review complete
- Reviewer: Peter McHale
- Review date: 2026-10-27
- Scope:
  - ✅ [normalize_joint_vcf.py:205–264 — phase metadata and VCF branches](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/normalize_joint_vcf.py#L205-L264)
  - ✅ [run_gtg_inheritance.py:166–220 — mapping and concordance commands](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/run_gtg_inheritance.py#L166-L220)
  - ✅ [ROADMAP.md:46–73 — real WDL/CEPH parity release gate](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/ROADMAP.md#L46-L73)
- Review questions: Are PS/PF removal, sample ordering, and complete-family
  filtering correct? Do mapping and concordance receive their intended VCF
  branches and thresholds? Are incomplete-family sites retained for annotation?
  Has real WDL/CEPH parity evidence been examined? Code review alone does not
  close that release gate.

### HR-003: Founder phasing, empty results, and duplicate policy — High scientific priority

- ✅ Review complete
- Reviewer: Peter McHale
- Review date: 2026-10-27
- Scope:
  - ✅ [phase_meth_to_founder_haps.py:296–320 — chromosome-wise reconciliation](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/phase_meth_to_founder_haps.py#L296-L320)
  - ✅ [phase_meth_to_founder_haps.py:34–62 — empty-map behavior](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/phase_meth_to_founder_haps.py#L34-L62)
  - ✅ [phase_meth_to_founder_haps.py:148–175 — BigWig duplicate handling](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/phase_meth_to_founder_haps.py#L148-L175)
- Review questions: Does batching preserve nested phase-set membership and
  chromosome boundaries? Do empty results represent unassignable founder phase
  correctly? Is collapsing identical BigWig values and rejecting conflicts
  consistent with the intended downstream ambiguity policy?

### HR-004: Transmission QC definitions and exclusions — High scientific priority

- ✅ Review complete
- Reviewer: Peter McHale
- Review date: 2026-10-27
- Scope:
  - ✅ [summarize_transmission_methylation.py:338–370 — pairing and exclusions](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/summarize_transmission_methylation.py#L338-L370)
  - ✅ [summarize_transmission_methylation.py:291–324 — metric formulas](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/summarize_transmission_methylation.py#L291-L324)
  - ✅ [summarize_transmission_methylation.py:129–171 — duplicate merging](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/summarize_transmission_methylation.py#L129-L171)
- Review questions: Are transmitted labels, mismatch exclusions, denominators,
  and specificity sign correct? Should a conflict on either haplotype exclude
  the entire coordinate? Are null merging and numeric tolerance appropriate?

### HR-005: QC scoring and displayed methylation summaries — High interpretation burden

- ✅ Review complete
- Reviewer: Peter McHale
- Review date: 2026-10-27
- Scope:
  - ✅ [plot_transmission_qc.py:70–101 — outlier scoring](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/plot_transmission_qc.py#L70-L101)
  - ✅ [plot_transmission_qc.py:146–182 — genome-wide aggregation](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/plot_transmission_qc.py#L146-L182)
  - ✅ [summarize_haplotype_methylation.py:164–188 — methylation means](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/summarize_haplotype_methylation.py#L164-L188)
- Review questions: Are cohort size, MAD scaling, fallback behavior, and related
  pedigree edges handled appropriately for interpretation? Are chromosome
  eligibility and CpG weighting defensible? Can duplicate rows overweight a
  coordinate in displayed means?

### HR-006: Input identity and Nextflow wiring — High operational priority

- [ ] Review complete
- Reviewer: not yet recorded
- Review date: not yet recorded
- Scope:
  - [ ] [tapestry_validate.py:587–638 — sample/artifact association](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/tapestry_validate.py#L587-L638)
  - [ ] [tapestry_validate.py:852–892 — artifact validation](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/tapestry_validate.py#L852-L892)
  - [ ] [main.nf:774–816 — stage wiring](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/main.nf#L774-L816)
- Review questions: Do array ordering and identity checks associate artifacts
  with the right samples? Are null artifacts handled correctly? Do consuming
  stages enforce record invariants omitted from preflight? Are positional tuple
  indexes and sample-keyed joins correct?

### HR-007: Shared legacy helpers — High regression risk

- [ ] Review complete
- Reviewer: not yet recorded
- Review date: not yet recorded
- Scope:
  - [ ] [util/hap_map.py:26–48 — shared hap-map writer](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/util/hap_map.py#L26-L48)
  - [ ] [hap_map_trio.py:98–117 — trio SNV-count column](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/hap_map_trio.py#L98-L117)
  - [ ] [util/write_data.py:122–135 — methylation writer](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/util/write_data.py#L122-L135)
- Review questions: Are writer/producer column contracts compatible in both
  pedigree and trio paths? After removing explicit shell sorting, do all callers
  provide sorted data? Are output names and indexes compatible with consumers?

### HR-008: Container, provenance, and verification — Focused technical review

- [ ] Review complete
- Reviewer: not yet recorded
- Review date: not yet recorded
- Scope:
  - [ ] [Dockerfile:57–67 — native dependency build](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/Dockerfile#L57-L67)
  - [ ] [Dockerfile:3–5 — UCSC binary pin](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/Dockerfile#L3-L5)
  - [ ] [capture_versions.py:15–18 — recorded checksum constant](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/src/capture_versions.py#L15-L18)
  - [ ] [run_nextflow_e2e.sh:104–112 — resume assertion](https://github.com/quinlan-lab/tapestry/blob/e2502e2c78cfd50009d10c301043c5be9e32c4a1/tests/run_nextflow_e2e.sh#L104-L112)
- Review questions: Do native builds work on the target FIPS-enabled Apptainer
  host? Does recorded provenance match the installed binary? Does the fixture
  exercise all intended stages and caching? Which scientific and cluster-runtime
  gates remain beyond this synthetic test?
