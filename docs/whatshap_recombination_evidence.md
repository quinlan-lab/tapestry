# WhatsHap recombination reports and transmitted haplotypes

Evidence reviewed 2026-10-30 against upstream **WhatsHap v2.8**, commit
`f07423a7ba0c7609e7d4d5c73fb5e4e240057830`. This is
source inspection and an upstream test example, not a new end-to-end run.

`whatshap phase --ped family.ped --recombination-list crossovers.txt ...`
can report putative crossovers and the parental haplotype switches across
them. **The report alone does not identify the transmitted haplotypes in
phase blocks with no reported crossover.** This is a limitation of the
event report, not a claim that WhatsHap lacks internal transmission states.

## Source-code evidence and upstream test

1. **The option requests events.** The
   [v2.8 CLI definition](https://github.com/whatshap/whatshap/blob/f07423a7ba0c7609e7d4d5c73fb5e4e240057830/whatshap/cli/phase.py#L1136-L1137)
   (`add_arguments`, pedigree-phasing options) describes
   `--recombination-list` as “Write putative recombination events to FILE.”
   Its [`write_recombination_list` function](https://github.com/whatshap/whatshap/blob/f07423a7ba0c7609e7d4d5c73fb5e4e240057830/whatshap/cli/phase.py#L992-L1027) writes child ID, chromosome,
   two flanking positions, four parental transmission fields, and
   recombination cost. The fields ending in `1` and `2` refer to the two
   positions: `transmitted_hap_father1`, `transmitted_hap_father2`,
   `transmitted_hap_mother1`, `transmitted_hap_mother2`. The writer [adds one](https://github.com/whatshap/whatshap/blob/f07423a7ba0c7609e7d4d5c73fb5e4e240057830/whatshap/cli/phase.py#L1013-L1023)
   to internal positions, so reported positions are one-based variant
   coordinates bracketing an inferred event, not an exact crossover base.

2. **Only changes produce events.** In
   [v2.8 `find_recombination`](https://github.com/whatshap/whatshap/blob/f07423a7ba0c7609e7d4d5c73fb5e4e240057830/whatshap/pedigree.py#L138-L176),
   positions are grouped by component/block. An event is appended only
   when consecutive inspected transmission states differ. Father states
   are decoded with `% 2`, mother states with `// 2`; these are binary
   haplotype indices. No initial-state or constant-transmission row is
   emitted. Consequently a constant-state block contributes no event,
   regardless of which parental homologs were transmitted.

3. **An upstream test demonstrates the switches.**
   [v2.8 `test_find_recombination`](https://github.com/whatshap/whatshap/blob/f07423a7ba0c7609e7d4d5c73fb5e4e240057830/tests/test_pedigree.py#L53-L78)
   supplies one block with transmission vector `[0, 0, 1, 1, 0]` and expects
   two events: father `0 → 1` between internal positions 5432 and 8307,
   then father `1 → 0` between 9000 and 9500. Mother remains `0 → 0` in
   both events. Thus each row records both parents' flanking states, even
   when only one parent switches. The CLI writer would report positions
   5433/8308 and 9001/9501.

## Upstream GitHub issues and pull requests

Discussions and change history checked 2026-10-30:

| Reference | Evidence and scope |
| --- | --- |
| [Issue #295: infer break point after phasing](https://github.com/whatshap/whatshap/issues/295) | Tobias Marschall explicitly recommends [`--recombination-list`](https://github.com/whatshap/whatshap/issues/295#issuecomment-801200899) for reporting pedigree recombination events. He cautions about false positives in segmental duplications and [suggests plotting the reported inheritance states](https://github.com/whatshap/whatshap/issues/295#issuecomment-804822832) to distinguish rapid switching from persistent switches. A [user example](https://github.com/whatshap/whatshap/issues/295#issuecomment-805020274) includes both parents' before/after transmission columns. |
| [Issue #96: Recombination not detected for two variants](https://github.com/whatshap/whatshap/issues/96) | Reports exactly the boundary behavior visible in the v2.8 source: the first pair of variants is not checked, and blocks containing only two variants are skipped. This corroborates why absence of a report row cannot establish absence of a crossover. |
| [Issue #394: RecombinationEvent sorting TypeError](https://github.com/whatshap/whatshap/issues/394) | Reports a crash with `--recombination-list` in WhatsHap 1.1. Marcel Martin [states that the fix will be in v1.5](https://github.com/whatshap/whatshap/issues/394#issuecomment-1223099047). This is a report-generation bug, not a request for complete transmission output. |
| [PR #512: Heuristic version for (Ped)MEC algorithm](https://github.com/whatshap/whatshap/pull/512) (merged 2024-01-26) | Related algorithm context: the description discusses internal transmission states and predicting recombination as a change of transmission. It warns that the heuristic can miss events in corner cases. It does not establish a complete block-wise transmission export or remove the event-only limitation described above. |

The issues directly support crossover reporting and its caveats; PR #512
provides related algorithm context. The conclusion about missing transmission
labels in event-free blocks rests on the report writer and event-selection
logic inspected above, rather than an explicit statement in these discussions.

## Interpretation for Tapestry

An event-free block could have any constant father/mother transmission
combination; the event list cannot distinguish them. Absence of a row also
does not prove biological absence of a crossover: these are inferred events,
and v2.8 skips blocks with at most two positions and begins comparisons at
index 2 (omitting the first adjacent pair). This follows directly from
the [block-length guard and comparison loop](https://github.com/whatshap/whatshap/blob/f07423a7ba0c7609e7d4d5c73fb5e4e240057830/whatshap/pedigree.py#L157-L160).

Transmission labels must also be interpreted within their phase-block
orientation; the binary indices are not automatically Tapestry's A/B/C/D
labels across disconnected blocks. Tapestry's trio workflow therefore still
needs its parent/child phased-allele comparison to assign parental homologs
in blocks without reported events. The current
[`run-whatshap.sh` invocation](https://github.com/quinlan-lab/tapestry/blob/438d3838b4a7725d5676611358008e3b63b8976f/run-whatshap.sh#L143-L152) does not request this optional report.
