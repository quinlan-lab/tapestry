# Trio-wise workflow (mostly stubbed)

This section will document tapestry's trio-wise workflow — pedMEC /
whatshap phasing across a parent–parent–child trio, with each
methylation measurement labelled by one of the parents' haplotypes.
Most of the section is *not yet built*. For the time being, the
canonical reference is the **Trio-wise workflow** section of the
top-level [`README.md`](../../README.md).

For motivation — *why* phase methylation in the first place — see
[the shared motivation page](../motivation/motivation.md).

## Planned page structure

```
trio_wise_workflow/
  index.md                   # this page
  pedmec_phasing/            # Step 1 — run-whatshap.sh                  (TODO)
  parent_haplotype_phasing/  # Step 3 — phase_meth_to_parent_haps.py     (TODO)
  all_cpg_expansion_trio/    # Step 4 — expand_to_all_cpgs.trio.sh       (TODO)
```

The trio-wise side will be much shorter than the pedigree-wise side
because the bit-vector concordance machinery is established in general
form on the [pedigree-wise founder-phased-methylation
page](../pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md);
the trio-wise version only has to explain the two differences:
(a) the phase source is pedMEC / whatshap rather than
`gtg-concordance`, and (b) the kid's haplotypes are labelled by
parental letters (A/B in dad, C/D in mom) rather than by
founder-of-the-pedigree letters. A/B/C/D are *fixed as dad's hap1/hap2
and mom's hap1/hap2* — they are not defined as transmitted vs
non-transmitted.
