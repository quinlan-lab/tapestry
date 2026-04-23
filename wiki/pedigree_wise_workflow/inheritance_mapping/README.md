# Inheritance mapping (vendored from upstream)

This section documents Step 1B of the pedigree-wise workflow
(`build-iht-based-haplotype-map-and-phase-variants.sh`), which is
implemented in the Rust binaries `gtg-ped-map` (`map_builder.rs`) and
`gtg-concordance` (`gtg_concordance.rs`).

The three sub-pages below are **vendored verbatim** from the upstream
[Platinum-Pedigree-Inheritance wiki](https://github.com/petermchale/Platinum-Pedigree-Inheritance/tree/main/wiki),
synced from commit
[`7448e5e946adbc7969ad5fd5e0730d7cace23a8d`](https://github.com/petermchale/Platinum-Pedigree-Inheritance/tree/7448e5e946adbc7969ad5fd5e0730d7cace23a8d/wiki).
This is the user's fork of the
[Platinum-Pedigree-Consortium](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance)
repository; the fork is the source of truth until the user's open PR
to the consortium repo is merged. Once merged, this README should be
updated to point at the consortium repo as the upstream.

## Pages

1. [Structural haplotype mapping in a nuclear family](nuclear_family/nuclear_family.md)
   — `gtg-ped-map`'s algorithm derived end-to-end on a two-founder,
   three-child pedigree with a single paternal recombinant. Read this
   first; the other two pages assume its algorithmic vocabulary.
2. [Structural haplotype mapping across three generations](three_generations/three_generations.md)
   — extends the nuclear-family pedigree by adding an outside-marriage
   spouse and two grandchildren; covers the depth-ordered ancestor-first
   walk and the ancestral-vs-de-novo crossover diagnostic.
3. [Phasing alleles consistently with the haplotype map](concordance/concordance.md)
   — `gtg-concordance`'s `2^F` founder-phase enumeration, which
   re-reads the VCF inside each IHT block and phases every variant
   using the block's letter map.

## Re-sync procedure

To pull in upstream changes:

```bash
# 1. Choose the upstream commit you want to vendor (here, current main of the fork).
UPSTREAM_REPO=petermchale/Platinum-Pedigree-Inheritance
UPSTREAM_SHA=$(gh api repos/${UPSTREAM_REPO}/commits/main --jq '.sha')

# 2. Clone the upstream at that SHA into a tmp dir.
TMP=$(mktemp -d)
git clone --depth 1 --branch main https://github.com/${UPSTREAM_REPO}.git "${TMP}/ppi"

# 3. Verify the SHA matches what was resolved in step 1.
git -C "${TMP}/ppi" rev-parse HEAD   # should equal $UPSTREAM_SHA

# 4. Rsync the three vendored folders into this directory.
rsync -a --delete \
  "${TMP}/ppi/wiki/nuclear_family" \
  "${TMP}/ppi/wiki/three_generations" \
  "${TMP}/ppi/wiki/concordance" \
  wiki/pedigree_wise_workflow/inheritance_mapping/

# 5. Re-apply the one local edit: the [wiki](../index.md) link at the top
#    of each vendored page must point at the tapestry top-level wiki, not
#    the upstream one.
for f in wiki/pedigree_wise_workflow/inheritance_mapping/*/*.md; do
  sed -i '' 's|\[wiki\](\.\./index\.md)|[wiki](../../../index.md)|g' "$f"
done

# 6. Update the pinned SHA at the top of this README to ${UPSTREAM_SHA}.

# 7. Clean up.
rm -rf "${TMP}"
```

The upstream `wiki/generate_wiki.py` is idempotent, so this procedure
produces a byte-stable result modulo the deliberate `[wiki]` link
rewrite in step 5. Rust permalinks inside the vendored pages continue
to point at upstream commit SHAs — leave them alone; they live in the
upstream repo, not in tapestry.

## Notes

- The vendored pages link to each other via relative sibling paths
  (e.g. `../nuclear_family/nuclear_family.md` from
  `three_generations/three_generations.md`). Those links continue to
  resolve correctly inside `inheritance_mapping/` because the three
  sub-folders retain their upstream sibling structure.
- The only link rewritten on vendoring is the `[wiki](...)` link at
  the top of each page, which originally pointed at the upstream
  `wiki/index.md` (`../index.md`) and now points at the tapestry
  top-level catalog (`../../../index.md`).
