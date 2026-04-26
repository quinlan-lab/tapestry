# Structural haplotype mapping across three generations

This page is part of the [wiki](../../../README.md) and extends the
[nuclear-family walkthrough](../nuclear_family/nuclear_family.md) by
marrying **Kid3** — the paternal recombinant from that page — to an
outside-marriage founder **Spouse** and adding two grandchildren
**GK1** and **GK2**. The point is not to re-run the per-site
mechanics from §3-4 of the nuclear-family page, which generalise
verbatim and are not repeated, but to elaborate two things that only
become visible once a third generation is in the pedigree:

1. `gtg-ped-map` does not have a separate "G2→G3 pass" bolted on top
   of a G1→G2 pass. The *one* operation on the nuclear-family page
   that really is a walk over the pedigree —
   [`track_alleles_through_pedigree`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L295), which
   performs [Step 1](../nuclear_family/nuclear_family.md#3-informative-site-detection-founder-letter-tagging-and-haplotype-inference-within-a-linkage-block)
   (flag, per (parent, spouse, child) triple, the sites at which the
   parent is heterozygous and the spouse is homozygous — the sites
   where the parent's two alleles are distinguishable in the
   children) and [Step 2](../nuclear_family/nuclear_family.md#3-informative-site-detection-founder-letter-tagging-and-haplotype-inference-within-a-linkage-block)
   (at each such site, write one of the parent's two letters onto
   whichever child-homolog carries the parent's unique allele) —
   is **a single ancestor-first walk** that considers G1, G2 and G3
   together, iterating over (parent, spouse, child) triples in
   **depth order** — an ordering that processes founders first,
   then their children, then their grandchildren, so that a
   non-founder parent's letters are always in place before she is
   asked to play the parent role in a later triple (§2 gives the
   precise definition and this pedigree's depth assignments). For
   this pedigree the walk visits two triples back-to-back, and by
   the time it reaches the second triple the G2 parent's letters
   have already been written into her slot pair — so they serve as
   "parent letters" for the G3 sub-problem without any special
   casing.

   It is worth emphasising that the rest of the nuclear-family
   pipeline is **not** part of this walk. Sibship backfilling
   ([`backfill_sibs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L804)) runs once per VCF
   record at the whole-family level, right after the walk, consuming
   its output. The flip pass
   ([`perform_flips_in_place`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L702)), block collapse
   ([`collapse_identical_iht`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L385)) and gap-fill
   ([`fill_missing_values`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L617) then
   [`fill_missing_values_by_neighbor`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L540)) are
   run later still — in that execution order, with
   `perform_flips_in_place` invoked three times in total, once
   before collapse and again after each gap-fill step — as
   *global passes over the pooled grid of every record × every
   individual*. None of them iterate triples or depths at all. So
   "single ancestor-first walk" describes the letter-assignment step
   specifically — not the full nuclear-family pipeline — and the
   rest of the post-processing simply sees a wider grid (more
   individuals, one more generation) in the G3 case than it did in
   the nuclear-family case.
2. Two qualitatively different kinds of crossover can appear in a
   grandchild's letter trace: an **ancestral** crossover, inherited
   unchanged from one of Kid3's homologs and therefore shared across
   every G3 descendant that inherits the affected segment, versus a
   **de novo** crossover, introduced in Kid3's own meiosis and
   therefore visible in only one grandchild. Telling them apart
   requires two grandchildren, because the diagnostic signature is
   that the ancestral transition appears in *both* while the de novo
   one appears in only *one*.

All line numbers refer to commit `7bb9f88`. As in the other
walkthrough pages, each function link is followed by its call site in
`main()` in [`map_builder.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L989) so you can step
through the driver source in parallel.

The toy simulation spans 8 VCF sites. Everything below
is reproducible by running

```
python wiki/generate_wiki.py --page three_generations
```

which regenerates both the figure PNGs referenced here and this
markdown file itself.

## 1. The extended pedigree

![Figure 1 — Three-generation pedigree with outside marriage](fig1.png)

The **G1** layer — `(Dad, Mom) → {Kid1, Kid2, Kid3}` — is exactly
the pedigree the [nuclear-family page](../nuclear_family/nuclear_family.md)
walks through. On this page Kid3 then marries **Spouse**, a fresh
founder, adding the **G2→G3** layer `(Kid3, Spouse) → {GK1, GK2}`.

Kid1 and Kid2 are not drawn in Figure 1, but they are still part
of the pedigree the algorithm processes, and their presence is
load-bearing: the [multi-child guard in
`backfill_sibs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L818) disables partition inference
for single-child families, so Kid3's `A`/`B` paternal labels and `C`
maternal label exist only because Kid1 and Kid2 anchor the partition
at each G1-informative site. Pulling them out of the pedigree would
leave Kid3's letters undetermined, and with them the whole G2→G3
letter trace.

Transmission in the toy simulation:

- **Kid3** carries Dad's ancestral paternal-homolog crossover from
  the nuclear-family page: paternal slot letter `A` on sites 0–3 and
  `B` on sites 4–7; maternal slot letter `C` throughout.
- **Spouse**, a founder, is handed the fresh letter pair `(E, F)` by
  [`Iht::new`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/iht.rs#L172) (driver calls at
  [`map_builder.rs:1059`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L1059) for the master
  template and [`map_builder.rs:1111`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L1111) per VCF
  site), identically to how Dad and Mom were handed `(A, B)` / `(C, D)`
  at the start of the nuclear-family walkthrough.
- **GK1** inherits Spouse's `E` homolog on the paternal slot, plus a
  Kid3 gamete with a crossover: Kid3's paternal homolog (letters
  `A,A,A,A,B,B`) on sites 0–5, then Kid3's maternal homolog (letter
  `C,C`) on sites 6–7.
- **GK2** inherits Spouse's `F` homolog on the paternal slot, plus
  Kid3's paternal homolog unrecombined (letters `A,A,A,A,B,B,B,B`).

## 2. One ancestor-first walk, two triples

The routine
[`track_alleles_through_pedigree`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L295) (driver call
at [`map_builder.rs:1116`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L1116)) is called once
per VCF record and walks every (parent, spouse, child) triple in
ancestor-first depth order given by
[`family.get_individual_depths()`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/ped.rs#L155).

*Depth order* here means: every founder has depth 0, and every
non-founder has depth one greater than the deeper of its two
parents — so a child of a founder and a non-founder sits one
level below the non-founder parent, not below the founder. The
function does a breadth-first sweep from the founders and returns
the list of individuals sorted by depth ascending (ties broken
alphabetically by ID). Iterating triples in
that order guarantees that when a triple `(P, S) → {C ...}` is
visited, both `P` and `S` have already been handled — either as
founders whose slot pairs were pre-filled by [`Iht::new`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/iht.rs#L172),
or as children of a shallower triple whose letters were written by
a previous iteration — so their slot pairs already hold the
letters the current triple needs to read. For this pedigree the
depths are Dad = Mom = Spouse = 0, Kid1 = Kid2 = Kid3 = 1, and
GK1 = GK2 = 2, and the walk visits two triples in this order:

1. **`(Dad, Mom) → {Kid1, Kid2, Kid3}`** — the nuclear-family
   iteration. Its output, for the Kid3 row in particular, is what
   §3.3 / §4.2 of the nuclear-family page already showed: paternal
   slot `A|B` with a block boundary between sites 3 and 4, maternal
   slot `C` throughout.
2. **`(Kid3, Spouse) → {GK1, GK2}`** — the G3 iteration. Same
   function, same code path, but Kid3 now fills the parent role.

The only thing that changes between the two triples is whose slot
pair `find_valid_char` / `get_iht_markers` reads when tagging
carriers. In the first triple it reads Dad's `(A, B)` and Mom's
`(C, D)` — static founder pairs that [`Iht::new`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/iht.rs#L172)
pre-filled at startup. In the second it reads *Kid3's* slot pair,
which [`track_alleles_through_pedigree`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L295) itself
wrote during the first triple. [`get_iht_markers`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L274)
(called from inside the walk at
[`map_builder.rs:328`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L328)) is what reads it back
for the second triple.

One point worth stating explicitly before moving on: **non-founder
parents are first-class** in this walk. Their letters just happen
to vary across sites, whereas founder letters are constant. Kid3
plays the parent role in triple 2 without any special casing beyond
"look up whichever letter sits on the relevant homolog at this
particular site." (§3 below takes up the relationship between the
`(paternal, maternal)` letter pair the walk writes and the
classical inheritance vector of Lander & Green.)

## 3. Relation to the Lander-Green inheritance vector

The per-individual `(paternal-slot letter, maternal-slot letter)`
pairs that §2's walk writes into the `Iht` grid are, up to
reshaping, an *inheritance vector* in the sense of Lander & Green
(*PNAS* 84:2363–2367, [1987](https://doi.org/10.1073/pnas.84.8.2363))
— not a different object. Concatenating the grid's per-individual
rows at a single site gives exactly the length-`2n` vector over
the founder-homolog alphabet that the Lander-Green Hidden Markov
Model manipulates and that tools like
[Merlin](https://doi.org/10.1038/ng786) (Abecasis et al.,
*Nat. Genet.* 30:97–101, 2002) maintain distributions over.

At **site 4** of the toy simulation — the first site after the
ancestral `A → B` crossover on Kid3's paternal homolog — the
`Iht` grid reads

| individual | pair at site 4 |
|---|---|
| Kid1 | `(A, C)` |
| Kid2 | `(B, D)` |
| Kid3 | `(B, C)` |
| GK1  | `(E, B)` |
| GK2  | `(F, B)` |

and stringing the five rows together in a fixed individual-order
produces the length-10 inheritance vector

```
( Kid1_pat, Kid1_mat, Kid2_pat, Kid2_mat, Kid3_pat, Kid3_mat,
  GK1_pat,  GK1_mat,  GK2_pat,  GK2_mat )
= ( A, C, B, D, B, C, E, B, F, B )
```

The Rust pipeline stores it per-individual rather than
concatenated because that lay-out makes the (parent, spouse,
child) triple a natural unit of work, and because the cross-site
post-processing passes
([`collapse_identical_iht`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L385), the two gap-fills,
and [`perform_flips_in_place`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L702)) run
per-individual-column across records. The underlying object is
the same.

Both Lander-Green and `gtg-ped-map` condition on the per-site
genotypes and both produce per-site inheritance vectors. The
difference is *how*.

**Lander-Green** treats the per-site inheritance vector as a
latent state in a hidden Markov model. The emission distribution
at a site is `P(observed genotypes | inheritance vector, founder
alleles)` — an indicator of Mendel-consistency under an error
model, with founder alleles integrated out under an
allele-frequency prior. The transition distribution between
adjacent sites is a recombination model parameterised by a
**genetic map**: each inheritance-vector bit flips across an
inter-site interval with probability tied to the recombination
fraction `θ` on that interval. Forward-backward then returns a
posterior `P(inheritance vector at site s | all observed
genotypes)` at every site, which downstream tools turn into
linkage LOD scores, maximum-likelihood phasing, or imputed
genotypes.

**`gtg-ped-map`**, by contrast, does not build a probabilistic
model. Per site, §2's triple walk deduces the inheritance vector
*deterministically* from Mendelian rules: the unique child
carrying the parent's rare allele is the carrier, its
informative-slot letter is fixed, and
[`backfill_sibs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L804) writes the parent's other
letter on the non-carriers. No emission probabilities; no
integration over founder alleles. Across sites, the
flip / collapse / gap-fill pipeline (listed at the end of §2)
deterministically *minimises the number of letter transitions
between adjacent sites*, which is equivalent to maximising
linkage-block length, or equivalently minimising the total number
of inferred recombinations. That is a Hamming-style parsimony
criterion, not a posterior; no genetic map is consulted.

## 4. Ancestral vs de novo crossovers

![Figure 2 — Ancestral vs de novo crossovers after block collapse](fig2.png)

Two maternal-row transitions appear, and the contrast between them
is the point of this page:

- **Ancestral `A → B` at sites 3/4, in *both* grandchildren's
  maternal rows.** Both GK1 and GK2 inherit Kid3's paternal homolog
  across sites 0–5 (GK1 up to her own crossover at 5/6, GK2
  throughout). The letter on Kid3's paternal homolog itself flips
  `A → B` between sites 3 and 4 — that's the G1-Dad crossover from
  the nuclear-family page — so both grandchildren inherit the flip
  along with the homolog. The transition is **not introduced at
  G2→G3**; it was already present on Kid3's row when the walk
  reached triple 2.
- **De novo `B → C` at sites 5/6, in GK1's maternal row only.**
  Between sites 5 and 6 GK1's gamete from Kid3 crosses over,
  switching from Kid3's paternal homolog to Kid3's maternal one.
  The letter flip is *produced* by that crossover — a "recombinant
  of a recombinant," since it layers on top of the ancestral
  paternal-homolog crossover Kid3 already carried. GK2's gamete
  from Kid3 has no such crossover, so her maternal letter stays on
  `B` through sites 6–7.

The same raw emission rule — one row of `{prefix}.recombinants.txt`
per letter transition per child per slot — governs both, via
[`summarize_child_changes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L673) (driver call at
[`map_builder.rs:1228`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7bb9f882dc0956934f6b4979fd69c04d84b56e98/code/rust/src/bin/map_builder.rs#L1228)). But that rule does
not know which meiosis produced which transition. For this pedigree
the counts come out:

| crossover | meiotic events | rows emitted |
|---|---|---|
| Ancestral `A→B` | 1 (in G1-Dad) | 2 (GK1 m, GK2 m) |
| De novo `B→C`   | 1 (in Kid3's G2→G3 meiosis) | 1 (GK1 m) |

Note that a pairwise comparison of GK1 and GK2's gametes at G3
alone cannot detect the ancestral `A → B` crossover: both
grandchildren inherit Kid3's paternal homolog across sites 0–5,
so their pair relation is unchanged across sites 3/4 even though
the letter on that homolog flips. The letter representation catches
the transition only by chaining back to Kid3's G1 row — an
asymmetry that the ancestor-first walk handles transparently but
that a pairwise pipeline along the lines of
[nuclear_family §5](../nuclear_family/nuclear_family.md#5-an-equivalent-pairwise-comparison-algorithm)
would have to reproduce explicitly.
