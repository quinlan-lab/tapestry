#!/usr/bin/env python
"""Draw a family pedigree with a painted haplotype pair under each individual.

The family tree comes from a 6-column PED file; the haplotype pair drawn under
each individual comes from a Platinum-Pedigree ``.iht.txt`` haplotype map (the
inheritance blocks from ``gtg-ped-map``) and is painted exactly as in the
consortium's ``code/plotting/plot-iht.R``:

  * each individual contributes two haplotype tracks (hap1 above hap2),
  * a track is a run of horizontal segments, one per haplotype block,
  * each segment is coloured *categorically* by the inherited allele /
    founder-haplotype label using matplotlib's ``turbo`` colormap,
  * founder genotypes are unphased (``A/B``); offspring are phased (``A|B``).

Unlike ``plot-iht.R`` -- which stacks every haplotype in one flat list -- this
script positions each individual's painted pair beneath its node in the pedigree
so that the flow of founder haplotypes through the family is visible.

Before drawing, if the pedigree has apex *trios* (a founder couple with a single
child who is itself a parent), the trio founders (G0) are dropped and the lone
children (G1) -- together with the same haplotype letters where they recur in
descendants -- are relabelled to founder colours A, B, C, D; married-in founders
take the next free letters. This is deliberate: gtg-ped-map cannot resolve
recombination in the lone child of a trio apex, so G1's raw tracks come out flat
and would misleadingly read as "no crossovers". Relabelling keeps G1's colours
connected to its descendants' recombination mosaics. See ``collapse_apex_trios``
for the full, source-line-referenced rationale.

Input contract (Platinum ``gtg-ped-map`` output, space-delimited)::

    #chrom start end marker_count len <ind_1> <ind_2> ... <ind_N> <trailing>

The trailing column is dropped, mirroring ``plot-iht.R``. Per-individual cells
hold a genotype such as ``H1|H2`` (or ``H1/H2`` for unphased founders).

PED columns (Platinum / standard): family, sample, father, mother, sex, pheno.
``sex`` 1=male (square), 2=female (circle). Parents that are referenced but have
no row of their own are added as implicit founders.

Usage::

    python plot_pedigree_haplotypes.py --ped trio.ped --iht trio.iht.txt \
        --out pedigree_haplotypes.pdf [--chrom chr20]

References (Platinum-Pedigree-Inheritance, pinned at e12aca6):
  - Haplotype-map building, .iht.txt format and the plotting method this script
    reproduces:
    https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/HAPLOTYPING.md#plotting-haplotype-maps
  - The original R plotting script (code/plotting/plot-iht.R):
    https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/plotting/plot-iht.R
"""

import argparse
import itertools
import string

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
_MISSING = {"0", "", ".", "na", "nan", "-"}


def parse_ped(path):
    """Return ``{sample: {"sex", "father", "mother"}}`` from a 6-column PED.

    Parents referenced but never given their own row are materialised as
    implicit founders (sex inferred from the parental role they play).
    """
    people = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            _, sample, father, mother, sex = fields[:5]
            father = None if father.lower() in _MISSING else father
            mother = None if mother.lower() in _MISSING else mother
            people[sample] = {"sex": sex, "father": father, "mother": mother}

    # Materialise implicit founders so every referenced parent has a node.
    for attrs in list(people.values()):
        for role, inferred_sex in (("father", "1"), ("mother", "2")):
            parent = attrs[role]
            if parent and parent not in people:
                people[parent] = {"sex": inferred_sex, "father": None, "mother": None}
    return people


def parse_iht(path, chrom=None):
    """Return ``(chrom, {sample: {"hap1": segs, "hap2": segs}}, allele_order)``.

    ``segs`` is a list of ``(start, end, allele)`` tuples. ``allele_order`` is
    the sorted set of allele labels used for a stable colour mapping.
    """
    df = pl.read_csv(
        path,
        separator=" ",
        has_header=True,
        comment_prefix=None,
    )
    # Mirror plot-iht.R: drop the trailing (support) column.
    df = df.select(df.columns[:-1])

    chrom_col = df.columns[0]  # "#chrom"
    fixed = {chrom_col, "start", "end", "marker_count", "len"}
    individuals = [c for c in df.columns if c not in fixed]

    if chrom is None:
        chrom = df[chrom_col][0]
    df = df.filter(pl.col(chrom_col) == chrom).sort("start")

    haplotypes = {ind: {"hap1": [], "hap2": []} for ind in individuals}
    alleles = set()
    for row in df.iter_rows(named=True):
        start, end = int(row["start"]), int(row["end"])
        for ind in individuals:
            cell = str(row[ind])
            parts = cell.replace("|", "/").split("/")
            if len(parts) != 2:
                continue
            a1, a2 = parts
            haplotypes[ind]["hap1"].append((start, end, a1))
            haplotypes[ind]["hap2"].append((start, end, a2))
            alleles.update((a1, a2))

    return chrom, haplotypes, sorted(alleles)


# --------------------------------------------------------------------------- #
# Display transform: collapse apex trios — drop G0, relabel G1 -> A,B,C,D
# --------------------------------------------------------------------------- #
#
# WHY THIS EXISTS — the "first generation never recombines" artifact
# ------------------------------------------------------------------
# In CEPH1463 the apex of the pedigree is two *trios*, not quads:
#
#     NA12889 x NA12890  ->  NA12877   (their ONLY child)
#     NA12891 x NA12892  ->  NA12878   (their ONLY child)
#
# gtg-ped-map (Platinum-Pedigree-Inheritance, the program that produced the
# .iht.txt this script paints) cannot resolve recombination in an individual
# who is the *lone child* of a founder couple. So NA12877 and NA12878 -- the G1
# row, i.e. the parents of the big sibship -- come out as two solid, seam-free
# haplotype tracks (a constant label per haplotype, e.g. NA12877 = E|G), even
# though crossovers certainly happened in the grandparental meioses that made
# them. Those flat G1 tracks are an algorithmic limitation, NOT biology, and are
# misleading in a figure.
#
# Crucially, gtg-ped-map labels each child by the *parent haplotype letter it
# inherited* (find_valid_char hands the child the parent's current label), so
# the sibship below NA12877/NA12878 carries exactly those same letters {E,G} and
# {I,K} as a recombination mosaic -- it is the *sibship* (>= 2 meioses) that
# reveals the crossovers the lone parent could not. The two haplotypes of each
# trio apex therefore already ARE founder haplotypes for everything beneath them.
#
# The display transform is therefore a relabel of G1's existing letters (not a
# fresh painting of G1), so G1's colours stay connected to its descendants'
# recombination mosaics. Given the apex trios exist, we:
#   1. identify the apex trios from the PED (a founder couple whose single shared
#      child is itself a parent);
#   2. relabel the letters carried by the G1 children (NA12877, NA12878) -- and
#      those same letters wherever they recur in descendants -- to A, B, C, D;
#   3. relabel the married-in founders' letters -- and their letters in
#      descendants -- to the next free letters (E, F, ...); then
#   4. draw, mapping the final letters to colours.
# G0 itself is dropped from the plot entirely. Because step 2 keeps G1's real
# (relabelled) segments, a true trio apex still renders as a clean solid A|B / C|D
# pair, but now its colours flow into the sibship's mosaics.
#
# Why gtg-ped-map cannot phase a trio apex (binary = code/rust/src/bin/map_builder.rs,
# all links pinned at commit e12aca6b49ee7208952467db4a2a9e2f79b98efb):
#
#   1. Founders get fixed two-letter labels A/B, C/D, ...; children start blank.
#      Iht::new
#      https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/iht.rs#L172-L198
#
#   2. A transmitted marker is only informative when the parent is heterozygous
#      and carries an allele the spouse lacks. unique_allele
#      https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L243-L268
#
#   3. The label handed to the child is the parent's *first* valid haplotype
#      letter -- find_valid_char returns the first non-'?' char, so a founder
#      (A,B) always contributes 'A'. The code never decides whether the
#      transmitted allele physically sits on haplotype A or B.
#      find_valid_char
#      https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L285-L293
#      flabel pick + child assignment
#      https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L328-L359
#      So a lone child's paternal slot is only ever 'A' or '?' (never 'B'), and
#      its maternal slot only ever 'C' or '?'. A constant label = no crossover.
#
#   4. The two mechanisms that would introduce the complementary founder allele
#      both require >= 2 children and skip a trio:
#        * backfill_sibs -- assigns the non-inherited founder allele (the 'B') to
#          siblings lacking the 'A' marker -- is gated on `children.len() > 1`.
#          https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L804-L818
#        * perform_flips_in_place -- the orientation-consistency pass -- operates
#          only on get_founders_with_multiple_children(), so trio founders are
#          excluded and never flip-corrected.
#          https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L1135-L1140
#          get_flipable_alleles / get_founders_with_multiple_children
#          https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/iht.rs#L554-L581
#          https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/ped.rs#L264-L281
#
#   5. A recombination is emitted only when a child's label *changes* between
#      consecutive blocks; a constant-or-missing track yields zero recombinants.
#      summarize_child_changes
#      https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/code/rust/src/bin/map_builder.rs#L673-L700
#
# This is why the HAPLOTYPING.md note says gtg-ped-map phases "quads and larger":
# resolving a founder's two transmitted haplotypes needs >= 2 observed meioses
# (>= 2 children). The 11 G3 siblings below satisfy this and DO show crossovers;
# NA12877/NA12878 cannot and do not.
def _label_stream():
    """Yield 'A','B',...,'Z','AA','AB',... -- an unbounded supply of labels."""
    for width in itertools.count(1):
        for combo in itertools.product(string.ascii_uppercase, repeat=width):
            yield "".join(combo)


def _children_of(people):
    """Map each individual to the list of individuals that name it as a parent."""
    kids = {pid: [] for pid in people}
    for pid, attrs in people.items():
        for parent in (attrs["father"], attrs["mother"]):
            if parent in kids:
                kids[parent].append(pid)
    return kids


def _distinct_labels(hap):
    """Allele labels in a {hap1,hap2} record, hap1 before hap2, de-duplicated."""
    seen = []
    for which in ("hap1", "hap2"):
        for _, _, allele in hap.get(which, []):
            if allele not in seen:
                seen.append(allele)
    return seen


def find_apex_trios(people):
    """Return ``(g0_apex, g1)`` for trios sitting at the top of the pedigree.

    An apex trio is a couple **both of whom are founders** (no parents) whose
    only shared child is a single individual that is itself a parent. ``g0_apex``
    is the set of those founder parents (to be dropped); ``g1`` is the list of
    their lone children (the row to be relabelled A/B, C/D, ...). For CEPH1463
    this finds the two trios NA12889xNA12890->NA12877 and NA12891xNA12892->NA12878.
    """
    founders = {pid for pid, a in people.items()
                if not a["father"] and not a["mother"]}
    kids = _children_of(people)

    g0_apex, g1, seen_couples = set(), [], set()
    for pid, attrs in people.items():
        father, mother = attrs["father"], attrs["mother"]
        if father not in founders or mother not in founders:
            continue
        couple = frozenset((father, mother))
        if couple in seen_couples:
            continue
        seen_couples.add(couple)
        shared = [c for c in kids[father] if c in kids[mother]]
        # Lone child who is itself a parent => an apex trio above the sibship.
        if len(shared) == 1 and kids[shared[0]]:
            g0_apex.update((father, mother))
            g1.append(shared[0])
    return g0_apex, sorted(g1)


def collapse_apex_trios(people, haplotypes):
    """Drop apex-trio founders (G0) and relabel G1 + descendants to A,B,C,D,...

    Returns ``(people, haplotypes, allele_order)``. If the PED has no apex trios
    (see :func:`find_apex_trios`) the inputs are returned unchanged (bar a fresh
    ``allele_order``), so non-trio pedigrees are drawn exactly as before.

    When apex trios exist, the relabel follows the four-step recipe documented
    above: G1's haplotype letters (and the same letters in descendants) become
    A,B,C,D,...; married-in founders' letters (and theirs in descendants) take
    the next free letters; G0 is removed entirely. The relabel is a single global
    bijection applied to every track, so a colour means the same haplotype
    everywhere it appears -- G1's founder colours flow down into the sibship.
    """
    def order(): return sorted({
        al for h in haplotypes.values()
        for which in ("hap1", "hap2") for _, _, al in h[which]
    })

    g0_apex, g1 = find_apex_trios(people)
    if not g1:
        return people, haplotypes, order()

    # 1. Drop G0; sever G1's now-dangling parent links so G1 render as founders.
    for pid in g0_apex:
        people.pop(pid, None)
        haplotypes.pop(pid, None)
    for attrs in people.values():
        for role in ("father", "mother"):
            if attrs[role] in g0_apex:
                attrs[role] = None

    # Build one global label bijection, consuming A,B,C,D,... in priority order.
    labels = _label_stream()
    remap = {}

    def reserve(pid):
        for allele in _distinct_labels(haplotypes.get(pid, {})):
            if allele not in remap:
                remap[allele] = next(labels)

    # 2. G1 letters first  -> A, B, C, D, ...
    for pid in g1:
        reserve(pid)
    # 3. then the married-in founders' letters -> next free letters.
    marryins = sorted(pid for pid, a in people.items()
                      if not a["father"] and not a["mother"] and pid not in g1)
    for pid in marryins:
        reserve(pid)
    # Defensive: any label only ever seen in descendants gets the next free slot
    # (in real gtg-ped-map output every descendant label traces to a G1 or
    # married-in founder, so this is a no-op there).
    for pid in sorted(haplotypes):
        reserve(pid)

    # 4. Apply the bijection to every track, then recolour from the new labels.
    for hap in haplotypes.values():
        for which in ("hap1", "hap2"):
            hap[which] = [(s, e, remap[al]) for s, e, al in hap[which]]

    return people, haplotypes, order()


# --------------------------------------------------------------------------- #
# Pedigree layout (generational rows + barycenter x-relaxation)
# --------------------------------------------------------------------------- #
def _couples(people):
    """Frozensets of co-parent pairs that share at least one child."""
    pairs = set()
    for attrs in people.values():
        f, m = attrs["father"], attrs["mother"]
        if f and m:
            pairs.add(frozenset((f, m)))
    return pairs


def assign_levels(people):
    """Generation index per person (founders = 0), with spouses aligned."""
    level = {pid: 0 for pid in people}
    couples = _couples(people)
    for _ in range(len(people) + 1):
        changed = False
        for pid, attrs in people.items():
            parents = [p for p in (attrs["father"], attrs["mother"]) if p]
            if parents:
                want = max(level[p] for p in parents) + 1
                if want != level[pid]:
                    level[pid], changed = want, True
        for couple in couples:  # keep mates on the same row
            a, b = tuple(couple)
            top = max(level[a], level[b])
            if level[a] != top or level[b] != top:
                level[a] = level[b] = top
                changed = True
        if not changed:
            break
    return level


def _initial_order(people, level):
    """A DFS preorder from founders, giving each row a sane starting order."""
    children = {pid: [] for pid in people}
    for pid, attrs in people.items():
        for parent in (attrs["father"], attrs["mother"]):
            if parent:
                children[parent].append(pid)

    spouse = {}
    for couple in _couples(people):
        a, b = tuple(couple)
        spouse.setdefault(a, b)
        spouse.setdefault(b, a)

    order = {lvl: [] for lvl in set(level.values())}
    seen = set()

    def visit(pid):
        if pid in seen:
            return
        seen.add(pid)
        order[level[pid]].append(pid)
        mate = spouse.get(pid)
        if mate and mate not in seen and level[mate] == level[pid]:
            seen.add(mate)
            order[level[mate]].append(mate)
        kids = children[pid] + (children.get(mate, []) if mate else [])
        for kid in sorted(set(kids)):
            visit(kid)

    for founder in sorted(p for p in people if level[p] == 0):
        visit(founder)
    for pid in people:  # anyone unreached (defensive)
        visit(pid)
    return order


def layout(people, slot=2.6, gap=1.1, row_height=4.0, iterations=40):
    """Return ``{pid: (x, y)}`` placing founders at the top.

    A light Sugiyama-style relaxation: alternately pull each node toward the
    mean x of its neighbours (parents, children, spouse), then re-pack each row
    left-to-right keeping spouse pairs glued together. Robust to the
    cross-marriages that defeat a plain tree layout.
    """
    level = assign_levels(people)
    order = _initial_order(people, level)

    children = {pid: [] for pid in people}
    for pid, attrs in people.items():
        for parent in (attrs["father"], attrs["mother"]):
            if parent:
                children[parent].append(pid)
    spouse = {}
    for couple in _couples(people):
        a, b = tuple(couple)
        spouse.setdefault(a, b)
        spouse.setdefault(b, a)

    x = {}
    for lvl, ids in order.items():
        for i, pid in enumerate(ids):
            x[pid] = i * slot

    def repack(ids):
        """Place ids left-to-right honouring min spacing, gluing spouse pairs."""
        blocks, used = [], set()
        for pid in ids:
            if pid in used:
                continue
            mate = spouse.get(pid)
            if mate in ids and mate not in used:
                used.update((pid, mate))
                blocks.append([pid, mate])
            else:
                used.add(pid)
                blocks.append([pid])
        blocks.sort(key=lambda b: np.mean([x[p] for p in b]))
        members = [pid for block in blocks for pid in block]
        desired = float(np.mean([x[pid] for pid in members]))  # barycenter to keep
        cursor = 0.0
        for block in blocks:
            for j, pid in enumerate(block):
                x[pid] = cursor + j * slot
            cursor += len(block) * slot + gap
        # Translate the whole (now overlap-free) row so its centre of mass is
        # unchanged: this preserves the cross-row barycenter alignment instead of
        # slamming every row against x=0, which would un-centre children.
        shift = desired - float(np.mean([x[pid] for pid in members]))
        for pid in members:
            x[pid] += shift

    levels_sorted = sorted(order)
    for it in range(iterations):
        sweep = levels_sorted if it % 2 == 0 else levels_sorted[::-1]
        for lvl in sweep:
            ids = order[lvl]
            for pid in ids:
                neigh = []
                neigh += [x[c] for c in children[pid] if c in x]
                for parent in (people[pid]["father"], people[pid]["mother"]):
                    if parent:
                        neigh.append(x[parent])
                if spouse.get(pid) in x:
                    neigh.append(x[spouse[pid]])
                if neigh:
                    x[pid] = float(np.mean(neigh))
            order[lvl] = sorted(ids, key=lambda p: x[p])
            repack(order[lvl])

    # Final top-down pass: anchor each individual under its parents (founders
    # over their children) so children end up exactly centred beneath the
    # marriage midpoint the parent->child connector drops from.
    for lvl in levels_sorted:
        for pid in order[lvl]:
            parents = [p for p in (people[pid]["father"], people[pid]["mother"])
                       if p]
            anchor = ([x[p] for p in parents] if parents
                      else [x[c] for c in children[pid] if c in x])
            if anchor:
                x[pid] = float(np.mean(anchor))
        order[lvl] = sorted(order[lvl], key=lambda p: x[p])
        repack(order[lvl])

    return {pid: (x[pid], -level[pid] * row_height) for pid in people}, level


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #
def _merge_runs(segs):
    """Collapse consecutive segments sharing an allele into single runs.

    The iht map stores one block per marker interval, so a stretch of the same
    inherited haplotype arrives as many adjacent same-allele segments. Merging
    them into one rectangle means a solid founder track becomes a single vector
    shape -- both seam-free and far fewer objects to wrangle in Illustrator.
    """
    runs = []
    for start, end, allele in segs:
        if runs and runs[-1][2] == allele and runs[-1][1] == start:
            runs[-1] = (runs[-1][0], end, allele)
        else:
            runs.append((start, end, allele))
    return runs


def _paint_pair(ax, cx, cy, haps, xmin, xspan, colors, paint_w, track_h, gap_h):
    """Draw the hap1/hap2 painted tracks centred at (cx, cy).

    Each iht block is a *vector* rectangle so the painting stays fully editable
    in Illustrator (each haplotype block is its own selectable, recolourable
    shape). To avoid the antialiasing hairlines that abutting vector fills leave
    along every shared edge -- visible in PDF viewers even inside a solid-colour
    founder track -- two things are done:

      * consecutive same-allele blocks are merged (``_merge_runs``), so a solid
        stretch is one rectangle with no internal edges at all, and
      * every rectangle is stroked in *its own fill colour*, so the stroke
        overpaints the seam at boundaries between differently-coloured blocks.
    """
    left = cx - paint_w / 2.0
    for k, hap in enumerate(("hap1", "hap2")):
        top = cy - k * (track_h + gap_h)
        ax.add_patch(Rectangle((left, top - track_h), paint_w, track_h,
                               facecolor="none", edgecolor="0.6", linewidth=0.4,
                               zorder=2))
        for start, end, allele in _merge_runs(haps[hap]):
            sx = left + (start - xmin) / xspan * paint_w
            w = max((end - start) / xspan * paint_w, paint_w * 0.004)
            color = colors[allele]
            # Same-colour stroke: fills the hairline at shared edges while
            # keeping the block a clean, individually editable vector shape.
            ax.add_patch(Rectangle((sx, top - track_h), w, track_h,
                                   facecolor=color, edgecolor=color,
                                   linewidth=0.3, zorder=3))


def _build_colors(allele_order, palette="auto"):
    """Map each allele to a colour.

    ``turbo`` (continuous, as in ``plot-iht.R``) is faithful but blends adjacent
    categories, so for many founder haplotypes ``auto`` switches to categorical
    ``tab20`` (its dark/light hue pairs keep a founder's two haplotypes apart);
    beyond 20 it extends with ``tab20b``/``tab20c`` (60 distinct colours).
    """
    n = len(allele_order)
    if palette == "turbo" or (palette == "auto" and n <= 10):
        cols = list(plt.cm.turbo(np.linspace(0.05, 0.95, n)))
    else:  # "tab20" or auto with many alleles -> categorical
        base = (list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors)
                + list(plt.cm.tab20c.colors))
        cols = [base[i % len(base)] for i in range(n)]
    return {a: cols[i] for i, a in enumerate(allele_order)}


def draw(people, pos, level, haplotypes, allele_order, chrom, out,
         title=None, palette="auto", mode_note=None):
    paint_w, track_h, gap_h = 1.9, 0.55, 0.12
    # Symbols are enlarged from a bare node marker to a label-bearing shape: the
    # individual's id is drawn *inside* the circle/square, so the radius is sized
    # to comfortably hold a sample id such as "NA12877" rather than a dot.
    sym_r = 1.05

    starts = [s for h in haplotypes.values() for s, _, _ in h["hap1"]]
    ends = [e for h in haplotypes.values() for _, e, _ in h["hap1"]]
    xmin, xmax = min(starts), max(ends)
    xspan = max(xmax - xmin, 1)

    colors = _build_colors(allele_order, palette)

    # Size the canvas from the actual layout extent (scales in both the widest
    # row and the number of generations) rather than a raw headcount, so deep
    # and/or wide pedigrees are neither cramped nor padded with whitespace.
    node_xs = [p[0] for p in pos.values()]
    node_ys = [p[1] for p in pos.values()]
    pair_height = sym_r + 0.35 + 2 * track_h + gap_h  # symbol bottom -> paint bottom
    width_units = (max(node_xs) - min(node_xs)) + paint_w + 2.0
    height_units = (max(node_ys) - min(node_ys)) + sym_r + pair_height + 2.0
    scale = 0.6  # inches per layout unit
    fig, ax = plt.subplots(figsize=(max(7.0, width_units * scale),
                                    max(5.0, height_units * scale)))

    # ---- connectors (drawn first, under the symbols) ----
    line_kw = dict(color="0.35", linewidth=1.0, zorder=1)
    for couple in _couples(people):
        a, b = tuple(couple)
        (xa, ya), (xb, yb) = pos[a], pos[b]
        ax.plot([xa, xb], [ya, yb], **line_kw)
        mx, my = (xa + xb) / 2.0, (ya + yb) / 2.0
        kids = [k for k, at in people.items()
                if {at["father"], at["mother"]} >= {a, b}]
        if not kids:
            continue
        # Drop the marriage line to a sibship bar that sits below the parents'
        # painted pair, then branch up to each child. The bar is widened to span
        # the marriage midpoint as well as every child, so the drop always meets
        # the children's risers even if a child is not perfectly centred.
        paint_bottom = my - sym_r - 0.35 - 2 * track_h - gap_h
        sib_y = paint_bottom - 0.4
        ax.plot([mx, mx], [my, sib_y], **line_kw)
        kid_xs = [pos[k][0] for k in kids]
        ax.plot([min(kid_xs + [mx]), max(kid_xs + [mx])], [sib_y, sib_y],
                **line_kw)
        for k in kids:
            kx, ky = pos[k]
            ax.plot([kx, kx], [sib_y, ky + sym_r], **line_kw)

    # Single-parent links (only one parent present in the family).
    for pid, attrs in people.items():
        parents = [p for p in (attrs["father"], attrs["mother"]) if p]
        if len(parents) == 1:
            px, py = pos[parents[0]]
            kx, ky = pos[pid]
            ax.plot([px, px], [py, ky + sym_r], **line_kw)

    # ---- individuals: sex symbol + painted haplotype pair ----
    for pid, (cx, cy) in pos.items():
        if people[pid]["sex"] == "1":
            ax.add_patch(Rectangle((cx - sym_r, cy - sym_r), 2 * sym_r, 2 * sym_r,
                                   facecolor="white", edgecolor="black",
                                   linewidth=1.2, zorder=4))
        else:
            ax.add_patch(Circle((cx, cy), sym_r, facecolor="white",
                                edgecolor="black", linewidth=1.2, zorder=4))
        ax.text(cx, cy, pid, ha="center", va="center",
                fontsize=9, zorder=5)
        _paint_pair(ax, cx, cy - sym_r - 0.35, haplotypes[pid], xmin, xspan,
                    colors, paint_w, track_h, gap_h)

    # ---- legend / cosmetics ----
    handles = [Line2D([0], [0], marker="s", linestyle="none", markersize=22,
                      markerfacecolor=colors[a], markeredgecolor="0.4", label=a)
               for a in allele_order]
    legend = ax.legend(handles=handles, title="Allele / founder haplotype",
                       loc="center left", bbox_to_anchor=(1.01, 0.5),
                       fontsize=18, labelspacing=1.0, handletextpad=0.8,
                       borderpad=1.2, frameon=False)
    legend.get_title().set_fontsize(21)

    span_mb = f"{xmin / 1e6:.2f}–{xmax / 1e6:.2f} Mb"
    note = f", {mode_note}" if mode_note else ""
    ax.set_title(title or f"Pedigree haplotype painting — {chrom} ({span_mb}{note})")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.margins(0.08)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ped", required=True, help="6-column PED file")
    ap.add_argument("--iht", required=True, help="Platinum .iht.txt haplotype map")
    ap.add_argument("--out", required=True, help="output image (.png/.pdf/.svg)")
    ap.add_argument("--chrom", default=None,
                    help="chromosome to plot (default: first in the iht file)")
    ap.add_argument("--palette", default="auto",
                    choices=["auto", "turbo", "tab20"],
                    help="allele colours: 'turbo' (as plot-iht.R), 'tab20' "
                         "(categorical), or 'auto' (turbo<=10 alleles, else tab20)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--no-collapse", dest="collapse", action="store_false",
                    help="draw the RAW gtg-ped-map output verbatim: keep G0 and "
                         "leave every label as-is (no apex-trio drop/relabel). Use "
                         "this to show the 'before' beside the default 'after'.")
    args = ap.parse_args()

    people = parse_ped(args.ped)
    chrom, haplotypes, allele_order = parse_iht(args.iht, args.chrom)

    # Restrict to individuals present in both the PED and the haplotype map.
    people = {p: a for p, a in people.items() if p in haplotypes}
    missing = [p for p in people if not haplotypes[p]["hap1"]]
    for p in missing:  # keep the node, but it simply has no painted segments
        haplotypes.setdefault(p, {"hap1": [], "hap2": []})

    if args.collapse:
        # If the PED has apex trios, drop their founders (G0) and relabel the lone
        # children (G1) -- and their letters in descendants -- to A,B,C,D, with the
        # married-in founders taking the next free letters. See collapse_apex_trios()
        # for the full rationale: gtg-ped-map cannot resolve recombination in the lone
        # child of a trio apex, so G1's flat tracks would otherwise read misleadingly
        # as "no crossovers".
        people, haplotypes, allele_order = collapse_apex_trios(people, haplotypes)

    # Make the two renders self-labelling so a before/after pair is unambiguous.
    mode_note = ("raw gtg-ped-map output" if not args.collapse
                 else "apex trios relabelled A,B,C,D")
    pos, level = layout(people)
    out = draw(people, pos, level, haplotypes, allele_order, chrom, args.out,
               title=args.title, palette=args.palette, mode_note=mode_note)
    print(f"wrote {out}  ({len(people)} individuals, {chrom}, "
          f"{len(allele_order)} alleles, "
          f"{'raw' if not args.collapse else 'collapsed'})")


if __name__ == "__main__":
    main()
