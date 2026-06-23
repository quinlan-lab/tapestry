# Apex-trio relabeling: before vs. after

This worked example shows what `src/plot_pedigree_haplotypes.py` does to the raw
output of [gtg-ped-map](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/HAPLOTYPING.md),
and why. The pedigree mirrors the apex of CEPH1463: two grandparent couples
(`GP1×GP2`, `GP3×GP4`) each have a **single** child (`DAD`, `MOM`) who together
parent a sibship (`KID1..KID3`); `MARRYIN` marries in and parents `GKID` with
`KID1`.

## The two figures

| file | command | what it shows |
|------|---------|---------------|
| `before_raw_gtg-ped-map.pdf` | `--no-collapse` | gtg-ped-map's output **verbatim** |
| `after_apex-trio_relabel.pdf` | (default) | after dropping G0 and relabeling G1 |

**Before** — `DAD` and `MOM` are each the lone child of a grandparent couple (an
apex *trio*). gtg-ped-map cannot resolve recombination in the lone child of a
trio, so `DAD` (`A|C`) and `MOM` (`E|G`) are painted as two flat, seam-free
tracks with no crossovers — even though crossovers certainly happened in the
grandparental meioses that produced them. Read literally, the first generation
looks recombination-free, which is misleading. (See the long comment in
`collapse_apex_trios` for the gtg-ped-map source lines that cause this.)

**After** — the transform:
1. identifies the apex trios from the PED;
2. relabels the G1 haplotype letters — and the same letters where they recur in
   descendants — to founder colors `A, B, C, D` (`DAD → A/B`, `MOM → C/D`);
3. relabels the married-in founder's letters (and theirs in descendants) to the
   next free letters (`E, F, …`);
4. drops G0 and draws.

Because it is a single global *relabel*, a color means the same
haplotype everywhere it appears: `DAD`'s `A`/`B` reappear in the sibship as a
recombination mosaic (e.g. `KID1` switches `A→B` along `chr1`), and `GKID`'s
lower track ties back through `KID1` to `DAD`. The sibship (≥2 meioses) is what
makes those crossovers visible — exactly the information a trio apex lacks.

## Regenerate

```bash
D=examples/apex_trio_relabeling
python src/plot_pedigree_haplotypes.py --ped $D/example.ped --iht $D/example.iht.txt \
    --no-collapse --out $D/before_raw_gtg-ped-map.pdf
python src/plot_pedigree_haplotypes.py --ped $D/example.ped --iht $D/example.iht.txt \
    --out $D/after_apex-trio_relabel.pdf
```

`example.iht.txt` is Mendelian-consistent by construction (every child allele is
carried by the named parent at every site); the founder labels `A/B … I/J` are
assigned exactly as gtg-ped-map would (founders sorted by ID).
