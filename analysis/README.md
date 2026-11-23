# Slide-deck analyses

Each analysis keeps its methods, reproduction command, results, and slide notes
beside its code. Generated figures and genomic outputs remain under the ignored
`trio_dev_data/output/` directory.

| Analysis | Slide deck | Generated outputs |
| --- | --- | --- |
| [Crossover filtering](crossover_experiments/README.md) | [Discovering crossovers for methylation phasing in trios](https://docs.google.com/presentation/d/1DNbL8Fv1w71YkX4RbH7z1vQOdtUXww-2XkzPDmkIMVo/edit) | `trio_dev_data/output/crossover-experiments/` |
| [Crossover-aware trio methylation](nextflow-trio-model/README.md) | [Phase blocks, crossovers, bit vectors, and child methylation](https://docs.google.com/presentation/d/1fMsmfMTXjtrsqQmHAlmu7VouI0DEg1I7fLOuYRM2h8c/edit) | `trio_dev_data/output/nextflow-trio-model/` |

These are sibling analyses, following the output branches. The original
`crossover_experiments/` source path is retained so existing commands and links
continue to work. `nextflow-trio/` and `nextflow-trio-distance/` are earlier local
output branches; the second deck uses the model run, not those intermediate runs.
