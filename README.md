# tapestry

<img src="images/Screenshot 2025-07-23 at 5.04.42 PM.png" alt="XXX" width="500"/>

## Dependencies 

We assume the following tools are in the user's PATH: 
* `bedGraphToBigWig`, `bgzip`, `tabix`
* `gtg-ped-map`, `gtg-concordance` (https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance)
* `hiphase` (https://github.com/PacificBiosciences/HiPhase)
* `aligned_bam_to_cpg_scores` (https://github.com/PacificBiosciences/pb-CpG-tools)

## Workflow 

1. Build an inheritance-based haplotype map and phase variants using the `build-iht-based-haplotype-map-and-phase-variants.sh` script. 
2. Phase variants using read-backed phasing with the `run-hiphase.sh` script.
3. Use `aligned_bam_to_cpg_scores.sh` to generate methylation levels from the haplotagged BAM files produced in step 2.
4. Phase methylation data to founder haplotypes using the `phase-meth-to-founder-haps.sh` script.

## TODO

- [ ] Create the CL tool called "phase-meth-to-founder-haps.sh" from functions defined in `assign-founder-haplotypes-to-methylation-levels.ipynb` 
- [ ] Use excalidraw to create a diagram showing how we get from bams to founder-phased methylation levels
- [ ] Convert manual workflow into a Snakemake workflow
- [ ] Write usage documentation