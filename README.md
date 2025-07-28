# tapestry

<img src="images/Screenshot 2025-07-23 at 5.04.42 PM.png" alt="XXX" width="500"/>

## Dependencies 

We assume the following command-line tools are in the user's PATH: 

* `bedGraphToBigWig`, `bgzip`, `tabix`
* `gtg-ped-map`, `gtg-concordance` (https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance)
* `hiphase` (https://github.com/PacificBiosciences/HiPhase)
* `aligned_bam_to_cpg_scores` (https://github.com/PacificBiosciences/pb-CpG-tools)

Install the python dependencies:

```
/path/to/python3.11 -m venv .venv 
source .venv/bin/activate 

# https://github.com/brentp/cyvcf2?tab=readme-ov-file#github-building-htslib-and-cyvcf2-from-source
cd $HOME 
git clone --recursive https://github.com/brentp/cyvcf2
cd cyvcf2
pip install -r requirements.txt 
CYVCF2_HTSLIB_MODE=EXTERNAL python setup.py install

cd /path/to/tapestry
pip install -r requirements.txt
```

## Workflow 

1. Build an inheritance-based haplotype map and phase variants using the `build-iht-based-haplotype-map-and-phase-variants.sh` script. 
2. Phase variants using read-backed phasing with the `run-hiphase.sh` script.
3. Use `aligned_bam_to_cpg_scores.sh` to generate methylation levels from the haplotagged BAM files produced in step 2.
4. Phase methylation data to founder haplotypes using the `phase-meth-to-founder-haps.sh` script.

## TODO

- [ ] Create the CL tool called "phase-meth-to-founder-haps.sh" from functions defined in `assign-founder-haplotypes-to-methylation-levels.ipynb` 
- [ ] Create a python module (called "tapestry") for analyzing founder-phased DNA methylation data across generations in a pedigree, e.g., to answer the question "given a haplotype in a child, has its methylation changed relative to the same haplotype in the parent" 
- [ ] Use excalidraw to create a diagram showing how we get from bams to founder-phased methylation levels
- [ ] Convert manual workflow into a Snakemake workflow
- [ ] Write usage documentation