# tapestry

<img src="images/tapestry.png" alt="XXX" width="900"/>

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

1. Phase variants:
   - Build an inheritance-based haplotype map and phase variants using the `build-iht-based-haplotype-map-and-phase-variants.sh` script.
   - Phase variants using read-backed phasing with the `run-hiphase.sh` script.
2. Use `aligned_bam_to_cpg_scores.sh` to generate methylation levels from the haplotagged BAM files produced in step 1.
3. Phase methylation data to founder haplotypes using the `phase_meth_to_founder_haps.sh` script, which uses the data produced in steps 1 and 2.

## TODO

- [ ] Create a python module (called "tapestry") for analyzing founder-phased DNA methylation data across generations in a pedigree, e.g., to answer the question "given a haplotype in a child, has its methylation changed relative to the same haplotype in the parent" 
- [ ] Use excalidraw to create a diagram showing how we get from bams to founder-phased methylation levels
- [ ] Convert manual workflow into a Snakemake workflow 
