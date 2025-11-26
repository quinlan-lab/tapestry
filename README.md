# tapestry

A pipeline to phase DNA methylation from HiFi reads in a human pedigree to the haplotypes of the pedigree's founders. 

## Dependencies

We assume the following command-line tools are in the user's PATH: 

* `bedGraphToBigWig`, `bgzip`, `tabix`, `bcftools`
* `gtg-ped-map`, `gtg-concordance` (https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance)
* `hiphase` (https://github.com/PacificBiosciences/HiPhase)
* `aligned_bam_to_cpg_scores` (https://github.com/PacificBiosciences/pb-CpG-tools)

## Installation 

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

## Pedigree-wise workflow 

1. Phase variants:
   - Phase variants using read-backed phasing with the `run-hiphase.sh` script.
   - Build an inheritance-based haplotype map and inheritance-phase variants using the `build-iht-based-haplotype-map-and-phase-variants.sh` script.
2. Use `aligned_bam_to_cpg_scores.sh` to generate count- and model-based methylation levels from the haplotagged BAM files produced in step 1.
3. Phase count- and model-based methylation data to founder haplotypes using the `phase_meth_to_founder_haps.sh` script, which uses the data produced in steps 1 and 2.
4. Use `expand_to_all_cpgs.sh` script to generalize tapestry's output to include all CpG sites in the reference and sample genome, and unphased count- and model-based methylation levels, where available. This step also uses heterozygous SNVs to flag CpG sites that are "allele-specific", i.e., exist in only one of the haplotypes in the sample in question. Uses output of steps 2 and 3.

## Step 3 of workflow: Phasing count- and model-based DNA methylation to founder haplotypes 

The `phase_meth_to_founder_haps.sh` script calls a CL tool called `src/phase_meth_to_founder_haps.py`: 

```
vcf_read_phased="${read_phased_dir}/${uid}.GRCh38.deepvariant.glnexus.phased.vcf.gz" # single-sample vcf from hiphase
tsv_read_phase_blocks="${read_phased_dir}/${uid}.GRCh38.hiphase.blocks.tsv" # single-sample tsv from hiphase
vcf_iht_phased="${iht_phased_dir}/CEPH1463.GRCh38.pass.sorted.vcf.gz" # joint-called multi-sample vcf from gtg-ped-map/gtg-concordance
txt_iht_blocks="${iht_phased_dir}/CEPH1463.GRCh38.iht.sorted.txt" # multi-sample iht blocks file from gtg-ped-map/gtg-concordance
bed_meth_count_hap1="${meth_count_read_phased_dir}/${uid}.GRCh38.haplotagged.hap1.bed.gz" # bed file of count-based methylation from aligned_bam_to_cpg_scores for hap1
bed_meth_count_hap2="${meth_count_read_phased_dir}/${uid}.GRCh38.haplotagged.hap2.bed.gz" # bed file of count-based methylation from aligned_bam_to_cpg_scores for hap2
bed_meth_model_hap1="${meth_model_read_phased_dir}/${uid}.GRCh38.haplotagged.hap1.bed.gz" # bed file of model-based methylation from aligned_bam_to_cpg_scores for hap1
bed_meth_model_hap2="${meth_model_read_phased_dir}/${uid}.GRCh38.haplotagged.hap2.bed.gz" # bed file of model-based methylation from aligned_bam_to_cpg_scores for hap2

nohup python src/phase_meth_to_founder_haps.py \
    --uid ${uid} \
    --vcf_read_phased ${vcf_read_phased} \
    --tsv_read_phase_blocks ${tsv_read_phase_blocks} \
    --vcf_iht_phased ${vcf_iht_phased} \
    --txt_iht_blocks ${txt_iht_blocks} \
    --bed_meth_count_hap1 ${bed_meth_count_hap1} \
    --bed_meth_count_hap2 ${bed_meth_count_hap2} \
    --bed_meth_model_hap1 ${bed_meth_model_hap1} \
    --bed_meth_model_hap2 ${bed_meth_model_hap2} \
    --output_dir ${output_dir} \
    > ${output_dir}/${uid}.log 2>&1 &

```

This will produce a log file that looks like: 

```
2025-10-06 14:41:36 - INFO - Got read-based phasing data: 2445352 rows, 8 columns
2025-10-06 14:41:37 - INFO - Got read-based phase blocks: 10017 rows, 7 columns
2025-10-06 14:44:45 - INFO - Got inheritance-based phasing data: 2192365 rows, 7 columns
2025-10-06 14:44:45 - INFO - Got inheritance-based phase blocks: 1329 rows, 5 columns
2025-10-06 14:44:52 - INFO - Got all phasing data: 2189885 rows, 15 columns
2025-10-06 14:44:59 - INFO - Got hap map: 8904 rows, 7 columns
2025-10-06 14:44:59 - INFO - Got sites: 2189885 rows, 5 columns
2025-10-06 14:44:59 - INFO - Got sites where read-based and inheritance-based bit vectors don't match: 33390 rows, 5 columns
2025-10-06 14:45:00 - INFO - Wrote paternal and maternal hap-map blocks for IGV visualization
2025-10-06 14:45:00 - INFO - Wrote hap-map blocks
2025-10-06 14:45:00 - INFO - Wrote sites of bit-vectors, and sites where bit vectors are mismatched, for IGV visualization
2025-10-06 14:45:31 - INFO - Got read-based phasing of count-based methylation levels: 26729958 rows, 7 columns
2025-10-06 14:45:58 - INFO - Got read-based phasing of model-based methylation levels: 26729958 rows, 7 columns
2025-10-06 14:47:43 - INFO - Phased count-based methylation levels to founder haplotypes: 26729958 rows, 14 columns
2025-10-06 14:49:26 - INFO - Phased model-based methylation levels to founder haplotypes: 26729958 rows, 14 columns
2025-10-06 14:49:52 - INFO - Combined count- and model-based methylation levels: 26729958 rows, 16 columns
2025-10-06 14:49:52 - INFO - Percentage of CpG sites that are phased to founder haplotypes: 93%
2025-10-06 14:49:52 - INFO - Percentage of CpG sites that are within 50bp of a mismatch site: 0.183%
2025-10-06 14:49:57 - INFO - Wrote methylation levels phased to founder haplotypes
2025-10-06 14:52:15 - INFO - Wrote bigwig file for pat count-based methylation levels
2025-10-06 14:54:37 - INFO - Wrote bigwig file for pat model-based methylation levels
2025-10-06 14:56:54 - INFO - Wrote bigwig file for mat count-based methylation levels
2025-10-06 14:59:16 - INFO - Wrote bigwig file for mat model-based methylation levels
2025-10-06 14:59:17 - INFO - Done running /scratch/ucgd/lustre-labs/quinlan/u6018199/tapestry/src/phase_meth_to_founder_haps.py
```

and phases count- and model-based DNA methylation to founder haplotypes. 
The tool also creates files that collectively enable visualization of phased DNA methylation in IGV, e.g., 

<img src="images/tapestry.png" alt="XXX" width="900"/>

## Step 4 of workflow: Expanding the dataset to include all CpG sites in the reference and the sample genome, and to include unphased DNA methylation, and to flag CpGs as "allele-specific" 

We would like to be aware of the existence of all CpG sites in the reference genome, irrespective of whether phased or unphased methylation levels are reported for those sites. Additionally, DNA methylation levels are computed at all CpG sites in the sample's haplotypes, even if those CpG sites are not in the reference genome (an example of a variant "creating" a CpG site in the sample). Conversely, CpG sites can be destroyed by mutation, and therefore a methylation level is not reported at that site on the corresponding haplotype. We call such CpG sites "allele-specific". Finally, this step of the workflow reports various QC statistics, e.g., the percentage of CpG sites (in reference and sample genomes, and on phasable chromosomes) at which count-based methylation is phased to at least one parental haplotype. See `https://github.com/quinlan-lab/tapestry/blob/main/src/expand_to_all_cpgs.200081.ipynb` for examples of all these things. 

### Bed file format

The bed file output by step 4 includes a header containing run metadata and column headings. The header format provides key/value pairs as "##{key}={value}", and a final column header line in the form: "#col1 col2 col3...". An example header is shown below:

```
##source='/scratch/ucgd/lustre-labs/quinlan/u6018199/tapestry/src/expand_to_all_cpgs.py with args {'bed_all_cpgs_in_reference': '/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased.all-cpgs.2/all_cpg_sites_in_reference.bed', 'bed_meth_count_unphased': '/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.count.read-backed-phased/200084.GRCh38.haplotagged.combined.bed.gz', 'bed_meth_model_unphased': '/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.model.read-backed-phased/200084.GRCh38.haplotagged.combined.bed.gz', 'bed_meth_founder_phased': '/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased/200084.dna-methylation.founder-phased.bed', 'bed_het_site_mismatches': '/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased/200084.bit-vector-sites-mismatches.bed', 'bed_meth_founder_phased_all_cpgs': '/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased.all-cpgs.2/200084.dna-methylation.founder-phased.all_cpgs.bed', 'uid': '200084', 'vcf_iht_phased': '/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38/CEPH1463.GRCh38.pass.sorted.vcf.gz'}'
#chrom  start_cpg       end_cpg total_read_count        methylation_level_count methylation_level_model start_hap_map_block
     end_hap_map_block       haplotype_concordance_in_hap_map_block  num_het_SNVs_in_hap_map_block   total_read_count_pat    total_read_count_mat    founder_haplotype_pat   founder_haplotype_mat   methylation_level_pat_count     methylation_level_mat_count     methylation_level_pat_model     methylation_level_mat_model     cpg_is_within_50bp_of_mismatch_site     cpg_overlaps_at_least_one_snv   snv_genotypes   cpg_is_allele_specific
```

Definitions of column headings: 

Column | Definition 
:--- | :---
chrom | chromosome on which CpG site is found 
start_cpg | start coordinate of CpG dinucleotide 
end_cpg | end coordinate of CpG dinucleotide 
total_read_count | total number of HiFi reads overlapping CpG site
methylation_level_count | methylation level: fraction of reads that are methylated (inferred from sequencing kinetics) at given CpG site, c.f., https://github.com/PacificBiosciences/pb-CpG-tools?tab=readme-ov-file#output-modes-and-option-details
methylation_level_model | methylation level: neural network maps sequencing kinetics at multiple CpG sites to methylation probability at a central CpG site, c.f., "Methods" section of https://www.pacb.com/wp-content/uploads/poster_saunders.pdf
start_hap_map_block | start of hap-map block (a genomic interval in which parental phasing of reads can be performed)
end_hap_map_block | end of hap-map block 
haplotype_concordance_in_hap_map_block | the degree of similarity (on a scale of 0 to 1, with 1 representing perfect similarity) of the bit vectors: sequences of alleles observed at heterozygous sites on haplotypes inferred by read-backed vs inheritance-backed phasing 
num_het_SNVs_in_hap_map_block | number of heterozygous sites in the given hap-map block 
total_read_count_pat | number of reads from the (deduced) paternal haplotype    
total_read_count_mat | number of reads from the (deduced) maternal haplotype
founder_haplotype_pat | a label indicating which of the haplotypes in a founder of the pedigree the paternal haplotype in the given individual was descended from 
founder_haplotype_mat | a similar label corresponding to the maternal haplotype of the sample 
methylation_level_pat_count | count-based methylation level on the paternal haplotype (fraction of paternal reads that are methylated)      
methylation_level_mat_count | count-based methylation level on the maternal haplotype (fraction of maternal reads that are methylated)          
methylation_level_pat_model | model-based methylation level on the paternal haplotype (deep-learning estimated methylation using only the paternal reads)     
methylation_level_mat_model | model-based methylation level on the maternal haplotype (deep-learning estimated methylation using only the maternal reads) 
cpg_is_within_50bp_of_mismatch_site | is the given CpG site within 50bp of a heterozygous site at which the read-based and inheritance-based bit vectors (allele sequences corresponding to each haplotype) are mismatched 
cpg_overlaps_at_least_one_snv | does the CpG dinucleotide overlap an SNV? 
snv_genotypes | are the SNVs (if any) that overlap CpG dinucleotide `hom` or `het`? 
cpg_is_allele_specific | does the CpG dinucleotide appear in the reads corresponding to one haplotype, but not the reads of the other haplotype? 

## TODO

- [ ] Convert manual workflow into a Snakemake workflow (see `Snakefile`), creating a subdirectory called, e.g., `CEPH1463.GRCh38` in `read-backed-phasing` (output of `run-hiphase.sh`)
- [ ] Add cmdline invocation of the workflow to the header of the data files it outputs (https://g.co/gemini/share/594a669386ba)
