#!/bin/bash

# Usage:
#   Production:  ./run-whatshap.sh
#   Dev mode:    ./run-whatshap.sh --dev-dir trio_dev_data

DEV_DIR=""

# Argument handling
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dev-dir|-d)
            DEV_DIR="${2%/}"
            shift 2
            ;;
        *)
            echo "Error: Unknown parameter passed: $1"
            echo "Usage: $0 [--dev-dir <DEV_DATA_DIR>]"
            exit 1
            ;;
    esac
done

source src/util/logging.sh 

# --- Default Configurations (Production) ---
trio_ped="/scratch/ucgd/lustre-labs/quinlan/u6018199/tapestry/trio_dev_data/input/trio.ped"
reference="/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz" 
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/pedmec-phasing"

# https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets?tab=readme-ov-file#accessing-controlled-samples
kid_id="NA12878"
dad_id="NA12891" 
mom_id="NA12892"

vcf_joint_called="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/GRCh38"

bam_kid="${palladium_bam_dir}/${kid_id}.GRCh38.haplotagged.bam" # not topped off
bam_dad="${palladium_bam_dir}/${dad_id}.GRCh38.haplotagged.bam" # not topped off
bam_mom="${palladium_bam_dir}/${mom_id}.GRCh38.haplotagged.bam" # not topped off

# --- Optional Dev Data Overrides ---
if [ -n "$DEV_DIR" ]; then
    log_info "DEV MODE ENABLED: Reading from and writing to ${DEV_DIR}"
    
	# Override input paths to point to the generated dev data
    trio_ped="${DEV_DIR}/input/trio.ped"
    vcf_joint_called="${DEV_DIR}/input/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz"
	bam_kid="${DEV_DIR}/input/${kid_id}.GRCh38.haplotagged.bam"
    bam_dad="${DEV_DIR}/input/${dad_id}.GRCh38.haplotagged.bam"
    bam_mom="${DEV_DIR}/input/${mom_id}.GRCh38.haplotagged.bam"
    
    source "$(dirname "$0")/trio_dev_data_config.sh"
    dev_chrom="${DEV_REGION%%:*}"
    reference="${DEV_DIR}/input/${dev_chrom}.fa"

    # Override output dir to write entirely inside the dev directory
    # (Appending '/output' to keep the generated files separate from the raw dev inputs)
    output_dir="${DEV_DIR}/output/pedmec-phasing"
    
    # In dev mode, only phase the single dev chromosome
    dev_chromosomes=("${dev_chrom}")

    # Fail fast if dev data doesn't exist
    if [ ! -f "$vcf_joint_called" ]; then
        echo "Error: Dev VCF not found. Did you run trio_dev_data_create.sh?"
        exit 1
    fi
fi

# Create the output directory (whether production or dev)
mkdir -p ${output_dir}

get_bam() {
    local id=$1
    if [ "$id" == "$kid_id" ]; then
        echo "$bam_kid"
    elif [ "$id" == "$dad_id" ]; then
        echo "$bam_dad"
    elif [ "$id" == "$mom_id" ]; then
        echo "$bam_mom"
    else
        echo "Error: Unknown sample ID: $id" >&2
        exit 1
    fi
}

# --- Main Workflow ---

log_info "Unphasing: '${vcf_joint_called}'" 

# Step 1: unphase the input VCF to avoid mixed phasing
vcf_joint_called_unphased="${output_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.unphased.vcf"
# whatshap unphase ${vcf_joint_called} \
#     > ${vcf_joint_called_unphased} 

# bgzip -f ${vcf_joint_called_unphased}
# tabix ${vcf_joint_called_unphased}.gz

# Step 2: pedigree-aware phasing (parallelized by chromosome)
vcf_joint_called_phased="${output_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
per_chrom_dir="${output_dir}/per_chrom_phased"
mkdir -p ${per_chrom_dir}

if [ ${#dev_chromosomes[@]} -gt 0 ] 2>/dev/null; then
    chromosomes=("${dev_chromosomes[@]}")
else
    chromosomes=(chr{1..22} chrX chrY chrM)
fi
PHASE_THREADS=13
HAPLOTAG_THREADS=24

log_info "pedMEC Phasing: ${kid_id} ${dad_id} ${mom_id}"
log_info "Phasing ${#chromosomes[@]} chromosomes in parallel (${PHASE_THREADS} at a time) in: '${per_chrom_dir}'"

# https://whatshap.readthedocs.io/en/latest/guide.html#phasing-pedigrees
phase_chrom() {
    local chrom=$1
    local chrom_input_vcf="${per_chrom_dir}/${chrom}.unphased.vcf.gz"
    local chrom_vcf="${per_chrom_dir}/${chrom}.phased.vcf.gz"
    local chrom_log="${per_chrom_dir}/${chrom}.phased.log"

    # Pre-split: extract this chromosome from the unphased VCF
    bcftools view -r ${chrom} ${vcf_joint_called_unphased}.gz -Oz -o ${chrom_input_vcf}
    tabix ${chrom_input_vcf}

    whatshap phase \
        --ped ${trio_ped} \
        --sample ${kid_id} \
        --sample ${dad_id} \
        --sample ${mom_id} \
        --reference ${reference} \
        --output ${chrom_vcf} \
        ${chrom_input_vcf} \
        ${bam_kid} ${bam_dad} ${bam_mom} \
        > ${chrom_log} 2>&1
}
export -f phase_chrom
export trio_ped kid_id dad_id mom_id reference vcf_joint_called_unphased bam_kid bam_dad bam_mom per_chrom_dir

# printf '%s\n' "${chromosomes[@]}" | xargs -P ${PHASE_THREADS} -I {} bash -c 'phase_chrom "$@"' _ {}

log_info "All per-chromosome phasing complete. Merging..."

# Index per-chromosome VCFs, then merge with --naive (no decompression since chromosomes are non-overlapping)
# chrom_vcfs=()
# for chrom in "${chromosomes[@]}"; do
#     chrom_vcf="${per_chrom_dir}/${chrom}.phased.vcf.gz"
#     tabix -f ${chrom_vcf}
#     chrom_vcfs+=("${chrom_vcf}")
# done

# bcftools concat --naive "${chrom_vcfs[@]}" -Oz -o ${vcf_joint_called_phased}

log_info "Indexing genome-wide phased VCF: '${vcf_joint_called_phased}'"
# tabix ${vcf_joint_called_phased}

log_info "Creating phasing statistics ..." 

# for id in ${kid_id} ${dad_id} ${mom_id}; do
# 	stats_log="${output_dir}/${id}.stats.log"
# 	blocks_gtf="${output_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${id}.blocks.gtf"
# 	blocks_tsv="${output_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${id}.blocks.tsv"
# 	log_info "Running whatshap stats for ${id}"
# 	whatshap stats \
# 		--gtf ${blocks_gtf} \
# 		--block-list ${blocks_tsv} \
# 		--sample ${id} \
# 		${vcf_joint_called_phased} \
# 		> ${stats_log} 2>&1
# 	log_info "Done whatshap stats for ${id}"
# done

log_info "Created haplotype blocks from: '${vcf_joint_called_phased}'"

# Step 3: haplotag the bams (parallelized by sample)
log_info "Haplotagging 3 samples in parallel (${HAPLOTAG_THREADS} threads each)"

haplotag_sample() {
    local id=$1

    local input_bam
    if [ "$id" == "${kid_id}" ]; then
        input_bam="${bam_kid}"
    elif [ "$id" == "${dad_id}" ]; then
        input_bam="${bam_dad}"
    elif [ "$id" == "${mom_id}" ]; then
        input_bam="${bam_mom}"
    fi

    local stripped_bam="${output_dir}/${id}.GRCh38.stripped.bam"
    log_info "Stripping HP/PS tags from BAM for ${id}"
    # samtools view -b -x HP -x PS ${input_bam} -o ${stripped_bam}
    # samtools index ${stripped_bam}
    log_info "Done stripping HP/PS tags from BAM for ${id}"

    local output_bam="${output_dir}/${id}.GRCh38.haplotagged.bam"
    local haplotag_stats_log="${output_dir}/${id}.haplotag.stats.log"
    python -c "
import hashlib, functools
hashlib.md5 = functools.partial(hashlib.md5, usedforsecurity=False)
from whatshap.__main__ import main
main(argv=[
    'haplotag',
    '--reference', '${reference}',
    '--sample', '${id}',
    '--output', '${output_bam}',
    '--output-threads', '${HAPLOTAG_THREADS}',
    '${vcf_joint_called_phased}',
    '${stripped_bam}',
])" > ${haplotag_stats_log} 2>&1
    samtools index ${output_bam}
}
export -f haplotag_sample log_info log
export kid_id dad_id mom_id bam_kid bam_dad bam_mom reference output_dir vcf_joint_called_phased HAPLOTAG_THREADS

# printf '%s\n' ${kid_id} ${dad_id} ${mom_id} | xargs -P 3 -I {} bash -c 'haplotag_sample "$@" > "${output_dir}/$1.haplotag.log" 2>&1' _ {}
printf '%s\n' ${kid_id} | xargs -P 3 -I {} bash -c 'haplotag_sample "$@" > "${output_dir}/$1.haplotag.log" 2>&1' _ {}

log_info "All samples haplotagged"

log_info "Done: run-whatshap.sh"