#!/bin/bash

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
trio_ped="trio.ped"
reference="/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz" 
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing"

# https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets?tab=readme-ov-file#accessing-controlled-samples
kid_id="NA12878"
dad_id="NA12891" 
mom_id="NA12892"

vcf_joint_called="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/GRCh38"

bam_kid="${palladium_bam_dir}/${kid_id}.GRCh38.haplotagged.bam"
bam_dad="${palladium_bam_dir}/${dad_id}.GRCh38.haplotagged.bam"
bam_mom="${palladium_bam_dir}/${mom_id}.GRCh38.haplotagged.bam"

# --- Optional Dev Data Overrides ---
if [ -n "$DEV_DIR" ]; then
    log_info "DEV MODE ENABLED: Reading from and writing to ${DEV_DIR}"
    
	# Override input paths to point to the generated dev data
    trio_ped="${DEV_DIR}/input/trio.ped"
    vcf_joint_called="${DEV_DIR}/input/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz"
	bam_kid="${DEV_DIR}/input/${kid_id}.GRCh38.haplotagged.bam"
    bam_dad="${DEV_DIR}/input/${dad_id}.GRCh38.haplotagged.bam"
    bam_mom="${DEV_DIR}/input/${mom_id}.GRCh38.haplotagged.bam"
    
    reference="${DEV_DIR}/input/dev_reference.fa"

    # Override output dir to write entirely inside the dev directory
    # (Appending '/output' to keep the generated files separate from the raw dev inputs)
    output_dir="${DEV_DIR}/output"
    
    # Fail fast if dev data doesn't exist
    if [ ! -f "$vcf_joint_called" ]; then
        echo "Error: Dev VCF not found. Did you run create_dev_data.sh?"
        exit 1
    fi
fi

# Create the output directory (whether production or dev)
mkdir -p ${output_dir}

# --- Main Workflow ---

log_info "Unphasing: '${vcf_joint_called}'" 

# Step 1: unphase the input VCF to avoid mixed phasing
vcf_joint_called_unphased="${output_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.unphased.vcf"
whatshap unphase ${vcf_joint_called} \
    > ${vcf_joint_called_unphased} \
    2> ${output_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.unphased.log 

bgzip -f ${vcf_joint_called_unphased}
tabix ${vcf_joint_called_unphased}.gz

# Step 2: pedigree-aware phasing
vcf_joint_called_phased="${output_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"

log_info "Phasing: ${kid_id}" 

whatshap phase \
    --ped ${trio_ped} \
    --sample ${kid_id} \
    --sample ${dad_id} \
    --sample ${mom_id} \
    --reference ${reference} \
    --output ${vcf_joint_called_phased} \
    ${vcf_joint_called_unphased}.gz \
    ${bam_kid} ${bam_dad} ${bam_mom} \
    2> ${output_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.${kid_id}.phased.log

log_info "Indexing: '${vcf_joint_called_phased}'" 
tabix ${vcf_joint_called_phased}

log_info "${kid_id} phased in: '${vcf_joint_called_phased}'"