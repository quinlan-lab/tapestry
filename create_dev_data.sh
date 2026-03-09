#!/bin/bash

# Check if an output directory was provided
if [ "$#" -ne 1 ]; then
    echo "Error: Target directory required."
    echo "Usage: $0 <DEV_DATA_OUTPUT_DIR>"
    exit 1
fi

DEV_DIR=$1
mkdir -p "${DEV_DIR}/dev_bams"

# --- Source Configurations ---
vcf_joint_called="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/GRCh38"

kid_id="NA12878"
dad_id="NA12891" 
mom_id="NA12892"

DEV_REGION="chr22:1-10700000"

echo "Creating dev data in ${DEV_DIR} for region ${DEV_REGION}..."

# Create PED file
dev_ped="${DEV_DIR}/trio.ped"
echo "Creating pedigree file..."
cat << 'EOF' > ${dev_ped}
# https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets?tab=readme-ov-file#accessing-controlled-samples
1463	NA12878	NA12891	NA12892	2	UW
EOF

# Subset VCF
dev_vcf="${DEV_DIR}/dev_input.vcf.gz"
echo "Subsetting VCF..."
bcftools view -r ${DEV_REGION} -s ${kid_id},${dad_id},${mom_id} ${vcf_joint_called} -Oz -o ${dev_vcf}
bcftools index ${dev_vcf}

# Subset BAMs function
subset_bam() {
    local in_bam=$1
    local out_bam="${DEV_DIR}/dev_bams/$(basename $in_bam)"
    echo "Subsetting ${in_bam}..."
    samtools view -b ${in_bam} ${DEV_REGION} > ${out_bam}
    samtools index ${out_bam}
}

subset_bam "${palladium_bam_dir}/${kid_id}.GRCh38.haplotagged.bam"
subset_bam "${palladium_bam_dir}/${dad_id}.GRCh38.haplotagged.bam"
subset_bam "${palladium_bam_dir}/${mom_id}.GRCh38.haplotagged.bam"

# Subset reference FASTA from UCSC
dev_ref="${DEV_DIR}/dev_reference.fa"
dev_chrom="${DEV_REGION%%:*}"
echo "Fetching reference for ${dev_chrom} from UCSC..."
curl -sS "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${dev_chrom}.fa.gz" | gunzip > "${dev_ref}"
samtools faidx "${dev_ref}"

echo "Dev data creation complete! You can now ingest this with the workflow script."