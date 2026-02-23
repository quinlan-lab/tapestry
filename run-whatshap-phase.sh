trio_ped="trio.ped"
reference="/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz" 
vcf_joint_called="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/GRCh38"
bam_kid=${palladium_bam_dir}/NA12883.GRCh38.haplotagged.bam
bam_dad=${palladium_bam_dir}/NA12877.GRCh38.haplotagged.bam
bam_mom=${palladium_bam_dir}/NA12878.GRCh38.haplotagged.bam
kid_sample_id="NA12883"
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing/"

mkdir -p ${output_dir}

# Step 1a: unphase the input VCF to avoid mixed phasing
# (whatshap phase --ped copies over variants it cannot phase, so if the input
# is already phased those variants would retain their original phase assignments
# rather than being assigned consistent pedigree-based phase)
vcf_unphased="${output_dir}/trio_unphased.vcf"
whatshap unphase ${vcf_joint_called} > ${vcf_unphased}

# Step 1b: pedigree-aware phasing
output_vcf="${output_dir}/trio_phased.vcf"

nohup whatshap phase \
    --ped ${trio_ped} \
    --sample ${kid_sample_id} \
    --use-ped-samples \
    --reference ${reference} \
    --output ${output_vcf} \
    ${vcf_unphased} \
    ${bam_kid} ${bam_dad} ${bam_mom} \
	> ${output_dir}/whatshap-phase.${kid_sample_id}.log 2>&1 &

bgzip ${output_vcf} && tabix ${output_vcf}.gz

echo "Done. Phased VCF written to ${output_vcf}.gz"
