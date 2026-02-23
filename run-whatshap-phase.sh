trio_ped="trio.ped"
reference="/path/to/reference.fa"
vcf_joint_called="/path/to/joint_called_trio.vcf.gz"
bam_kid="/path/to/kid_aligned.bam"
bam_dad="/path/to/dad_aligned.bam"
bam_mom="/path/to/mom_aligned.bam"
kid_sample_id="kid"
output_dir="/path/to/whatshap-phasing"

mkdir -p ${output_dir}

# Step 1a: unphase the input VCF to avoid mixed phasing
# (whatshap phase --ped copies over variants it cannot phase, so if the input
# is already phased those variants would retain their original phase assignments
# rather than being assigned consistent pedigree-based phase)
vcf_unphased="${output_dir}/trio_unphased.vcf"
whatshap unphase ${vcf_joint_called} > ${vcf_unphased}

# Step 1b: pedigree-aware phasing
output_vcf="${output_dir}/trio_phased.vcf"

whatshap phase \
    --ped ${trio_ped} \
    --sample ${kid_sample_id} \
    --use-ped-samples \
    --reference ${reference} \
    --output ${output_vcf} \
    ${vcf_unphased} \
    ${bam_kid} ${bam_dad} ${bam_mom}

bgzip ${output_vcf} && tabix ${output_vcf}.gz

echo "Done. Phased VCF written to ${output_vcf}.gz"
