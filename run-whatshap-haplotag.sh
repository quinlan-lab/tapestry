reference="/path/to/reference.fa"
vcf_trio_phased="/path/to/whatshap-phasing/trio_phased.vcf.gz"
bam_kid="/path/to/kid_aligned.bam"
kid_sample_id="kid"
output_dir="/path/to/whatshap-phasing"

mkdir -p ${output_dir}

output_bam="${output_dir}/${kid_sample_id}.haplotagged.bam"

whatshap haplotag \
    --reference ${reference} \
    --sample ${kid_sample_id} \
    --output ${output_bam} \
    --output-threads 4 \
    ${vcf_trio_phased} \
    ${bam_kid}

samtools index ${output_bam}

echo "Done. Haplotagged BAM written to ${output_bam}"
