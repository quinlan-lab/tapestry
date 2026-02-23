reference="/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz" 
vcf_trio_phased="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing/trio_phased.vcf.gz"
palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/GRCh38"
bam_kid=${palladium_bam_dir}/NA12883.GRCh38.haplotagged.bam
kid_sample_id="NA12883"
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing"

mkdir -p ${output_dir}

output_bam="${output_dir}/${kid_sample_id}.haplotagged.bam"

nohup whatshap haplotag \
    --reference ${reference} \
    --sample ${kid_sample_id} \
    --output ${output_bam} \
    --output-threads 4 \
    ${vcf_trio_phased} \
    ${bam_kid} \
	> ${output_dir}/whatshap-haplotag.${kid_sample_id}.log 2>&1 &

samtools index ${output_bam}

echo "Done. Haplotagged BAM written to ${output_bam}"
