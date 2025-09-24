# See umail exchange with Matt Holt on July 7th, 2025 
# Also see exchange with Tom Sasani: https://quinlangroup.slack.com/archives/D9LFRMXV3/p1751395550909439  

palladium_vcf_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/hiphase"
palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/GRCh38"
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing"
reference="/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz" 

bin_dir="/uufs/chpc.utah.edu/common/HIPAA/u6018199/hiphase-v1.5.0-x86_64-unknown-linux-gnu" 
export PATH=${bin_dir}:$PATH

# TODO: add tandem-repeat vcfs
# I've excluded them for now as the filenames follow an unknown/less-than-obvious naming scheme 

prefixes=$(python src/util/get_palladium_prefixes.py)
echo "Done loading prefixes..."

for prefix in $prefixes; do
	nohup hiphase \
		--threads 16 \
		--io-threads 4 \
		--sample-name ${prefix} \
		--vcf ${palladium_vcf_dir}/${prefix}.GRCh38.deepvariant.glnexus.phased.vcf.gz \
		--output-vcf ${output_dir}/${prefix}.GRCh38.deepvariant.glnexus.phased.vcf.gz \
		--vcf ${palladium_vcf_dir}/${prefix}.GRCh38.pbsv.phased.vcf.gz \
		--output-vcf ${output_dir}/${prefix}.GRCh38.pbsv.phased.vcf.gz \
		--bam ${palladium_bam_dir}/${prefix}.GRCh38.haplotagged.bam \
		--output-bam ${output_dir}/${prefix}.GRCh38.haplotagged.bam \
		--reference ${reference} \
		--summary-file ${output_dir}/${prefix}.GRCh38.hiphase.stats.tsv \
		--blocks-file ${output_dir}/${prefix}.GRCh38.hiphase.blocks.tsv \
		> ${output_dir}/${prefix}.log 2>&1 &
    echo "Started ${prefix}..."
done

echo "All hiphase jobs started."
echo "Check ${output_dir} for logs and outputs."
