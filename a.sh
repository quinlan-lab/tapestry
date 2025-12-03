# TODO 
# 1. 

bcftools query -r chr1:146169000-146171000 -s 200081 -f '[%CHROM\t%POS\t%SAMPLE\t%REF\t%ALT\t%GT\n]' /scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/hiphase/200081.GRCh38.deepvariant.glnexus.phased.vcf.gz
bcftools query -r chr1:146169000-146171000 -s 200081 -f '[%CHROM\t%POS\t%SAMPLE\t%REF\t%ALT\t%GT\n]' /scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing/200081.GRCh38.deepvariant.glnexus.phased.vcf.gz

bcftools query -r chr1:146169000-146171000 -s 200081 -f '[%CHROM\t%POS\t%SAMPLE\t%REF\t%ALT\t%GT\n]' /scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz
bcftools query -r chr1:146169000-146171000 -s 200081 -f '[%CHROM\t%POS\t%SAMPLE\t%GT\n]' /scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38/CEPH1463.GRCh38.pass.sorted.vcf.gz
