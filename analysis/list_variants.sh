# TODO 

SAMPLE="$1"
REGION="$2"

SINGLE_SAMPLE_VCF="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/hiphase/${SAMPLE}.GRCh38.deepvariant.glnexus.phased.vcf.gz"
SINGLE_SAMPLE_VCF_READ_PHASED="/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing/${SAMPLE}.GRCh38.deepvariant.glnexus.phased.vcf.gz" 

JOINT_CALLED_VCF="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
JOINT_CALLED_VCF_IHT_PHASED="/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38/CEPH1463.GRCh38.pass.sorted.vcf.gz"

FORMAT='[%CHROM\t%POS\t%REF\t%ALT\t%GT\n]'

echo "Variants from single-sample VCF:"
bcftools query -s ${SAMPLE} -r ${REGION} -f ${FORMAT} ${SINGLE_SAMPLE_VCF}
echo ""

echo "Variants from single-sample VCF (read-phased):"
bcftools query -s ${SAMPLE} -r ${REGION} -f ${FORMAT} ${SINGLE_SAMPLE_VCF_READ_PHASED}
echo ""

echo "Variants from joint-called VCF:"
bcftools query -s ${SAMPLE} -r ${REGION} -f ${FORMAT} ${JOINT_CALLED_VCF}
echo ""

echo "Variants from joint-called VCF (inheritance-phased):"
bcftools query -s ${SAMPLE} -r ${REGION} -f ${FORMAT} ${JOINT_CALLED_VCF_IHT_PHASED}
