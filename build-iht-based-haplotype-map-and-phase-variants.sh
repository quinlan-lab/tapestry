# https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/main/HAPLOTYPING.md#building-haplotype-maps-step-one

repo=/scratch/ucgd/lustre-labs/quinlan/u6018199/Platinum-Pedigree-Inheritance

ped=${repo}/data/CEPH1463.ped
vcf=/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz
output_dir=/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38
prefix=${output_dir}/CEPH1463.GRCh38

mkdir -p ${output_dir}

# nohup ${repo}/code/rust/target/release/gtg-ped-map \
#     --ped ${ped} \
#     --vcf ${vcf} \
#     --prefix ${prefix} \
#     --verbose \
#     > ${output_dir}/gtg-ped-map.log 2>&1 &
# echo "gtg-ped-map launched"

nohup ${repo}/code/rust/target/release/gtg-concordance \
    --ped ${ped} \
    --inheritance ${prefix}.iht.txt \
    --vcf ${vcf} \
    --prefix ${prefix} \
    --verbose \
    > ${output_dir}/gtg-concordance.log 2>&1 &
echo "gtg-concordance launched"