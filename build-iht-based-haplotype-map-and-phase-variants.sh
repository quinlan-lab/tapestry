# https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/main/HAPLOTYPING.md#building-haplotype-maps-step-one

repo=/scratch/ucgd/lustre-labs/quinlan/u6018199/Platinum-Pedigree-Inheritance

ped=${repo}/data/CEPH1463.ped
vcf=/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz
output_dir=/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38
prefix=${output_dir}/CEPH1463.GRCh38

mkdir -p ${output_dir}

bin_dir="${repo}/code/rust/target/release" 
export PATH=${bin_dir}:$PATH

# nohup gtg-ped-map \
#     --ped ${ped} \
#     --vcf ${vcf} \
#     --prefix ${prefix} \
#     --verbose \
#     > ${output_dir}/gtg-ped-map.log 2>&1 &
# echo "gtg-ped-map launched"

nohup gtg-concordance \
    --ped ${ped} \
    --inheritance ${prefix}.iht.txt \
    --vcf ${vcf} \
    --prefix ${prefix} \
    --verbose \
    > ${output_dir}/gtg-concordance.log 2>&1 &
echo "gtg-concordance launched"

# sort and compress and index vcf 
bcftools sort ${prefix}.pass.vcf -o ${prefix}.pass.sorted.vcf
src/util/compress-index-vcf --name ${prefix}.pass.sorted

# sort and compress markers 
(head -n 1 ${prefix}.markers.txt \
  && tail -n +2 ${prefix}.markers.txt \
  | sort -k1,1V -k2,2n) \
  > ${prefix}.markers.sorted.txt

# sort and compress iht blocks 
(head -n 1 ${prefix}.iht.txt \
  && tail -n +2 ${prefix}.iht.txt \
  | sort -k1,1V -k2,2n) \
  > ${prefix}.iht.sorted.txt
