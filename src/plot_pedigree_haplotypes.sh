repo=/scratch/ucgd/lustre-labs/quinlan/u6018199/Platinum-Pedigree-Inheritance

ped=${repo}/data/CEPH1463.ped
output_dir=/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38
iht=${output_dir}/CEPH1463.GRCh38.iht.sorted.txt

python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
  	--out images/tapestry.pedigree_haplotypes.CEPH1463.GRCh38.chr1.pdf