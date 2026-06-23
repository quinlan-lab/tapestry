repo=/scratch/ucgd/lustre-labs/quinlan/u6018199/Platinum-Pedigree-Inheritance

ped=${repo}/data/CEPH1463.ped
output_dir=/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38
iht=${output_dir}/CEPH1463.GRCh38.iht.sorted.txt

out=images/tapestry.pedigree_haplotypes.CEPH1463.GRCh38.chr1

# "after": drop the two apex trios (G0) and relabel G1 -> A,B,C,D (the default).
python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
  	--out ${out}.pdf

# "before": raw gtg-ped-map output verbatim, for side-by-side comparison so users
# can see exactly what the relabel changes (and why the trio apex looks flat).
python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
	--no-collapse \
  	--out ${out}.raw.pdf