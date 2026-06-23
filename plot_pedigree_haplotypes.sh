repo=/scratch/ucgd/lustre-labs/quinlan/u6018199/Platinum-Pedigree-Inheritance

ped=${repo}/data/CEPH1463.ped
output_dir=/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38
iht=${output_dir}/CEPH1463.GRCh38.iht.sorted.txt

chrom=chr1

out=images/tapestry.pedigree_haplotypes.CEPH1463.GRCh38.${chrom}

# "before": raw gtg-ped-map output verbatim, for side-by-side comparison so users
# can see exactly what the relabel changes (and why the trio apex looks flat).
python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
	--chrom ${chrom} \
	--no-collapse \
  	--out ${out}.raw.pdf

# "after": drop the two apex trios (G0) and relabel G1 -> A,B,C,D (the default),
# with the default vertical gap between hap1 and hap2.
python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
	--chrom ${chrom} \
	--hap-gap 0.12 \
  	--out ${out}.trios-removed.with-gap.pdf

# Same "after", but with no gap between hap1 and hap2 (the contiguous Fig 4A/B
# two-stripe look).
python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
	--chrom ${chrom} \
	--hap-gap 0 \
  	--out ${out}.trios-removed.no-gap.pdf