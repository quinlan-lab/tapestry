repo=/scratch/ucgd/lustre-labs/quinlan/u6018199/Platinum-Pedigree-Inheritance

ped=${repo}/data/CEPH1463.ped
output_dir=/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38
iht=${output_dir}/CEPH1463.GRCh38.iht.sorted.txt

# The default chrom/start_mb/end_mb below (chr22, 10-50 Mb) are chosen so the
# output can be validated by manual inspection against Zev Kronenberg's pedigree
# haplotype visualization of the same region:
#   https://docs.google.com/document/d/1HGBhjqk4cqZaWIkh0GeAaBC_fcZYWRRnYO5cjiifKEY/edit?usp=sharing
# The same visualization also appears in the Platinum-Pedigree-Inheritance docs:
#   https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/e12aca6b49ee7208952467db4a2a9e2f79b98efb/HAPLOTYPING.md
chrom=chr22

# Optional region to visualize, in Mb. Leave both empty to plot the whole
# chromosome (no clipping), exactly as the script did before a region was added.
# Set both to zoom in; blocks are clipped to the window and the Mb range is
# recorded in the filename.
start_mb=10
end_mb=50

out=images/tapestry.pedigree_haplotypes.CEPH1463.GRCh38.${chrom}
region_args=""
if [ -n "${start_mb}" ] && [ -n "${end_mb}" ]; then
	start_bp=$(awk "BEGIN{printf \"%d\", ${start_mb}*1000000}")
	end_bp=$(awk "BEGIN{printf \"%d\", ${end_mb}*1000000}")
	region_args="--start ${start_bp} --end ${end_bp}"
	out=${out}.${start_mb}-${end_mb}Mb
fi

# "before": raw gtg-ped-map output verbatim, for side-by-side comparison so users
# can see exactly what the relabel changes (and why the trio apex looks flat).
python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
	--chrom ${chrom} \
	${region_args} \
	--no-collapse \
  	--out ${out}.raw.pdf

# "after": drop the two apex trios (G0) and relabel G1 -> A,B,C,D (the default),
# with the default vertical gap between hap1 and hap2.
python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
	--chrom ${chrom} \
	${region_args} \
	--hap-gap 0.12 \
  	--out ${out}.trios-removed.with-gap.pdf

# Same "after", but with no gap between hap1 and hap2 (the contiguous Fig 4A/B
# two-stripe look).
python src/plot_pedigree_haplotypes.py \
	--ped ${ped} \
	--iht ${iht} \
	--chrom ${chrom} \
	${region_args} \
	--hap-gap 0 \
  	--out ${out}.trios-removed.no-gap.pdf