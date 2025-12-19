sample=$1
locus=$2
bed="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased.all-cpgs/${sample}.dna-methylation.founder-phased.all_cpgs.sorted.bed.gz"

tabix -h ${bed} ${locus} \
	| tail -n +2 \
	| cut -f1-3,4-5,11-12,15-16,20-22 \
	| column -s $'\t' -t