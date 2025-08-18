# c.f., /scratch/ucgd/lustre-labs/quinlan/u6018199/medgenome-hifi-k1375-hackathon/dna-methylation/aligned_bam_to_cpg_scores.sh

PB_CPG_TOOL_MODE="model" 

palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing"
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.${PB_CPG_TOOL_MODE}.read-backed-phased" 

mkdir -p ${output_dir}

bin_dir="/uufs/chpc.utah.edu/common/HIPAA/u6018199/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin" 
export PATH=${bin_dir}:$PATH

# aligned_bam_to_cpg_scores --help

echo "Getting Palladium prefixes..."
prefixes=$(python src/util/get_palladium_prefixes.py)
echo "... Got prefixes" 

for prefix in $prefixes; do
    nohup aligned_bam_to_cpg_scores \
        --bam ${palladium_bam_dir}/${prefix}.GRCh38.haplotagged.bam \
        --output-prefix ${output_dir}/${prefix}.GRCh38.haplotagged \
        --threads 8 \
        --min-coverage 10 \
        --min-mapq 1 \
        --pileup-mode ${PB_CPG_TOOL_MODE} \
        > ${output_dir}/${prefix}.2.log 2>&1 &
    echo "Started ${prefix}"
done

