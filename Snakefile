# [WIP] complete this snakemake pipeline 

# Resources: 
# https://snakemake.readthedocs.io/en/stable/index.html
# https://github.com/quinlan-lab/Snakemake_Tutorial

# Installation Problem: 
# https://uofu.service-now.com/it?sys_id=6273777683206290ba4da250ceaad3c3&view=sp&id=ticket_history&table=incident
# https://quinlangroup.slack.com/archives/D9LFRMXV3/p1743466049592919

# Installation Resolution: 
# https://github.com/snakemake/snakemake/releases/tag/v9.2.1

# Command to run:
# snakemake --cores 16 --snakefile dna-methylation/Snakefile

input_dir = "/scratch/ucgd/lustre-core/UCGD_Datahub/Mosaic/1654/UCGD/GRCh38/LongRead/Data/PolishedCrams"
reference = "/scratch/ucgd/lustre-core/common/data/Reference/homo_sapiens/GRCh38/primary_assembly_decoy_phix.fa" 
output_dir = "/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/medgenome-model-mode-test-snakemake"
bin_dir = "/uufs/chpc.utah.edu/common/HIPAA/u6018199/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin/"

sys.path.append(f'/scratch/ucgd/lustre-labs/quinlan/u6018199/medgenome-hifi-k1375-hackathon/dna-methylation/util') 

from get_id_to_paths import get_id_to_paths_medgenome_pilot, get_uids

experiment = get_id_to_paths_medgenome_pilot()
samples = get_uids(experiment) 
for sample in samples: 
    print(sample)

rule all:
    input: 
        expand(
            "{output_dir}/{sample}.GRCh38.haplotagged.combined.bed.gz.tbi", 
            output_dir=output_dir, 
            sample=samples
        )
    shell: 
        """
        echo {input}
        """

rule compute_methylation_level: 
    log:
        # TODO: is expand necessary, to indicate that output_dir is not a wildcard
        # or could I make output_dir a snakemake param? 
        stderr = expand(
            "{output_dir}/snakemake-logs/{sample}-compute_methylation_level.log",
            output_dir=output_dir,
            sample=samples
        )
    input:
        bam = expand(
            "{input_dir}/{sample}.cram",
            input_dir=input_dir,
            sample=samples
        ),
        ref = reference
    output:
        # TODO: 
        # does this output have to match the input of "rule all"
        tbi = expand(
            "{output_dir}/{sample}.GRCh38.haplotagged.combined.bed.gz.tbi",
            output_dir=output_dir,
            sample=samples
        )
    shell:
        """
        echo "hello"
        touch {output.tbi}
        """
        # TODO: make bin_dir a snakemake param? 

        # echo $bin_dir
        # aligned_bam_to_cpg_scores --help
        # export PATH=${bin_dir}:$PATH
        # aligned_bam_to_cpg_scores \
        #     --bam {input.bam} \
        #     --ref {input.ref} \
        #     --output-prefix {output.prefix} \
        #     --threads 8 \
        #     --min-coverage 10 \
        #     --min-mapq 1 \
        #     --pileup-mode model \
        #     2> {log.stderr}
