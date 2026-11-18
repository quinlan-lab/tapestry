import java.nio.file.Files


def requiredFile(value, option) {
    if (!value) {
        error "Missing required trio argument: --${option} <path>"
    }
    def path = file(value, glob: false, checkIfExists: true)
    if (!Files.isRegularFile(path)) {
        error "--${option} must be a file: ${path}"
    }
    return path
}


def findIndex(data, candidates) {
    for (def name : candidates) {
        def index = data.resolveSibling(name)
        if (Files.isRegularFile(index)) {
            return index
        }
    }
    error "Missing index for ${data}. Expected one of: ${candidates.join(', ')}"
}


def indexedBam(value, option) {
    def bam = requiredFile(value, option)
    def name = bam.name
    if (!name.endsWith('.bam')) {
        error "--${option} must be a BAM file: ${bam}"
    }
    def stem = name - ~/\.bam$/
    return [bam, findIndex(bam, ["${name}.bai", "${stem}.bai", "${name}.csi", "${stem}.csi"])]
}


def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}


process PHASE_TRIO_CHROMOSOME {
    tag "${target}"
    label 'trio_phasing'
    publishDir "${publish_root}/per_chrom_phased", mode: 'copy', overwrite: true

    input:
    tuple val(chrom), val(target)
    path ped, stageAs: 'trio.ped'
    path joint_files, stageAs: 'joint/*'
    path reference_files, stageAs: 'reference/*'
    path child_files, stageAs: 'child/*'
    path father_files, stageAs: 'father/*'
    path mother_files, stageAs: 'mother/*'
    path prefilter
    path shell_helpers
    val publish_root

    output:
    tuple val(chrom), path("${chrom}.phased.vcf.gz"), path("${chrom}.phased.vcf.gz.tbi"),
        path("${chrom}.recombination.txt"), emit: phased

    script:
    def joint = shellQuote(joint_files[0])
    def reference = shellQuote(reference_files[0])
    def bams = [child_files, father_files, mother_files].collect { shellQuote(it[0]) }.join(' ')
    """
    set -euo pipefail
    source ./logging.sh
    source ./trio_ped.sh
    read_trio_ids "${ped}"

    bcftools view --regions "${target}" --regions-overlap 0 \\
        -Oz -o selected.vcf.gz ${joint}
    python3 "${prefilter}" --input selected.vcf.gz --output filtered.vcf.gz \\
        --samples "\$kid_id" "\$dad_id" "\$mom_id"
    python3 -m whatshap unphase filtered.vcf.gz | bgzip -c > unphased.vcf.gz
    bcftools index --tbi unphased.vcf.gz

    python3 -m whatshap phase \\
        --ped "${ped}" --reference ${reference} --chromosome "${chrom}" \\
        --sample "\$kid_id" --sample "\$dad_id" --sample "\$mom_id" \\
        --recombination-list "${chrom}.recombination.txt" \\
        --output "${chrom}.phased.vcf.gz" unphased.vcf.gz ${bams}
    bcftools index --tbi "${chrom}.phased.vcf.gz"
    """
}


def phasingTargets(chromosomeList, region) {
    if (!chromosomeList) {
        error 'Trio phasing requires --chromosomes <chr1,chr2,...>'
    }
    def chromosomes = chromosomeList.toString().split(',', -1).collect { it.trim() }
    if (!chromosomes.every { it ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/ } ||
            chromosomes.size() != chromosomes.toSet().size()) {
        error '--chromosomes must be a comma-separated list of unique chromosome names'
    }
    if (region) {
        def interval = region.toString() =~ /^([^:]+):([1-9][0-9]*)-([1-9][0-9]*)$/
        if (!interval.matches() || chromosomes != [interval[0][1]] ||
                interval[0][2].toBigInteger() > interval[0][3].toBigInteger()) {
            error '--region must be CHROM:START-END with one-based inclusive coordinates on the only requested chromosome'
        }
        return [tuple(chromosomes.first(), region.toString())]
    }
    def autosomes = chromosomes.findAll { it ==~ /chr[0-9]+/ }
    def other_chromosomes = chromosomes - autosomes
    return (autosomes.sort { it.substring(3).toBigInteger() } + other_chromosomes.sort())
        .collect { chrom -> tuple(chrom, chrom) }
}


workflow TRIO_PHASING {
    take:
    inputs

    main:
    def targets = phasingTargets(inputs.chromosomes, inputs.region)
    if (!inputs.outdir) {
        error 'Trio phasing requires --outdir <directory>'
    }
    def ped = requiredFile(inputs.ped, 'ped')
    def joint = requiredFile(inputs.jointVcf, 'joint-vcf')
    def joint_files = [joint, findIndex(joint, ["${joint.name}.csi", "${joint.name}.tbi"])]
    def reference = requiredFile(inputs.reference, 'reference-fasta')
    def reference_files = [reference, findIndex(reference, ["${reference.name}.fai"])]
    if (reference.name.endsWith('.gz') || reference.name.endsWith('.bgz')) {
        reference_files.add(findIndex(reference, ["${reference.name}.gzi"]))
    }
    def child = indexedBam(inputs.childBam, 'child-bam')
    def father = indexedBam(inputs.fatherBam, 'father-bam')
    def mother = indexedBam(inputs.motherBam, 'mother-bam')
    def prefilter = file("${projectDir}/src/prefilter_whatshap_vcf.py", checkIfExists: true)
    def shell_helpers = ['logging.sh', 'trio_ped.sh'].collect {
        file("${projectDir}/src/util/${it}", checkIfExists: true)
    }
    def publish_root = file(inputs.outdir, glob: false).toAbsolutePath()

    PHASE_TRIO_CHROMOSOME(
        Channel.fromList(targets),
        ped, joint_files, reference_files, child, father, mother,
        prefilter, shell_helpers, publish_root.toString()
    )

    emit:
    phased = PHASE_TRIO_CHROMOSOME.out.phased
}
