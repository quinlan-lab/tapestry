def requiredFile(value, option) {
    if (!value) {
        error "Missing required trio argument: --${option} <path>"
    }
    def path = file(value, glob: false, checkIfExists: true)
    if (!java.nio.file.Files.isRegularFile(path)) {
        error "--${option} must be a file: ${path}"
    }
    return path
}


def findIndex(data, candidates) {
    def index = candidates.collect { name -> data.resolveSibling(name) }
        .find { path -> java.nio.file.Files.isRegularFile(path) }
    if (index != null) {
        return index
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


def annotationTarget(target) {
    // Include the G base when a CpG's C is the last target base.
    return target.contains(':') ? target.replaceFirst(/-([0-9]+)$/, { all, end -> "-${end.toLong() + 1}" }) : target
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


process SCREEN_TRIO_CROSSOVERS {
    tag "${chrom}"
    label 'trio_crossover_screening'
    publishDir "${publish_root}/per_chrom_phased", mode: 'copy', overwrite: true

    input:
    tuple val(chrom), path(recombination)
    path screener
    val distance_bp
    val publish_root

    output:
    tuple val(chrom), path("${chrom}.crossovers.tsv"), path("${chrom}.crossovers.bed"),
        emit: crossover_candidates

    script:
    """
    python3 "${screener}" --input "${recombination}" --chromosome "${chrom}" \\
        --distance-bp "${distance_bp}" --tsv "${chrom}.crossovers.tsv" \\
        --bed "${chrom}.crossovers.bed"
    """
}


process MATCH_TRIO_BIT_VECTORS {
    tag "${chrom}"
    label 'trio_bit_vector_matching'
    publishDir "${publish_root}/per_chrom_hap_maps/${chrom}", mode: 'copy', overwrite: true

    input:
    tuple val(chrom), path(vcf), path(index), path(crossovers)
    path ped, stageAs: 'trio.ped'
    path source_dir
    val distance_bp
    val publish_root

    output:
    tuple val(chrom),
        path('*.hap-map-blocks.{paternal,maternal}.bed*', arity: '6'),
        path('*.bit-vector-sites-mismatches.{paternal,maternal}.*', arity: '8'),
        path('*.hap-map.tsv', arity: '2'),
        emit: parental_maps

    script:
    """
    set -euo pipefail
    source "${source_dir}/util/logging.sh"
    source "${source_dir}/util/trio_ped.sh"
    read_trio_ids "${ped}"

    for sample in "\$kid_id" "\$dad_id" "\$mom_id"; do
        python3 -m whatshap stats --sample "\$sample" --block-list "\$sample.blocks.tsv" \\
            "${vcf}" > "\$sample.stats.tsv" 2> "\$sample.stats.log"
    done

    export PYTHONDONTWRITEBYTECODE=1
    export PYTHONPATH="${source_dir}:${source_dir}/util"
    export POLARS_MAX_THREADS="${task.cpus}"
    python3 "${source_dir}/phase_meth_to_parent_haps.py" \\
        --kid_id "\$kid_id" --dad_id "\$dad_id" --mom_id "\$mom_id" \\
        --vcf_pedmec_phased "${vcf}" --crossovers-tsv "${crossovers}" \\
        --crossover-distance-bp "${distance_bp}" \\
        --blocks_tsv_kid "\$kid_id.blocks.tsv" --blocks_tsv_dad "\$dad_id.blocks.tsv" \\
        --blocks_tsv_mom "\$mom_id.blocks.tsv" --output_dir . --skip-methylation
    """
}


process CALL_TRIO_MODEL {
    tag "${target}"
    label 'trio_model_calling'
    publishDir "${publish_root}/per_chrom_model/${chrom}", mode: 'copy', overwrite: true

    input:
    tuple val(chrom), path(vcf), path(index), val(target)
    path ped, stageAs: 'trio.ped'
    path reference_files, stageAs: 'reference/*'
    path child_files, stageAs: 'child/*'
    path father_files, stageAs: 'father/*'
    path mother_files, stageAs: 'mother/*'
    path shell_helpers
    val min_coverage
    val min_mapq
    val publish_root

    output:
    tuple val(chrom), val(target), path('model'), emit: models
    tuple val(chrom), path('haplotagged/*'), emit: haplotagged

    script:
    def reference = shellQuote(reference_files[0])
    def bam_paths = [child_files, father_files, mother_files].collect { shellQuote(it[0]) }.join(' ')
    """
    set -euo pipefail
    source ./logging.sh
    source ./trio_ped.sh
    read_trio_ids "${ped}"
    mkdir model haplotagged
    ids=("\$kid_id" "\$dad_id" "\$mom_id")
    roles=(kid dad mom)
    bams=(${bam_paths})
    for i in 0 1 2; do
        role="\${roles[\$i]}"
        # Preserve MM/ML while removing inherited HP/PS labels before re-tagging.
        samtools view -b -x HP -x PS -@ ${task.cpus} \\
            -o stripped.bam "\${bams[\$i]}" "${target}"
        samtools index -@ ${task.cpus} stripped.bam
        python3 -m whatshap haplotag --reference ${reference} \\
            --sample "\${ids[\$i]}" --output-threads ${task.cpus} \\
            --output "haplotagged/\$role.haplotagged.bam" "${vcf}" stripped.bam
        samtools index -@ ${task.cpus} "haplotagged/\$role.haplotagged.bam"
        aligned_bam_to_cpg_scores --bam "haplotagged/\$role.haplotagged.bam" \\
            --output-prefix "model/\$role" --threads ${task.cpus} \\
            --min-coverage ${min_coverage} --min-mapq ${min_mapq} --pileup-mode model
    done
    """
}

process EXPAND_TRIO_MODEL {
    tag "${target}"
    label 'trio_model_expansion'
    publishDir "${publish_root}/per_chrom_methylation/${chrom}", mode: 'copy', overwrite: true,
        saveAs: { name -> name.tokenize('/').last() }

    input:
    tuple val(chrom), val(target), path(models), path(maps), path(mismatches), path(map_tables)
    path ped, stageAs: 'trio.ped'
    path joint_files, stageAs: 'joint/*'
    path reference_files, stageAs: 'reference/*'
    path source_dir
    val publish_root

    output:
    tuple val(chrom), path('results/*.bed.gz'), path('results/*.bed.gz.tbi'),
        path('results/*.bw', arity: '9'), path('results/*.json'), emit: methylation

    script:
    def reference = shellQuote(reference_files[0])
    def joint = shellQuote(joint_files[0])
    def annotation_target = annotationTarget(target)
    """
    set -euo pipefail
    source "${source_dir}/util/logging.sh"
    source "${source_dir}/util/trio_ped.sh"
    read_trio_ids "${ped}"
    bcftools view --regions "${annotation_target}" --regions-overlap 1 \\
        -Ov -o original-trio.vcf ${joint}
    export PYTHONPATH="${source_dir}:${source_dir}/util"
    export POLARS_MAX_THREADS="${task.cpus}"
    python3 "${source_dir}/trio_model_methylation.py" --model-dir "${models}" \\
        --map-dir . --reference ${reference} --annotation-vcf original-trio.vcf \\
        --target "${target}" --kid-id "\$kid_id" --dad-id "\$dad_id" --mom-id "\$mom_id" \\
        --output-dir results
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


workflow TRIO {
    take:
    inputs

    main:
    def targets = phasingTargets(inputs.chromosomes, inputs.region)
    if (!inputs.outdir) {
        error 'Trio phasing requires --outdir <directory>'
    }
    def distance_bp = inputs.crossoverDistanceBp.toString()
    if (!(distance_bp ==~ /[0-9]+/)) {
        error '--crossover-distance-bp must be a non-negative integer'
    }
    [inputs.minCoverage, inputs.minMapq].each { setting ->
        if (!(setting.toString() ==~ /[0-9]+/)) {
            error '--min-coverage and --trio-min-mapq must be non-negative integers'
        }
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
    def screener = file("${projectDir}/src/report_whatshap_crossovers.py", checkIfExists: true)
    def shell_helpers = ['logging.sh', 'trio_ped.sh'].collect {
        file("${projectDir}/src/util/${it}", checkIfExists: true)
    }
    def publish_root = file(inputs.outdir, glob: false).toAbsolutePath()

    PHASE_TRIO_CHROMOSOME(
        Channel.fromList(targets),
        ped, joint_files, reference_files, child, father, mother,
        prefilter, shell_helpers, publish_root.toString()
    )
    SCREEN_TRIO_CROSSOVERS(
        PHASE_TRIO_CHROMOSOME.out.phased.map { chrom, vcf, index, report -> tuple(chrom, report) },
        screener, distance_bp, publish_root.toString()
    )
    matching_inputs = PHASE_TRIO_CHROMOSOME.out.phased
        .map { chrom, vcf, index, report -> tuple(chrom, vcf, index) }
        .join(
            SCREEN_TRIO_CROSSOVERS.out.crossover_candidates.map { chrom, tsv, bed -> tuple(chrom, tsv) },
            failOnDuplicate: true, failOnMismatch: true
        )
    MATCH_TRIO_BIT_VECTORS(
        matching_inputs, ped, file("${projectDir}/src", checkIfExists: true),
        distance_bp, publish_root.toString()
    )

    CALL_TRIO_MODEL(
        PHASE_TRIO_CHROMOSOME.out.phased.map { chrom, vcf, index, report -> tuple(chrom, vcf, index) }
            .join(Channel.fromList(targets), failOnDuplicate: true, failOnMismatch: true),
        ped, reference_files, child, father, mother, shell_helpers,
        inputs.minCoverage, inputs.minMapq, publish_root.toString()
    )
    EXPAND_TRIO_MODEL(
        CALL_TRIO_MODEL.out.models.join(MATCH_TRIO_BIT_VECTORS.out.parental_maps,
            failOnDuplicate: true, failOnMismatch: true),
        ped, joint_files, reference_files, file("${projectDir}/src", checkIfExists: true),
        publish_root.toString()
    )

    emit:
    phased = PHASE_TRIO_CHROMOSOME.out.phased
    crossover_candidates = SCREEN_TRIO_CROSSOVERS.out.crossover_candidates
    parental_maps = MATCH_TRIO_BIT_VECTORS.out.parental_maps
    methylation = EXPAND_TRIO_MODEL.out.methylation
    haplotagged = CALL_TRIO_MODEL.out.haplotagged
}
