import java.nio.file.Paths
import org.yaml.snakeyaml.Yaml

nextflow.enable.dsl = 2

params['run-config'] = null
params.container = null


def resolvedPath(baseDir, value) {
    def path = Paths.get(value.toString())
    return (path.isAbsolute() ? path : baseDir.resolve(path)).normalize().toAbsolutePath()
}


def collectManifestMounts(node, key, baseDir, roots) {
    def pathKeys = [
        'vcf', 'index', 'bed', 'bam', 'phase_blocks', 'phase_stats',
        'phase_haplotags', 'outputs_manifest', 'metadata'
    ] as Set
    if (node instanceof Map) {
        node.each { childKey, childValue ->
            collectManifestMounts(childValue, childKey.toString(), baseDir, roots)
        }
    }
    else if (node instanceof Collection) {
        node.each { child -> collectManifestMounts(child, key, baseDir, roots) }
    }
    else if (node instanceof String && pathKeys.contains(key)) {
        roots.add(resolvedPath(baseDir, node).parent)
    }
}


def validationMountOptions(configFile, config, containerEngine) {
    def roots = new LinkedHashSet()
    roots.add(configFile.parent)
    try {
        ['ped'].each { key ->
            if (config.pedigree?.get(key)) {
                roots.add(resolvedPath(configFile.parent, config.pedigree[key]).parent)
            }
        }
        ['fasta', 'fasta_index', 'fasta_gzi'].each { key ->
            if (config.reference?.get(key)) {
                roots.add(resolvedPath(configFile.parent, config.reference[key]).parent)
            }
        }
        if (config.upstream?.manifest) {
            def manifestPath = resolvedPath(configFile.parent, config.upstream.manifest)
            roots.add(manifestPath.parent)
            if (manifestPath.toFile().isFile()) {
                def manifest = new groovy.json.JsonSlurper().parse(manifestPath.toFile())
                collectManifestMounts(manifest, '', manifestPath.parent, roots)
            }
        }
    }
    catch (Exception ignored) {
        // The Python validator reports schema/path errors with precise field context.
    }
    if (containerEngine == 'docker') {
        return roots.collect { root ->
            "--mount type=bind,source=${root},target=${root},readonly"
        }.join(' ')
    }
    if (containerEngine in ['apptainer', 'singularity']) {
        return roots.collect { root ->
            "--bind ${root}:${root}:ro"
        }.join(' ')
    }
    return ''
}


def resolveRunSettings() {
    def runConfig = params['run-config']
    if (!runConfig) {
        error "Missing required pipeline argument: --run-config <run.yaml|run.json>"
    }

    def configFile = file(runConfig, checkIfExists: true).toAbsolutePath()
    def config
    if (configFile.name.toLowerCase().endsWith('.json')) {
        config = new groovy.json.JsonSlurper().parseText(configFile.text)
    }
    else {
        config = new Yaml().load(configFile.text)
    }
    if (!(config instanceof Map) || !(config.project instanceof Map) || !config.project.outdir) {
        error "Run config must define project.outdir; full validation is performed by VALIDATE_INPUTS"
    }

    def configuredOutdir = Paths.get(config.project.outdir.toString())
    if (!configuredOutdir.isAbsolute()) {
        configuredOutdir = configFile.parent.resolve(configuredOutdir)
    }
    return [
        configFile.toString(),
        configuredOutdir.normalize().toAbsolutePath().toString(),
        validationMountOptions(
            configFile,
            config,
            workflow.containerEngine?.toString()
        )
    ]
}


process VALIDATE_INPUTS {
    tag 'run contract'
    containerOptions { container_mounts }
    publishDir "${publish_root}/pipeline_info", mode: 'copy', overwrite: true

    input:
    env TAPESTRY_RUN_CONFIG
    val publish_root
    path validator
    path schema_dir
    val container_mounts

    output:
    stdout emit: validation_summary
    path 'resolved-run.json', emit: resolved_run
    path 'resolved-manifest.json', emit: resolved_manifest
    path 'normalized.ped', emit: normalized_ped
    path 'selected-samples.tsv', emit: selected_samples
    path 'selected-artifacts.json', emit: selected_artifacts
    path 'validation-report.json', emit: validation_report
    path 'config-fingerprint.txt', emit: config_fingerprint
    path 'validation.success', emit: validation_success

    script:
    """
    set -euo pipefail
    python3 "${validator}" \\
      --run-config "\${TAPESTRY_RUN_CONFIG}" \\
      --output-dir . \\
      --schema-dir "${schema_dir}"
    """
}


process FILTER_MODEL_BEDS {
    tag "${sample_id}"
    publishDir "${publish_root}/samples/${sample_id}/model_inputs", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id),
          path(combined_bed), path(combined_index),
          path(hap1_bed), path(hap1_index),
          path(hap2_bed), path(hap2_index),
          val(min_coverage), val(regions)
    path validation_success
    val publish_root
    path model_filter

    output:
    tuple val(sample_id),
          path('*.model.combined.bed.gz'), path('*.model.combined.bed.gz.tbi'),
          path('*.model.hap1.bed.gz'), path('*.model.hap1.bed.gz.tbi'),
          path('*.model.hap2.bed.gz'), path('*.model.hap2.bed.gz.tbi'),
          path('*.model-filter-qc.json'),
          emit: filtered_models

    script:
    """
    set -euo pipefail
    test -s validation.success
    python3 "${model_filter}" \\
      --sample-id "${sample_id}" \\
      --combined-bed "${combined_bed}" \\
      --combined-index "${combined_index}" \\
      --hap1-bed "${hap1_bed}" \\
      --hap1-index "${hap1_index}" \\
      --hap2-bed "${hap2_bed}" \\
      --hap2-index "${hap2_index}" \\
      --min-coverage "${min_coverage}" \\
      --regions "${regions}" \\
      --output-dir .
    """
}


process CAPTURE_RUNTIME_VERSIONS {
    tag 'runtime provenance'
    publishDir "${publish_root}/pipeline_info", mode: 'copy', overwrite: true

    input:
    path resolved_manifest
    path validation_success
    val publish_root
    path version_capture

    output:
    path 'versions.json', emit: versions

    script:
    """
    set -euo pipefail
    test -s validation.success
    python3 "${version_capture}" \\
      --resolved-manifest "${resolved_manifest}" \\
      --container "${task.container ?: 'none'}" \\
      --nextflow-version "${workflow.nextflow.version}" \\
      --output versions.json
    """
}


process NORMALIZE_JOINT_VCF {
    tag 'joint VCF'
    containerOptions { container_mounts }
    publishDir "${publish_root}/inheritance", mode: 'copy', overwrite: true

    input:
    path resolved_run
    path resolved_manifest
    path normalized_ped
    path validation_success
    val publish_root
    path normalizer
    val container_mounts

    output:
    path '*.all-sites.vcf.gz', emit: all_sites_vcf
    path '*.all-sites.vcf.gz.tbi', emit: all_sites_index
    path '*.complete-family-map.vcf.gz', emit: map_vcf
    path '*.complete-family-map.vcf.gz.tbi', emit: map_index
    path '*.vcf-normalization.json', emit: normalization_report

    script:
    """
    set -euo pipefail
    test -s validation.success
    python3 "${normalizer}" \\
      --resolved-run "${resolved_run}" \\
      --resolved-manifest "${resolved_manifest}" \\
      --ped "${normalized_ped}" \\
      --output-dir .
    """
}


process RUN_GTG_INHERITANCE {
    tag 'gtg map and concordance'
    publishDir "${publish_root}/inheritance", mode: 'copy', overwrite: true

    input:
    path resolved_run
    path normalized_ped
    path all_sites_vcf
    path all_sites_index
    path map_vcf
    path map_index
    path normalization_report
    path validation_success
    val publish_root
    path gtg_runner

    output:
    path '*.iht.sorted.txt', emit: iht
    path '*.markers.sorted.txt', emit: markers
    path '*.recombinants.sorted.txt', emit: recombinants
    path '*.pass.vcf.gz', emit: pass_vcf
    path '*.pass.vcf.gz.tbi', emit: pass_index
    path '*.fail.vcf.gz', emit: fail_vcf
    path '*.fail.vcf.gz.tbi', emit: fail_index
    path '*.filtering-stats.txt', emit: filtering_stats
    path '*.inheritance-qc.json', emit: inheritance_qc
    path '*.gtg-ped-map.log', emit: map_log
    path '*.gtg-concordance.log', emit: concordance_log

    script:
    """
    set -euo pipefail
    test -s validation.success
    python3 "${gtg_runner}" \\
      --resolved-run "${resolved_run}" \\
      --normalized-ped "${normalized_ped}" \\
      --map-vcf "${map_vcf}" \\
      --all-sites-vcf "${all_sites_vcf}" \\
      --normalization-report "${normalization_report}" \\
      --output-dir .
    """
}


process GENERATE_REFERENCE_CPGS {
    tag 'configured autosomes'
    publishDir "${publish_root}/reference", mode: 'copy', overwrite: true

    input:
    tuple path(reference_fasta), path(reference_fai), val(regions), val(reference_name)
    path validation_success
    val publish_root
    path generator

    output:
    path "${reference_name}.autosomes.cpgs.bed.gz", emit: cpg_bed
    path "${reference_name}.autosomes.cpgs.bed.gz.tbi", emit: cpg_index
    path 'reference-cpgs-qc.json', emit: reference_qc

    script:
    """
    set -euo pipefail
    test -s validation.success
    python3 "${generator}" \\
      --fasta "${reference_fasta}" \\
      --fasta-index "${reference_fai}" \\
      --regions "${regions}" \\
      --reference-name "${reference_name}" \\
      --output-dir .
    """
}


process PHASE_MODEL_TO_FOUNDERS {
    tag "${sample_id}"
    publishDir "${publish_root}/samples/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id),
          path(combined_bed), path(combined_index),
          path(hap1_bed), path(hap1_index),
          path(hap2_bed), path(hap2_index),
          path(model_filter_qc),
          path(read_phased_vcf), path(read_phased_index),
          path(phase_blocks), path(reference_fai),
          val(reference_name), val(regions), val(min_coverage),
          val(mismatch_window_bp), val(bigwig_enabled),
          path(inheritance_vcf), path(inheritance_index), path(iht)
    path validation_success
    val publish_root
    path source_dir

    output:
    tuple val(sample_id),
          path("${sample_id}.dna-methylation.bed"),
          path("${sample_id}.dna-methylation.bed.header"),
          path("${sample_id}.hap-map-blocks.bed"),
          path("${sample_id}.hap-map-blocks.bed.header"),
          path("${sample_id}.bit-vector-sites-mismatches.bed"),
          path("${sample_id}.bit-vector-sites-mismatches.bed.header"),
          path("${sample_id}.bit-vector-sites-mismatches.vcf.gz"),
          path("${sample_id}.bit-vector-sites-mismatches.vcf.gz.tbi"),
          path("${sample_id}.phasing-qc.json"),
          emit: founder_results
    path "${sample_id}.dna-methylation.pat.model.${reference_name}.bw", optional: true, emit: paternal_bigwig
    path "${sample_id}.dna-methylation.mat.model.${reference_name}.bw", optional: true, emit: maternal_bigwig
    path "${sample_id}.hap-map-blocks.paternal.bed.gz"
    path "${sample_id}.hap-map-blocks.paternal.bed.gz.tbi"
    path "${sample_id}.hap-map-blocks.paternal.bed.header"
    path "${sample_id}.hap-map-blocks.maternal.bed.gz"
    path "${sample_id}.hap-map-blocks.maternal.bed.gz.tbi"
    path "${sample_id}.hap-map-blocks.maternal.bed.header"

    script:
    def bigwig_option = bigwig_enabled ? '' : '--no-bigwig'
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    test -s validation.success
    export PYTHONPATH="${source_dir}:${source_dir}/util"
    python3 "${source_dir}/phase_meth_to_founder_haps.py" \\
      --uid "${sample_id}" \\
      --vcf_read_phased "${read_phased_vcf}" \\
      --tsv_read_phase_blocks "${phase_blocks}" \\
      --vcf_iht_phased "${inheritance_vcf}" \\
      --txt_iht_blocks "${iht}" \\
      --bed_meth_model_hap1 "${hap1_bed}" \\
      --bed_meth_model_hap2 "${hap2_bed}" \\
      --reference_fai "${reference_fai}" \\
      --reference_name "${reference_name}" \\
      --output_dir . ${bigwig_option}
    """
}


process EXPAND_MODEL_TO_ALL_CPGS {
    tag "${sample_id}"
    label 'all_cpg'
    publishDir "${publish_root}/samples/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id),
          path(founder_bed), path(founder_header),
          path(hap_map_bed), path(hap_map_header),
          path(mismatch_bed), path(mismatch_header),
          path(mismatch_vcf), path(mismatch_index),
          path(phasing_qc),
          path(combined_model_bed), path(combined_model_index),
          val(regions), val(min_coverage), val(mismatch_window_bp),
          val(reference_name), val(config_fingerprint),
          path(reference_cpgs), path(reference_cpgs_index),
          path(joint_vcf), path(joint_vcf_index)
    path validation_success
    val publish_root
    path expander

    output:
    tuple val(sample_id),
          path("${sample_id}.dna-methylation.all-cpgs.bed.gz"),
          path("${sample_id}.dna-methylation.all-cpgs.bed.gz.tbi"),
          path("${sample_id}.all-cpgs-qc.json"),
          emit: all_cpg_results

    script:
    """
    set -euo pipefail
    test -s validation.success
    python3 "${expander}" \\
      --reference-cpgs "${reference_cpgs}" \\
      --combined-model-bed "${combined_model_bed}" \\
      --founder-bed "${founder_bed}" \\
      --founder-header "${founder_header}" \\
      --mismatch-bed "${mismatch_bed}" \\
      --mismatch-header "${mismatch_header}" \\
      --phasing-qc "${phasing_qc}" \\
      --joint-vcf "${joint_vcf}" \\
      --sample-id "${sample_id}" \\
      --regions "${regions}" \\
      --mismatch-window-bp "${mismatch_window_bp}" \\
      --reference-name "${reference_name}" \\
      --min-coverage "${min_coverage}" \\
      --config-fingerprint "${config_fingerprint}" \\
      --output-dir .
    """
}


process WRITE_RESULTS_MANIFEST {
    tag 'published outputs'
    publishDir "${publish_root}", mode: 'copy', overwrite: true

    input:
    path resolved_run
    path selected_samples
    path inheritance_qc
    path reference_qc
    path runtime_versions
    path sample_qc
    path validation_success
    val publish_root
    path manifest_writer

    output:
    stdout emit: completion_summary
    path 'results-manifest.json', emit: results_manifest

    script:
    def sample_qc_args = sample_qc.collect { qc -> "--sample-qc \"${qc}\"" }.join(' ')
    """
    set -euo pipefail
    test -s validation.success
    python3 "${manifest_writer}" \\
      --resolved-run "${resolved_run}" \\
      --selected-samples "${selected_samples}" \\
      ${sample_qc_args} \\
      --output results-manifest.json
    """
}


workflow RUN_VALIDATION {
    def settings = resolveRunSettings()
    VALIDATE_INPUTS(
        Channel.value(settings[0]),
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src/tapestry_validate.py", checkIfExists: true)),
        Channel.value(file("${projectDir}/schemas", checkIfExists: true)),
        Channel.value(settings[2])
    )
    VALIDATE_INPUTS.out.validation_summary.view { it.trim() }

    emit:
    resolved_run = VALIDATE_INPUTS.out.resolved_run
    resolved_manifest = VALIDATE_INPUTS.out.resolved_manifest
    normalized_ped = VALIDATE_INPUTS.out.normalized_ped
    selected_samples = VALIDATE_INPUTS.out.selected_samples
    selected_artifacts = VALIDATE_INPUTS.out.selected_artifacts
    validation_report = VALIDATE_INPUTS.out.validation_report
    config_fingerprint = VALIDATE_INPUTS.out.config_fingerprint
    validation_success = VALIDATE_INPUTS.out.validation_success
}


workflow validate {
    RUN_VALIDATION()
}


workflow {
    // Every scientific process takes validation.success as an explicit dependency.
    def settings = resolveRunSettings()
    RUN_VALIDATION()
    validation_gate = RUN_VALIDATION.out.validation_success
    selected_records = RUN_VALIDATION.out.selected_artifacts
        .splitJson(path: 'samples')
    model_artifacts = selected_records
        .map { sample ->
            tuple(
                sample.id,
                file(sample.cpg_model.combined.bed, checkIfExists: true),
                file(sample.cpg_model.combined.index, checkIfExists: true),
                file(sample.cpg_model.hap1.bed, checkIfExists: true),
                file(sample.cpg_model.hap1.index, checkIfExists: true),
                file(sample.cpg_model.hap2.bed, checkIfExists: true),
                file(sample.cpg_model.hap2.index, checkIfExists: true),
                sample.min_coverage,
                sample.regions.join(',')
            )
        }
    phase_artifacts = selected_records
        .map { sample ->
            tuple(
                sample.id,
                file(sample.phased_small_variants.vcf, checkIfExists: true),
                file(sample.phased_small_variants.index, checkIfExists: true),
                file(sample.phase_blocks, checkIfExists: true),
                file(sample.reference_fai, checkIfExists: true),
                sample.reference_name,
                sample.regions.join(','),
                sample.min_coverage,
                sample.mismatch_window_bp,
                sample.bigwig,
                sample.config_fingerprint
            )
        }
    reference_artifact = selected_records
        .map { sample ->
            tuple(
                file(sample.reference_fasta, checkIfExists: true),
                file(sample.reference_fai, checkIfExists: true),
                sample.regions.join(','),
                sample.reference_name
            )
        }
        .unique()
        .first()
    FILTER_MODEL_BEDS(
        model_artifacts,
        validation_gate,
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src/filter_model_beds.py", checkIfExists: true))
    )
    CAPTURE_RUNTIME_VERSIONS(
        RUN_VALIDATION.out.resolved_manifest,
        validation_gate,
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src/capture_versions.py", checkIfExists: true))
    )
    NORMALIZE_JOINT_VCF(
        RUN_VALIDATION.out.resolved_run,
        RUN_VALIDATION.out.resolved_manifest,
        RUN_VALIDATION.out.normalized_ped,
        validation_gate,
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src/normalize_joint_vcf.py", checkIfExists: true)),
        Channel.value(settings[2])
    )
    RUN_GTG_INHERITANCE(
        RUN_VALIDATION.out.resolved_run,
        RUN_VALIDATION.out.normalized_ped,
        NORMALIZE_JOINT_VCF.out.all_sites_vcf,
        NORMALIZE_JOINT_VCF.out.all_sites_index,
        NORMALIZE_JOINT_VCF.out.map_vcf,
        NORMALIZE_JOINT_VCF.out.map_index,
        NORMALIZE_JOINT_VCF.out.normalization_report,
        validation_gate,
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src/run_gtg_inheritance.py", checkIfExists: true))
    )

    GENERATE_REFERENCE_CPGS(
        reference_artifact,
        validation_gate,
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src/generate_reference_cpgs.py", checkIfExists: true))
    )

    filter_phase_inputs = FILTER_MODEL_BEDS.out.filtered_models.join(
        phase_artifacts, by: 0
    )
    inheritance_bundle = RUN_GTG_INHERITANCE.out.pass_vcf
        .combine(RUN_GTG_INHERITANCE.out.pass_index)
        .combine(RUN_GTG_INHERITANCE.out.iht)
    founder_metadata = filter_phase_inputs.map { row ->
        tuple(
            row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7],
            row[8], row[9], row[10], row[11], row[12], row[13], row[14],
            row[15], row[16]
        )
    }
    founder_inputs = founder_metadata.combine(inheritance_bundle)
    PHASE_MODEL_TO_FOUNDERS(
        founder_inputs,
        validation_gate,
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src", checkIfExists: true))
    )

    all_cpg_metadata = filter_phase_inputs.map { row ->
        tuple(
            row[0], row[1], row[2], row[13], row[14], row[15], row[12], row[17]
        )
    }
    founder_with_model = PHASE_MODEL_TO_FOUNDERS.out.founder_results.join(
        all_cpg_metadata, by: 0
    )
    reference_bundle = GENERATE_REFERENCE_CPGS.out.cpg_bed.combine(
        GENERATE_REFERENCE_CPGS.out.cpg_index
    )
    joint_bundle = NORMALIZE_JOINT_VCF.out.all_sites_vcf.combine(
        NORMALIZE_JOINT_VCF.out.all_sites_index
    )
    all_cpg_inputs = founder_with_model
        .combine(reference_bundle)
        .combine(joint_bundle)
    EXPAND_MODEL_TO_ALL_CPGS(
        all_cpg_inputs,
        validation_gate,
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src/expand_model_to_all_cpgs.py", checkIfExists: true))
    )

    completed_sample_qc = EXPAND_MODEL_TO_ALL_CPGS.out.all_cpg_results
        .map { sample_id, bed, index, qc -> qc }
        .collect()
    WRITE_RESULTS_MANIFEST(
        RUN_VALIDATION.out.resolved_run,
        RUN_VALIDATION.out.selected_samples,
        RUN_GTG_INHERITANCE.out.inheritance_qc,
        GENERATE_REFERENCE_CPGS.out.reference_qc,
        CAPTURE_RUNTIME_VERSIONS.out.versions,
        completed_sample_qc,
        validation_gate,
        Channel.value(settings[1]),
        Channel.value(file("${projectDir}/src/write_results_manifest.py", checkIfExists: true))
    )
    WRITE_RESULTS_MANIFEST.out.completion_summary.view { it.trim() }
}
