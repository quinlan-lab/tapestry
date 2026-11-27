include { TRIO_PHASING } from './workflows/trio'
include { PEDIGREE } from './workflows/pedigree'

nextflow.enable.dsl = 2

params.outputsJson = null
params.ped = null
params.referenceFasta = null
params.referenceIndex = null
params.referenceGzi = null
params.outdir = null
params.projectId = null
params.samples = null
params.regions = null
params.mapMinQual = 20
params.mapMinDepth = 10
params.minRunMarkers = 10
params.concordanceMinQual = 20
params.concordanceMinDepth = 5
params.minCoverage = 10
params.mismatchWindowBp = 50
params.qcDiscordanceThreshold = 0.4
params.qcMinPairedCpgs = 100
params.bigwig = true
params.container = null
params.jointVcf = null
params.childBam = null
params.fatherBam = null
params.motherBam = null
params.chromosomes = null
params.region = null
params.crossoverDistanceBp = 10000
params.trioMinMapq = 1


workflow {
    TRIO_PHASING([
        ped: params.ped,
        jointVcf: params.jointVcf,
        reference: params.referenceFasta,
        childBam: params.childBam,
        fatherBam: params.fatherBam,
        motherBam: params.motherBam,
        chromosomes: params.chromosomes,
        region: params.region,
        crossoverDistanceBp: params.crossoverDistanceBp,
        minCoverage: params.minCoverage,
        minMapq: params.trioMinMapq,
        outdir: params.outdir
    ])
}


workflow pedigree {
    PEDIGREE([
        outputsJson: params.outputsJson,
        ped: params.ped,
        referenceFasta: params.referenceFasta,
        referenceIndex: params.referenceIndex,
        referenceGzi: params.referenceGzi,
        outdir: params.outdir,
        projectId: params.projectId,
        samples: params.samples,
        regions: params.regions,
        mapMinQual: params.mapMinQual,
        mapMinDepth: params.mapMinDepth,
        minRunMarkers: params.minRunMarkers,
        concordanceMinQual: params.concordanceMinQual,
        concordanceMinDepth: params.concordanceMinDepth,
        minCoverage: params.minCoverage,
        mismatchWindowBp: params.mismatchWindowBp,
        qcDiscordanceThreshold: params.qcDiscordanceThreshold,
        qcMinPairedCpgs: params.qcMinPairedCpgs,
        bigwig: params.bigwig
    ])
}
