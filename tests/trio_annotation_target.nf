nextflow.enable.dsl = 2

include { annotationTarget } from '../workflows/trio'

workflow {
    // C at 200 is in the target; a variant at its G at 201 must be annotated.
    assert annotationTarget('chr22:100-200') == 'chr22:100-201'
    assert annotationTarget('chr22:200-200') == 'chr22:200-201'
    assert annotationTarget('chr22') == 'chr22'
    println 'PASS: CpG annotation includes the following G and preserves chromosome targets'
}
