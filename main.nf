nextflow.enable.dsl = 2

include { SANGER } from './workflows/sanger.nf'
include { MUTECT2 } from './workflows/mutect2.nf'
include { GRIDSS } from './workflows/gridss.nf'

params.each { key, value -> if (value == null) { exit 1, "Parameter $key not specified!" } }

workflow {

    file("${params.TMP_DIR}").mkdir()
    file("${params.REF_DIR}").mkdir()
    file("${params.OUTPUT_DIR}").mkdir()

    /* SANGER() */
    /* MUTECT2() */
    GRIDSS()

}
