nextflow.enable.dsl=2

include { getRefs; getNames } from '../modules/helpers'

process gridss {
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(TUMOR_NAME), val(NORMAL_NAME)
    val GET_REFS_DONE

    output:
    path("*.vcf.gz") , emit: vcf, optional:true
    path("*.assembly.bam") , emit: assembly, optional:true
    val true, emit: done

    script:
    def AVAIL_MEM = 3
    if (!task.memory) {
        log.info '[GRIDSS] Available memory not known - defaulting to 3GB'
    } else {
        avail_mem = (task.memory.toGiga() - 1)
    }
    """
    ln -s ${params.REF_DIR}/genome.fa .;\\
    ln -s ${params.REF_DIR}/genome.fa.fai .;\\
    gridss \\
        --output gridss_${TUMOR_NAME}_${NORMAL_NAME}.vcf.gz \\
        --reference genome.fa \\
        --threads ${task.cpus} \\
        --jvmheap ${AVAIL_MEM}g \\
        --otherjvmheap ${AVAIL_MEM}g \\
        ${TUMOR} ${NORMAL};
    """
}

workflow GRIDSS {

    sample_files = Channel
        .fromPath( params.INPUT_CSV )
        .splitCsv()
        .map { it -> tuple(it[0], it[1]) }

    getRefs()

    getNames(
        sample_files
    )

    gridss(
        sample_files,
        getNames.out.sample_names,
        getRefs.out.done
    )

}
