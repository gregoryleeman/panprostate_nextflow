nextflow.enable.dsl=2

include { getRefs; getNames } from '../modules/helpers'

process mutect2 {

    container 'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(TUMOR_NAME), val(NORMAL_NAME)
    val GET_REFS_DONE

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi"), emit: tbi
    tuple val(meta), path("*.stats"), emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    val true, emit: done

    script:
    def AVAIL_MEM = 3072
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB.'
    } else {
        AVAIL_MEM = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${AVAIL_MEM}M -XX:-UsePerfData" \\
        Mutect2 \\
        --input ${TUMOR} --input ${NORMAL} \\
        --output mutect_${TUMOR_NAME}_${NORMAL_NAME}.vcf.gz \\
        --reference ${params.REF_DIR}/genome.fa \\
        --panel-of-normals ${params.REF_DIR}/TODO \\ # TODO vcf file to be used as a panel of normals. (*vcf.gz and *vcf.gz.tbi)
        --germline-resource ${params.REF_DIR}/TODO # TODO Population vcf of germline sequencing, containing allele fractions. (*vcf.gz and *vcf.gz.tbi)
        --intervals ${params.REF_DIR}/TODO \\ # TODO Specify region the tools is run on. (.bed)
        --tmp-dir . \\
    """
}

workflow MUTECT2 {

    sample_files = Channel
        .fromPath( params.INPUT_CSV )
        .splitCsv()
        .map { it -> tuple(it[0], it[1]) }

    getRefs()

    getNames(
        sample_files
    )

    mutect2(
        sample_files,
        getNames.out.sample_names,
        getRefs.out.done
    )

}
