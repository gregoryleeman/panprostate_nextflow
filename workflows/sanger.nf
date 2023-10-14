nextflow.enable.dsl=2

include { getRefs; getNames } from '../modules/helpers'

process cnInit {
    errorStrategy 'finish'

    input:
    val GET_REFS_DONE

    output:
    val true, emit: done

    script:
    """
    mkdir -p ${params.TMP_DIR}
    awk -v FS='\t' -v OFS='\t' '{print \$1,0,\$2, 2}' ${params.REF_DIR}/genome.fa.fai > ${params.TMP_DIR}/norm.cn.bed 
    awk -v FS='\t' -v OFS='\t' '{print \$1,0,\$2, 2}' ${params.REF_DIR}/genome.fa.fai > ${params.TMP_DIR}/tum.cn.bed
    """
}

process cavemanSetup {
    errorStrategy 'finish'

    memory "${task.memory}"
    cpus "${task.cpus}"

    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val CN_INIT_DONE

    output:
    val true, emit: done

    script:
    """
    caveman.pl \
     -r ${params.REF_DIR}/genome.fa.fai \
     -ig ${params.REF_DIR}/caveman/HiDepth.tsv \
     -b ${params.REF_DIR}/caveman/flagging \
     -ab ${params.REF_DIR}/vagrent \
     -u ${params.REF_DIR}/caveman \
     -s ${params.SPECIES} \
     -sa ${params.ASSEMBLY} \
     -t ${task.cpus} \
     -st ${params.PROTOCOL} \
     -tc ${params.TMP_DIR}/tum.cn.bed \
     -nc ${params.TMP_DIR}/norm.cn.bed \
     -td 5 -nd 2 \
     -tb ${TUMOR} \
     -nb ${NORMAL} \
     -c ${params.REF_DIR}/caveman/flag.vcf.config.WGS.ini \
     -f ${params.REF_DIR}/caveman/flagging/flag.to.vcf.convert.ini \
     -e 350000 \
     -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
     -x ${params.CONTIG_EXCLUDE} \
     -p setup
     """
}

process compareBamGenotypes {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val GET_REFS_DONE

    output:
    val true, emit: done

    script:
    """
    compareBamGenotypes.pl \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/genotyped \
    -nb ${NORMAL} \
    -j ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/genotyped/result.json \
    -tb ${TUMOR} \
    -s ${params.REF_DIR}/general.tsv \
    -g ${params.REF_DIR}/gender.tsv
    """
}

process verifyBamHomChkNormal {
    errorStrategy 'finish'

    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val GET_REFS_DONE

    script:
    """
    verifyBamHomChk.pl -d 25 \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_WT}/contamination \
    -b ${NORMAL} \
    -t ${task.cpus} \
    -j ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_WT}/contamination/result.json \
    -s ${params.REF_DIR}/verifyBamID_snps.vcf.gz
    """
}

process cavemanSplit {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val CAVEMAN_SETUP_DONE

    output:
    val true, emit: done

    script:
    """
    echo "CWD=\${PWD}" 1<> ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/tmpCaveman/caveman.cfg.ini
    caveman.pl \
     -r ${params.REF_DIR}/genome.fa.fai \
     -ig ${params.REF_DIR}/caveman/HiDepth.tsv \
     -b ${params.REF_DIR}/caveman/flagging \
     -ab ${params.REF_DIR}/vagrent \
     -u ${params.REF_DIR}/caveman \
     -s ${params.SPECIES} \
     -sa ${params.ASSEMBLY} \
     -t ${task.cpus} \
     -st ${params.PROTOCOL} \
     -tc ${params.TMP_DIR}/tum.cn.bed \
     -nc ${params.TMP_DIR}/norm.cn.bed \
     -td 5 -nd 2 \
     -tb ${TUMOR} \
     -nb ${NORMAL} \
     -c ${params.REF_DIR}/caveman/flag.vcf.config.WGS.ini \
     -f ${params.REF_DIR}/caveman/flagging/flag.to.vcf.convert.ini \
     -e 350000 \
     -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
     -x ${params.CONTIG_EXCLUDE} \
     -p split
     """
}

process ascat {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val GET_REFS_DONE

    output:
    val "${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat/${NAME_MT}.samplestatistics.txt", emit: sample_stats
    val "${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat/${NAME_MT}.copynumber.caveman.csv", emit: CN
    val true, emit: done

    script:
    """
    ASCAT_ADD_ARGS=''

    if [ ! -z ${params.ASCAT_PURITY} ]; then
    ASCAT_ADD_ARGS="\$ASCAT_ADD_ARGS -pu ${params.ASCAT_PURITY}"
    fi
    if [ ! -z ${params.ASCAT_PLOIDY} ]; then
    ASCAT_ADD_ARGS="\$ASCAT_ADD_ARGS -pi ${params.ASCAT_PLOIDY}"
    fi

    ascat.pl \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat \
    -t ${TUMOR} \
    -n ${NORMAL} \
    -sg ${params.REF_DIR}/ascat/SnpGcCorrections.tsv \
    -r ${params.REF_DIR}/genome.fa \
    -q 20 \
    -g L \
    -l ${params.REF_DIR}/gender.tsv \
    -rs ${params.SPECIES} \
    -ra ${params.ASSEMBLY} \
    -pr ${params.PROTOCOL} \
    -pl ILLUMINA \
    -c ${task.cpus} \
    -force \
    \$ASCAT_ADD_ARGS
    """
}

process brassInput {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val GET_REFS_DONE

    output:
    val true, emit: done

    script:
    """
    brass.pl -j 4 -k 4 -c ${task.cpus} \
    -d ${params.REF_DIR}/brass/HiDepth.bed.gz \
    -f ${params.REF_DIR}/brass/brass_np.groups.gz \
    -g ${params.REF_DIR}/genome.fa \
    -s ${params.SPECIES} -as ${params.ASSEMBLY} -pr ${params.PROTOCOL} -pl ILLUMINA \
    -g_cache ${params.REF_DIR}/vagrent/vagrent.cache.gz \
    -vi ${params.REF_DIR}/brass/viral.genomic.fa.2bit \
    -mi ${params.REF_DIR}/brass/all_ncbi_bacteria \
    -b ${params.REF_DIR}/brass/500bp_windows.gc.bed.gz \
    -ct ${params.REF_DIR}/brass/CentTelo.tsv \
    -cb ${params.REF_DIR}/brass/cytoband.txt \
    -t ${TUMOR} \
    -n ${NORMAL} \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/brass \
    -p input
    """
}

process brassCover {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val BRASS_INPUT_DONE

    output:
    val true, emit: done

    script:
    """
    brass.pl -j 4 -k 4 -c ${task.cpus} \
    -d ${params.REF_DIR}/brass/HiDepth.bed.gz \
    -f ${params.REF_DIR}/brass/brass_np.groups.gz \
    -g ${params.REF_DIR}/genome.fa \
    -s ${params.SPECIES} -as ${params.ASSEMBLY} -pr ${params.PROTOCOL} -pl ILLUMINA \
    -g_cache ${params.REF_DIR}/vagrent/vagrent.cache.gz \
    -vi ${params.REF_DIR}/brass/viral.genomic.fa.2bit \
    -mi ${params.REF_DIR}/brass/all_ncbi_bacteria \
    -b ${params.REF_DIR}/brass/500bp_windows.gc.bed.gz \
    -ct ${params.REF_DIR}/brass/CentTelo.tsv \
    -cb ${params.REF_DIR}/brass/cytoband.txt \
    -t ${TUMOR} \
    -n ${NORMAL} \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/brass \
    -p cover
    """
}

process prepAscat {
    errorStrategy 'finish'

    input:
    tuple val(NAME_MT), val(NAME_WT)
    val ASCAT_CN
    val MT_SAMPLE_STATS

    output:
    tuple val("${params.TMP_DIR}/${NAME_MT}_vs_${NAME_WT}/norm.cn.bed"), val("${params.TMP_DIR}/${NAME_MT}_vs_${NAME_WT}/tum.cn.bed"), emit: CN_files
    env NORM_CONTAM, emit: contam
    val true, emit: done

    script:
    """
    mkdir -p ${params.TMP_DIR}/${NAME_MT}_vs_${NAME_WT}/
    perl -ne '@F=(split q{,}, \$_)[1,2,3,4]; \$F[1]-1; print join("\t",@F)."\n";' < $ASCAT_CN > ${params.TMP_DIR}/${NAME_MT}_vs_${NAME_WT}/norm.cn.bed
    perl -ne '@F=(split q{,}, \$_)[1,2,3,6]; \$F[1]-1; print join("\t",@F)."\n";' < $ASCAT_CN > ${params.TMP_DIR}/${NAME_MT}_vs_${NAME_WT}/tum.cn.bed
    NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){print 1-\$1;}' ${MT_SAMPLE_STATS}`
    """
}

process pindel {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val GET_REFS_DONE

    output:
    val "${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.germline.bed", emit: germline
    val true, emit: done

    script:
    """
    pindel.pl \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel \
    -r ${params.REF_DIR}/genome.fa \
    -t ${TUMOR} \
    -n ${NORMAL} \
    -s ${params.REF_DIR}/pindel/simpleRepeats.bed.gz \
    -u ${params.REF_DIR}/pindel/pindel_np.gff3.gz \
    -f ${params.REF_DIR}/pindel/${params.PROTOCOL}_Rules.lst \
    -g ${params.REF_DIR}/vagrent/codingexon_regions.indel.bed.gz \
    -st ${params.PROTOCOL} \
    -as ${params.ASSEMBLY} \
    -sp ${params.SPECIES} \
    -e ${params.CONTIG_EXCLUDE} \
    -b ${params.REF_DIR}/pindel/HiDepth.bed.gz \
    -c ${task.cpus} \
    -sf ${params.REF_DIR}/pindel/softRules.lst
    """
}

process cavemanContam {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    tuple val(NORMAL_BED), val(TUMOR_BED)
    val NORM_CONTAM
    val CAVEMAN_SPLIT_DONE

    output:
    val true, emit: done

    script:
    """
    echo "CWD=\${PWD}" 1<> ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/tmpCaveman/caveman.cfg.ini
    caveman.pl \
    -r ${params.REF_DIR}/genome.fa.fai \
    -ig ${params.REF_DIR}/caveman/HiDepth.tsv \
    -b ${params.REF_DIR}/caveman/flagging \
    -ab ${params.REF_DIR}/vagrent \
    -u ${params.REF_DIR}/caveman \
    -s ${params.SPECIES} \
    -sa ${params.ASSEMBLY} \
    -t ${task.cpus} \
    -st ${params.PROTOCOL} \
    -tc ${TUMOR_BED} \
    -nc ${NORMAL_BED} \
    -td 5 -nd 2 \
    -tb ${TUMOR} \
    -nb ${NORMAL} \
    -c ${params.SNVFLAG} \
    -f ${params.REF_DIR}/caveman/flagging/flag.to.vcf.convert.ini \
    -e ${params.CAVESPLIT} \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
    -x ${params.CONTIG_EXCLUDE} \
    -k ${NORM_CONTAM} \
    -no-flagging -noclean
    """
}

process indexPindelGermline {
    errorStrategy 'finish'
    
    input:
    val GERMLINE_BED 

    output:
    val "${GERMLINE_BED}.gz"
    val true, emit: done

    script:
    """
    if [ -f ${GERMLINE_BED} ]; then
    # need to sort and index pindel germline
    sort -k1,1 -k2,2n -k3,3n ${GERMLINE_BED} | bgzip -c > ${GERMLINE_BED}.gz
    tabix -p bed ${GERMLINE_BED}.gz
    rm -f ${GERMLINE_BED}
    fi
    """
}

process addBrass {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    tuple val(MT_SAMPLE_STATS), val(CN)
    val BRASS_COVER_DONE

    output:
    val true, emit: done

    script:
    """
    brass.pl -j 4 -k 4 -c ${task.cpus} \
    -d ${params.REF_DIR}/brass/HiDepth.bed.gz \
    -f ${params.REF_DIR}/brass/brass_np.groups.gz \
    -g ${params.REF_DIR}/genome.fa \
    -s ${params.SPECIES} -as ${params.ASSEMBLY} -pr ${params.PROTOCOL} -pl ILLUMINA \
    -g_cache ${params.REF_DIR}/vagrent/vagrent.cache.gz \
    -vi ${params.REF_DIR}/brass/viral.genomic.fa.2bit \
    -mi ${params.REF_DIR}/brass/all_ncbi_bacteria \
    -b ${params.REF_DIR}/brass/500bp_windows.gc.bed.gz \
    -ct ${params.REF_DIR}/brass/CentTelo.tsv \
    -cb ${params.REF_DIR}/brass/cytoband.txt \
    -t ${TUMOR} \
    -n ${NORMAL} \
    -ss ${MT_SAMPLE_STATS} \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/brass
    """
}

process pindelAnnotation {
    errorStrategy 'finish'

    input:
    tuple val(NAME_MT), val(NAME_WT)
    val PINDEL_DONE

    script:
    """
    rm -f ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.annot.vcf.gz*

    AnnotateVcf.pl -t -c ${params.REF_DIR}/vagrent/vagrent.cache.gz \
    -i ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.flagged.vcf.gz \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.annot.vcf
    """
}

process cavemanFlag {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    tuple val(NORMAL_BED), val(TUMOR_BED)
    val NORM_CONTAM
    val GERMLINE_BED
    val CAVEMAN_CONTAM_DONE

    output:
    val "${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.flagged.muts.vcf.gz", emit: vcf
    val true, emit: done

    script:
    """
    echo "CWD=\${PWD}" 1<> ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/tmpCaveman/caveman.cfg.ini
    caveman.pl \
    -r ${params.REF_DIR}/genome.fa.fai \
    -ig ${params.REF_DIR}/caveman/HiDepth.tsv \
    -b ${params.REF_DIR}/caveman/flagging \
    -ab ${params.REF_DIR}/vagrent \
    -u ${params.REF_DIR}/caveman \
    -s ${params.SPECIES} \
    -sa ${params.ASSEMBLY} \
    -t ${task.cpus} \
    -st ${params.PROTOCOL} \
    -tc ${TUMOR_BED} \
    -nc ${NORMAL_BED} \
    -td 5 -nd 2 \
    -tb ${TUMOR} \
    -nb ${NORMAL} \
    -c ${params.SNVFLAG} \
    -f ${params.REF_DIR}/caveman/flagging/flag.to.vcf.convert.ini \
    -e ${params.CAVESPLIT} \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
    -x ${params.CONTIG_EXCLUDE} \
    -k ${NORM_CONTAM} \
    -in ${GERMLINE_BED} \
    -p flag
    """
}

process cavemanAnnotation {
    errorStrategy 'finish'
    
    input:
    tuple val(NAME_MT), val(NAME_WT)
    val CAVEMAN_FLAG

    script:
    """
    rm -f ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.annot.muts.vcf.gz*

    AnnotateVcf.pl -t -c ${params.REF_DIR}/vagrent/vagrent.cache.gz \
    -i ${CAVEMAN_FLAG} \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.annot.muts.vcf
    """
}

process verifyBamHomChkTumor {
    errorStrategy 'finish'

    input:
    tuple val(TUMOR), val(NORMAL)
    tuple val(NAME_MT), val(NAME_WT)
    val ASCAT_CN

    script:
    """
    verifyBamHomChk.pl -d 25 \
    -o ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}/contamination \
    -b ${TUMOR} \
    -t ${task.cpus} \
    -a ${ASCAT_CN} \
    -j ${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}/contamination/result.json \
    -s ${params.REF_DIR}/verifyBamID_snps.vcf.gz
    """
}

workflow SANGER {

    sample_files = Channel
        .fromPath( params.INPUT_CSV )
        .splitCsv()
        .map { it -> tuple(it[0], it[1]) }

    sample_files.view()

    getRefs()

    getNames(
        sample_files
    )

    getNames.out.sample_names.view()

    cnInit(
        getRefs.out.done
    )

    cavemanSetup(
        sample_files,
        getNames.out.sample_names,
        cnInit.out.done
    )


    cavemanSplit(
        sample_files,
        getNames.out.sample_names,
        cavemanSetup.out.done
    )


    ascat(
        sample_files,
        getNames.out.sample_names,
        getRefs.out.done
    )

    compareBamGenotypes(
        sample_files,
        getNames.out.sample_names,
        getRefs.out.done
    )

    prepAscat(
        sample_files,
        ascat.out.CN,
        ascat.out.sample_stats
    )

    cavemanContam(
        sample_files,
        getNames.out.sample_names,
        prepAscat.out.CN_files,
        prepAscat.out.contam,
        cavemanSplit.out.done
    )

    brassInput(
        sample_files,
        getNames.out.sample_names,
        getRefs.out.done
    )

    brassCover(
        sample_files,
        getNames.out.sample_names,
        brassInput.out.done
    )

    addBrass(
        sample_files,
        getNames.out.sample_names,
        ascat.out.sample_stats,
        brassCover.out.done
    )

    pindel(
        sample_files,
        getNames.out.sample_names,
        getRefs.out.done
    )

    indexPindelGermline(
        pindel.out.germline
    )

    pindelAnnotation(
        getNames.out.sample_names,
        pindel.out.done
    )

    cavemanFlag(
        sample_files,
        getNames.out.sample_names,
        prepAscat.out.CN_files,
        prepAscat.out.contam,
        indexPindelGermline,
        cavemanContam.out.done
    )

    cavemanAnnotation(
        sample_files,
        cavemanFlag.out.vcf
    )

    verifyBamHomChkTumor(
        sample_files,
        getNames.out.sample_names,
        ascat.out.CN
    )

    verifyBamHomChkNormal(
        sample_files,
        getNames.out.sample_names,
        getRefs.out.done
    )

}
