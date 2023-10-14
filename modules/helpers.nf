nextflow.enable.dsl=2

process getNames {
    errorStrategy 'finish'
    
    input:
    tuple val(TUMOR), val(NORMAL)

    output:
    tuple env(TUMOR_NAME), env(NORMAL_NAME), emit: sample_names
    val true, emit: done

    script:
    """
    TUMOR_NAME=`samtools view -H ${TUMOR} | grep "^@RG" | sed -e 's/^.*SM:\\(\\S*\\).*/\\1/g' | uniq`
    NORMAL_NAME=`samtools view -H ${NORMAL} | grep "^@RG" | sed -e 's/^.*SM:\\(\\S*\\).*/\\1/g' | uniq`
    """
}

process getRefs {

    output:
    val true, emit: done

    script:
    """
    echo "DONE";
    """
}

/*

process getRefs {

    // TODO gatk sequence dictionary

	output:
    val true, emit: done

    script:
    """
    if  [ ! -f "${params.REF_DIR}/genome.fa" ] || \
        [ ! -f "${params.REF_DIR}/genome.fa.fai" ] || \
        [ ! -f "${params.REF_DIR}/genome.fa.dict" ]; then 
        wget -qO- ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz | tar -C ${params.REF_DIR} -xzv --strip-components=1;
    fi;        

    if  [ ! -f "${params.REF_DIR}/gender.tsv" ] || \
        [ ! -f "${params.REF_DIR}/general.tsv" ] || \
        [ ! -f "${params.REF_DIR}/verifyBamID_snps.vcf.gz" ]; then
        wget -qO- ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/qcGenotype_GRCh38_hla_decoy_ebv.tar.gz | tar -C ${params.REF_DIR} -xzv --strip-components=1;
    fi;        

    if  [ ! -f "${params.REF_DIR}/ascat/SnpGcCorrections.tsv" ] || \
        [ ! -f "${params.REF_DIR}/brass/500bp_windows.gc.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/brass/500bp_windows.gc.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/brass/all_ncbi_bacteria.1.fa.2bit" ] || \
        [ ! -f "${params.REF_DIR}/brass/all_ncbi_bacteria.2.fa.2bit" ] || \
        [ ! -f "${params.REF_DIR}/brass/all_ncbi_bacteria.3.fa.2bit" ] || \
        [ ! -f "${params.REF_DIR}/brass/all_ncbi_bacteria.4.fa.2bit" ] || \
        [ ! -f "${params.REF_DIR}/brass/brass_np.groups.gz" ] || \
        [ ! -f "${params.REF_DIR}/brass/brass_np.groups.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/brass/CentTelo.tsv" ] || \
        [ ! -f "${params.REF_DIR}/brass/cytoband.txt" ] || \
        [ ! -f "${params.REF_DIR}/brass/HiDepth.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/brass/viral.genomic.fa.2bit" ]; then
        wget -qO- ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz | tar -C ${params.REF_DIR} -xzv --strip-components=1;
    fi;        

    if  [ ! -f "${params.REF_DIR}/caveman/flagging/centromeric_repeats.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/centromeric_repeats.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/flag.to.vcf.convert.ini" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/flag.vcf.config.ini" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/hi_seq_depth.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/hi_seq_depth.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/simple_repeats.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/simple_repeats.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/snps.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flagging/snps.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flag.vcf.config.WGS.ini" ] || \
        [ ! -f "${params.REF_DIR}/caveman/flag.vcf.config.WXS.ini" ] || \
        [ ! -f "${params.REF_DIR}/caveman/HiDepth.tsv" ] || \
        [ ! -f "${params.REF_DIR}/caveman/unmatchedNormal.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/caveman/unmatchedNormal.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/pindel/HiDepth.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/pindel/HiDepth.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/pindel/pindel_np.gff3.gz" ] || \
        [ ! -f "${params.REF_DIR}/pindel/pindel_np.gff3.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/pindel/simpleRepeats.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/pindel/simpleRepeats.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/pindel/softRules.lst" ] || \
        [ ! -f "${params.REF_DIR}/pindel/WGS_Rules.lst" ] || \
        [ ! -f "${params.REF_DIR}/pindel/WXS_Rules.lst" ]; then
        wget -qO- ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz | tar -C ${params.REF_DIR} -xzv --strip-components=1;
    fi;        

    if  [ ! -f "${params.REF_DIR}/vagrent/codingexon_regions.indel.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/codingexon_regions.indel.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/codingexon_regions.sub.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/codingexon_regions.sub.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/exon_regions.indel.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/exon_regions.indel.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/exon_regions.sub.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/exon_regions.sub.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/gene_regions.bed.gz" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/gene_regions.bed.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/vagrent.cache.gz" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/vagrent.cache.gz.tbi" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/vagrent.fa" ] || \
        [ ! -f "${params.REF_DIR}/vagrent/vagrent.fa.fai" ]; then
        wget -qO- ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz | tar -C ${params.REF_DIR} -xzv --strip-components=1;
    fi;        

    """

}
*/
