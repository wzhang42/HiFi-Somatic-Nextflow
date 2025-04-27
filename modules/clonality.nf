process AMBER {
    publishDir "${params.outdir}/${params.project}/AMBER/${SampleName}", mode: 'copy', overwrite: true, pattern: "*.*" 

    input:
      tuple(
            val(SampleName),
            path(normal_bam),
            path(normal_bam_bai),
            path(tumor_bam),
            path(tumor_bam_bai)
      )
      val(Ref_Fa)  //path(Ref_Fa)
    
    output:
      tuple(
            val(SampleName),
            path("${SampleName}.Normal.amber.homozygousregion.tsv"),
            path("${SampleName}.Tumor.amber.baf.pcf"),
            path("${SampleName}.Tumor.amber.baf.tsv.gz"),
            path("${SampleName}.Tumor.amber.contamination.tsv"),
            path("${SampleName}.Tumor.amber.contamination.vcf.gz"),
            path("${SampleName}.Tumor.amber.contamination.vcf.gz.tbi"),
            path("${SampleName}.Tumor.amber.qc")
       )

      tuple(
            val(SampleName),
            path(AMBER_OutDir),
            emit: To_Purple
      )
    
    //-ref_genome_version GRCh38   mkdir AMBER_OutDir
    // -loci ensembl_data_dir/copy_number/GermlineHetPon.*.vcf.gz
    // -ref_genome_version  ${params.AMBER.REF_GENOME_VER} \
    // mkdir AMBER_OutDir cp ./*.* AMBER_OutDir
    // echo foo task path: \$PWD
    // export TMPDIR=${params.TMP_DIR}
    // export R_LIBS=${params.R_4_2_LIBS}
    // export _JAVA_OPTIONS=-Djava.io.tmpdir=${params.TMP_DIR}
    script:
    """
      export TMPDIR=${params.TMP_DIR}
      export _JAVA_OPTIONS=-Djava.io.tmpdir=${params.TMP_DIR}

      java -Xmx24g \
      -jar ${params.AMBER.BIN_JAR} \
      -reference ${SampleName}.Normal \
      -reference_bam ${normal_bam} \
      -tumor ${SampleName}.Tumor \
      -tumor_bam ${tumor_bam} \
      -output_dir ./ \
      -threads ${task.cpus} \
      -ref_genome ${Ref_Fa} \
      -ref_genome_version V38 \
      -loci ${params.ENSEMBLE_DATA_DIR}/copy_number/GermlineHetPon.*.vcf.gz
      
      ln -s "\${PWD}" AMBER_OutDir
    """
}

process COBALT {
    publishDir "${params.outdir}/${params.project}/COBALT/${SampleName}", mode: 'copy', overwrite: true, pattern: "*.*"

    input:
      tuple(
            val(SampleName),
            path(normal_bam),
            path(normal_bam_bai),
            path(tumor_bam),
            path(tumor_bam_bai)
      )
      val(Ref_Fa)  //path(Ref_Fa)
    
    output:
      tuple(
            val(SampleName),
            path("${SampleName}.Normal.cobalt.gc.median.tsv"),
            path("${SampleName}.Normal.cobalt.ratio.median.tsv"),
            path("${SampleName}.Normal.cobalt.ratio.pcf"),
            path("${SampleName}.Tumor.cobalt.gc.median.tsv"),
            path("${SampleName}.Tumor.cobalt.ratio.pcf"),
            path("${SampleName}.Tumor.cobalt.ratio.tsv.gz")
       )

       tuple(
            val(SampleName),
            path(COBALT_OutDir),
            emit: To_Purple
       )
    
     // -gc_profile ensembl_data_dir/copy_number/GC_profile.*.cnp
     // export TMPDIR=${params.TMP_DIR}
     // export _JAVA_OPTIONS=-Djava.io.tmpdir=${params.TMP_DIR}
     // export R_LIBS=${params.R_4_2_LIBS}
    script:
    """
       export TMPDIR=${params.TMP_DIR}
       export _JAVA_OPTIONS=-Djava.io.tmpdir=${params.TMP_DIR}

       java -Xmx24g \
       -jar ${params.COBALT.BIN_JAR} \
       -reference ${SampleName}.Normal \
       -reference_bam ${normal_bam} \
       -tumor ${SampleName}.Tumor \
       -tumor_bam ${tumor_bam} \
       -ref_genome ${Ref_Fa} \
       -output_dir ./ \
       -threads ${task.cpus} \
       -pcf_gamma 1000 \
       -gc_profile ${params.ENSEMBLE_DATA_DIR}/copy_number/GC_profile.*.cnp
       
       ln -s "\${PWD}" COBALT_OutDir
    """
}

process PURPLE {
    publishDir "${params.outdir}/${params.project}/PURPLE/${SampleName}", mode: 'copy', overwrite: true

    input:
      tuple(
            val(SampleName),
            path(AMBER_OutDir),
            path(COBALT_OutDir)
      )
      val(Ref_Fa)  //path(Ref_Fa)
    
    output:
      tuple(
            val(SampleName),
            path("${SampleName}.Tumor.purple.cnv.gene.tsv"),
            path("${SampleName}.Tumor.purple.cnv.somatic.tsv"),
            path("${SampleName}.Tumor.purple.germline.deletion.tsv"),
            path("${SampleName}.Tumor.purple.purity.range.tsv"),
            path("${SampleName}.Tumor.purple.purity.tsv"),
            path("${SampleName}.Tumor.purple.qc"),
            path("${SampleName}.Tumor.purple.segment.tsv"),
            path("${SampleName}.Tumor.purple.somatic.clonality.tsv"),
            path(PLOTS),
            path(CIRCOS)          
       )
    
    // -ensembl_data_dir ensembl_data_dir/common/ensembl_data
    // export TMPDIR=${params.TMP_DIR}
    // export _JAVA_OPTIONS=-Djava.io.tmpdir=${params.TMP_DIR}
    // export R_LIBS=${params.R_4_2_LIBS}
    script:
    """
       export TMPDIR=${params.TMP_DIR}
       export _JAVA_OPTIONS=-Djava.io.tmpdir=${params.TMP_DIR}

       java -Xmx24g \
       -jar ${params.PURPLE.BIN_JAR} \
       -reference ${SampleName}.Normal \
       -tumor ${SampleName}.Tumor \
       -output_dir . \
       -amber ${AMBER_OutDir} \
       -cobalt ${COBALT_OutDir} \
       -gc_profile ${params.ENSEMBLE_DATA_DIR}/copy_number/GC_profile.*.cnp \
       -ref_genome ${Ref_Fa} \
       -ref_genome_version ${params.PURPLE.REF_GENOME_VER} \
       -ensembl_data_dir ${params.ENSEMBLE_DATA_DIR}/common/ensembl_data \
       -highly_diploid_percentage ${params.PURPLE.highlyDiploidPercentage} \
       -somatic_min_purity_spread ${params.PURPLE.somaticMinPuritySpread} \
       -threads ${task.cpus} \
       -circos ${params.PURPLE.CIRCOS_BIN} \
       -max_purity ${params.PURPLE.max_purity} \
       -min_purity ${params.PURPLE.min_purity} \
       -min_ploidy ${params.PURPLE.min_ploidy} \
       -max_ploidy ${params.PURPLE.max_ploidy}
       
       ln -s "\${PWD}/plot" PLOTS
       ln -s "\${PWD}/circos" CIRCOS
    """
}
