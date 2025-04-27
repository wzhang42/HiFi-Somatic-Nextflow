process Group_UnAligned_BAMs {
  input:
       tuple(
            val(SampleName),
            val(Normal_Tumor),
            path(HiFi_Unaligned_BAM)
        )

  output:
       tuple(
           val(SampleName),
           val(Normal_Tumor),
           env(LabelGroup_Tag),
           path(HiFi_Unaligned_BAM),
           emit: To_GroupBAM
       )

   script:
   """
     LabelGroup_Tag=${SampleName}_${Normal_Tumor}
   """
}

process PBMM2_Align {
    publishDir "${params.outdir}/${params.project}/PBMM2_Align/${SampleNames.first()}", mode: 'copy', overwrite: true
    input:
        tuple(
            val(SampleNames),
            val(Normal_Tumors),
            val(MergeGroup_Tag),
            path(HiFi_Unaligned_BAMs)
        )
        val(Ref_Fa) // path(Ref_Fa)
    
    output:
        tuple(
            env(SampleName),
            env(Normal_Tumor),
            path("${SampleNames.first()}.${Normal_Tumors.first()}.aln.bam"),
            path("${SampleNames.first()}.${Normal_Tumors.first()}.aln.bam.bai"),
            emit: To_Pair_Normal_Tumor
        )
        
        tuple(
            env(SampleName),
            env(Normal_Tumor),
            path("${SampleNames.first()}.${Normal_Tumors.first()}.aln.bam"),
            path("${SampleNames.first()}.${Normal_Tumors.first()}.aln.bam.bai"),
            emit: To_BAM_QC_mosdepth
        )

        tuple(
            env(SampleName),
            env(Normal_Tumor),
            path("${SampleNames.first()}.${Normal_Tumors.first()}.aln.bam"),
            path("${SampleNames.first()}.${Normal_Tumors.first()}.aln.bam.bai"),
            emit: To_BAM_QC_bam_coverage
       	)

        tuple(
            env(SampleName),
	    env(Normal_Tumor),
            path("${SampleNames.first()}.${Normal_Tumors.first()}.aln.bam"),
            path("${SampleNames.first()}.${Normal_Tumors.first()}.aln.bam.bai"),
            emit: To_Clair3_SNV
        )

    script:
    def SAMPLENAME=SampleNames.first()
    def NORMAL_TUMOR=Normal_Tumors.first()
    """
    SampleName=${SAMPLENAME}
    Normal_Tumor=${NORMAL_TUMOR}
    export TMPDIR=${params.TMP_DIR}

    echo ${HiFi_Unaligned_BAMs.collect{"$it "}.join()} | tr ' ' '\n'> ${SAMPLENAME}.${NORMAL_TUMOR}.fofn
    pbmm2 align \
    ${Ref_Fa} \
    ${SAMPLENAME}.${NORMAL_TUMOR}.fofn \
    ${SAMPLENAME}.${NORMAL_TUMOR}.aln.bam \
    --sample ${SAMPLENAME} \
    --sort \
    -j ${task.cpus} \
    -J 4 \
    --unmapped \
    --preset HIFI \
    --log-level INFO --log-file pbmm2.log    
    """
}

process Pair_Normal_Tumor_Bams {
  input: 
       tuple(
            val(SampleName),
            val(Normal_Tumors),
            path(PBMM2_AlignedBAMs),
            path(PBMM2_AlignedBAM_bais)
        )

  output:
       tuple(
           val(SampleName),
           path("${SampleName}_normal_aln.bam"),
           path("${SampleName}_normal_aln.bam.bai"),
           path("${SampleName}_tumor_aln.bam"),
           path("${SampleName}_tumor_aln.bam.bai"),           
           emit: To_DeepSomatic_SNV
       )

        tuple(
           val(SampleName),
           path("${SampleName}_normal_aln.bam"),
           path("${SampleName}_normal_aln.bam.bai"),
           path("${SampleName}_tumor_aln.bam"),
           path("${SampleName}_tumor_aln.bam.bai"),
           emit: To_Severus_SV
       )
       
       tuple(
           val(SampleName),
           path("${SampleName}_normal_aln.bam"),
           path("${SampleName}_normal_aln.bam.bai"),
           path("${SampleName}_tumor_aln.bam"),
           path("${SampleName}_tumor_aln.bam.bai"),
           emit: To_AMBER
       )
       
       tuple(
           val(SampleName),
           path("${SampleName}_normal_aln.bam"),
           path("${SampleName}_normal_aln.bam.bai"),
           path("${SampleName}_tumor_aln.bam"),
           path("${SampleName}_tumor_aln.bam.bai"),
           emit: To_COBALT
       )
       
       tuple(
           val(SampleName),
           path("${SampleName}_normal_aln.bam"),
           path("${SampleName}_normal_aln.bam.bai"),
           emit: To_HiPhasing_Normal_Bam
       )

        tuple(
           val(SampleName),
           path("${SampleName}_tumor_aln.bam"),
           path("${SampleName}_tumor_aln.bam.bai"),
           emit: To_HiPhasing_Tumor_Bam
       )
       
       tuple(
           val(SampleName),
           path("${SampleName}_normal_aln.bam"),
           path("${SampleName}_normal_aln.bam.bai"),
           path("${SampleName}_tumor_aln.bam"),
           path("${SampleName}_tumor_aln.bam.bai"),
           emit: To_CNVKit_Somatic
       )

   script:
   if ( Normal_Tumors[0]== 'Tumor')
   """
      ln -s ${PBMM2_AlignedBAMs[1]} ${SampleName}_normal_aln.bam 
      ln -s ${PBMM2_AlignedBAM_bais[1]} ${SampleName}_normal_aln.bam.bai 
      ln -s ${PBMM2_AlignedBAMs[0]} ${SampleName}_tumor_aln.bam 
      ln -s ${PBMM2_AlignedBAM_bais[0]} ${SampleName}_tumor_aln.bam.bai 
   """
   else  
   """
      ln -s ${PBMM2_AlignedBAMs[0]} ${SampleName}_normal_aln.bam
      ln -s ${PBMM2_AlignedBAM_bais[0]} ${SampleName}_normal_aln.bam.bai
      ln -s ${PBMM2_AlignedBAMs[1]} ${SampleName}_tumor_aln.bam
      ln -s ${PBMM2_AlignedBAM_bais[1]} ${SampleName}_tumor_aln.bam.bai
   """
}
