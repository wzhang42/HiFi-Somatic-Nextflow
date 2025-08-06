process DeepSomatic_SNV {
    publishDir "${params.outdir}/${params.project}/Deepsomatic_SNV/${SampleName}", mode: 'copy', overwrite: true

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
            path("${SampleName}_deepsomatic.vcf.gz"),
            path("${SampleName}_deepsomatic.vcf.gz.tbi"),
            path("${SampleName}_deepsomatic.g.vcf.gz"),
            path("${SampleName}_deepsomatic.g.vcf.gz.tbi"),
            path("${SampleName}_deepsomatic.visual_report.html")
        )
    
        tuple(
            val(SampleName),
            path("${SampleName}_deepsomatic.vcf.gz"),
            path("${SampleName}_deepsomatic.vcf.gz.tbi"),
            emit: To_MutationalPattern
         )

        tuple(
            val(SampleName),
            path("${SampleName}_deepsomatic.vcf.gz"),
            path("${SampleName}_deepsomatic.vcf.gz.tbi"),
            emit: To_PredictHRD
        )

        tuple(
            val(SampleName),
            path("${SampleName}_deepsomatic.vcf.gz"),
            path("${SampleName}_deepsomatic.vcf.gz.tbi"),
            emit: To_Hiphasing_Somatic_Vcf
        )
     
        tuple(
            val(SampleName),
            path("${SampleName}_deepsomatic.vcf.gz"),
            path("${SampleName}_deepsomatic.vcf.gz.tbi"),
            emit: To_Annot_Somatic
        )

    // --sample_name_tumor=${SampleName}.tumor \
    // --sample_name_normal=${SampleName}.normal \
    
    script:
    """
    export TMPDIR=${params.TMP_DIR}
    ${params.Deepsomatic.CMD_Prefix} \
    /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
    --model_type=PACBIO \
    --ref=${Ref_Fa} \
    --reads_normal=${normal_bam} \
    --reads_tumor=${tumor_bam} \
    --output_vcf=${SampleName}_deepsomatic.vcf.gz \
    --output_gvcf=${SampleName}_deepsomatic.g.vcf.gz \
    --sample_name_tumor=${SampleName}.Tumor \
    --sample_name_normal=${SampleName}.Normal \
    --num_shards=${task.cpus} \
    --logging_dir=${SampleName}_logs
    """
}

process Clair3_SNV {
   publishDir "${params.outdir}/${params.project}/Clair3_SNV/${SampleName}", mode: 'copy', overwrite: true
   
   input:
      tuple(
           val(SampleName),
           val(Normal_Tumor),
           path(aligned_bam),
           path(aligned_bam_bai)
      )
      val(Ref_Fa)  //path(Ref_Fa)
   
   output:
      tuple(
            val(SampleName),
            val(Normal_Tumor),
            path("${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz"),
            path("${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz.tbi")
       )
      
       tuple(
            val(SampleName),
            val(Normal_Tumor),
            path("${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz"),
            path("${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz.tbi"),
            emit: To_Pair_Normal_Tumor
       )

   /*    tuple(
            val(SampleName),
            val(Normal_Tumor),
            path("${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz"),
            path("${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz.tbi"),
            emit Clair3_SNV_Ch1
       )
   */

   // gunzip -c clair3_${SampleName}/merge_output.vcf.gz |awk '{if(\$7!="LowQual") {print \$0}}' OFS='\t' | bgzip -c > ${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz 
   // tabix -p vcf ${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz 

   script:
   """
        export TMPDIR=${params.TMP_DIR}
        /opt/bin/run_clair3.sh \
        --bam_fn=${aligned_bam} \
        --ref_fn=${Ref_Fa} \
        --threads=${task.cpus} \
        --platform=${params.Clair3_SNV.clair_platform} \
        --model_path="/opt/models/${params.Clair3_SNV.clair_model}" \
        --output=clair3_${SampleName}_${Normal_Tumor} \
        --sample_name=${SampleName}.${Normal_Tumor}
       
        gunzip -c clair3_${SampleName}_${Normal_Tumor}/merge_output.vcf.gz |awk 'BEGIN{OFS="\t"};{if(\$7!="LowQual") print \$0}' | bgzip -c > ${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz
        tabix -p vcf ${SampleName}.${Normal_Tumor}.clair3.small_variants.vcf.gz
    """
}

process Pair_Normal_Tumor_Vcfs {
  input:
       tuple(
            val(SampleName),
            val(Normal_Tumors),
            path(Vcf_gzs),
            path(Vcf_gz_tbis)
        )

  output:
       tuple(
           val(SampleName),
           path("${SampleName}_normal.vcf.gz"),
           path("${SampleName}_normal.vcf.gz.tbi"),
           path("${SampleName}_tumor.vcf.gz"),
           path("${SampleName}_tumor.vcf.gz.tbi")
       )

       tuple(
           val(SampleName),
           path("${SampleName}_normal.vcf.gz"),
           path("${SampleName}_normal.vcf.gz.tbi"),
      	   path("${SampleName}_tumor.vcf.gz"),
           path("${SampleName}_tumor.vcf.gz.tbi"),
           emit: To_CNVKit_Merge_Germline_VCFs
       )

       tuple(
           val(SampleName),
           path("${SampleName}_normal.vcf.gz"),
           path("${SampleName}_normal.vcf.gz.tbi"),
           emit: To_HiPhasing_Normal_Vcf
       )
       
       tuple(
           val(SampleName),
           path("${SampleName}_tumor.vcf.gz"),
           path("${SampleName}_tumor.vcf.gz.tbi"),
           emit: To_HiPhasing_Tumor_Vcf
       )
       
   script:
   if ( Normal_Tumors[0]== 'Tumor')
   """
      ln -s ${Vcf_gzs[1]} ${SampleName}_normal.vcf.gz
      ln -s ${Vcf_gz_tbis[1]} ${SampleName}_normal.vcf.gz.tbi
      ln -s ${Vcf_gzs[0]} ${SampleName}_tumor.vcf.gz
      ln -s ${Vcf_gz_tbis[0]} ${SampleName}_tumor.vcf.gz.tbi
   """
   else
   """
      ln -s ${Vcf_gzs[0]} ${SampleName}_normal.vcf.gz
      ln -s ${Vcf_gz_tbis[0]} ${SampleName}_normal.vcf.gz.tbi
      ln -s ${Vcf_gzs[1]} ${SampleName}_tumor.vcf.gz
      ln -s ${Vcf_gz_tbis[1]} ${SampleName}_tumor.vcf.gz.tbi
   """
} 

process ClairS_SNV {
   publishDir "${params.outdir}/${params.project}/ClairS_SNV/${SampleName}", mode: 'copy', overwrite: true
   
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
            path("${SampleName}_snv.vcf.gz"),
            path("${SampleName}_indel.vcf.gz"),
            path("${SampleName}_snv.vcf.gz.tbi"),
            path("${SampleName}_indel.vcf.gz.tbi"),
            path("${SampleName}_clair3_tumor_germline_output.vcf.gz"),
            path("${SampleName}_clair3_normal_germline_output.vcf.gz")
        )

        tuple(
            val(SampleName),
            path("${SampleName}_snv.vcf.gz"),
            path("${SampleName}_indel.vcf.gz"),
            path("${SampleName}_snv.vcf.gz.tbi"),
            path("${SampleName}_indel.vcf.gz.tbi"),
            path("${SampleName}_clair3_tumor_germline_output.vcf.gz"),
            path("${SampleName}_clair3_normal_germline_output.vcf.gz"),
            emit: ClairS_Ch1
         )

    script:
    """
      /opt/bin/run_clairs \
      --normal_bam_fn ${normal_bam} \
      --tumor_bam_fn ${tumor_bam} \
      --ref_fn ${Ref_Fa} \
      --threads ${task.cpus} \
      --platform {params.ClairS_SNV.platform} \
      --output_dir . \
      --conda_prefix /opt/conda/envs/clairs \
      --enable_indel_calling \
      --bed_fn ${params.ClairS_SNV.contig} \
      --enable_clair3_germline_output

      mv snv.vcf.gz ${SampleName}_snv.vcf.gz
      mv indel.vcf.gz ${SampleName}_indel.vcf.gz
      mv snv.vcf.gz.tbi ${SampleName}_snv.vcf.gz.tbi
      mv indel.vcf.gz.tbi ${SampleName}_indel.vcf.gz.tbi
      mv clair3_tumor_germline_output.vcf.gz ${SampleName}_clair3_tumor_germline_output.vcf.gz
      mv clair3_normal_germline_output.vcf.gz ${SampleName}_clair3_normal_germline_output.vcf.gz
    """
    
}
