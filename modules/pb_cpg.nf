process pb_cpg_score {
      publishDir "${params.outdir}/${params.project}/pb_cpg_score/$SampleName", mode: 'copy', overwrite: true 

      input:
      tuple(
            val(SampleName),
            val(Normal_Tumor),
            path(AlignedBAM),
            path(AlignedBAM_bai)
      )
      val(Ref_Fa)
        
      output:
      tuple(
           val(SampleName),
           path("${SampleName}_${Normal_Tumor}.combined.bed.gz"),
           path("${SampleName}_${Normal_Tumor}.hap1.bed.gz"),
           path("${SampleName}_${Normal_Tumor}.hap2.bed.gz")
       )
       
       tuple(
           val(SampleName),
           path("${SampleName}_${Normal_Tumor}.combined.bed.gz"),
           emit: To_DMR
       )
      
      //--model "$PILEUP_MODEL_DIR"/pileup_calling_model.v1.tflite
      script: 
      """
        export TMPDIR=${params.TMP_DIR} 
      
        aligned_bam_to_cpg_scores \
        --threads ${task.cpus} \
        --bam ${AlignedBAM} \
        --ref ${Ref_Fa} \
        --output-prefix ${SampleName}_${Normal_Tumor} \
        --min-mapq ${params.pb_cpg_score.minmapq} \
        --min-coverage ${params.pb_cpg_score.mincov} \
        --model ${params.pb_cpg_score.model}

        gzip ${SampleName}_${Normal_Tumor}.combined.bed
        gzip ${SampleName}_${Normal_Tumor}.hap1.bed
        gzip ${SampleName}_${Normal_Tumor}.hap2.bed
      """
}
