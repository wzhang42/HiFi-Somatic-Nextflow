process CNVKit_Somatic {
    publishDir "${params.outdir}/${params.project}/CNVKit_Somatic/${SampleName}", mode: 'copy', overwrite: true

    input:
        tuple(
            val(SampleName),
            path(normal_bam),
            path(normal_bam_bai),
            path(tumor_bam),
            path(tumor_bam_bai)
        )
        val(Ref_Fa) // path(Ref_Fa)        

     output:
        tuple(
            val(SampleName),
            path("${SampleName}_cnvkit"),
         )

         tuple(
            val(SampleName),
            path("${SampleName}_cnvkit"),
            emit: To_CNVKit_Recall
         )

     script:
     """
      cnvkit.py batch \
        ${tumor_bam} \
        --normal ${normal_bam} \
        --annotate ${params.CNVKit.refFlat} \
        -f ${Ref_Fa} \
        --target-avg-size ${params.CNVKit.target_avg_size} \
        -m wgs \
        -p ${task.cpus} \
        --diagram --scatter \
        --output-dir ${SampleName}_cnvkit \
        --segment-method cbs
     """
}
