 process DSS_DMR {
 publishDir "${params.outdir}/${params.project}/DSS_DMR/${SampleName}", mode: 'copy', overwrite: true
   input:
   tuple(
    val(SampleName),
    path(tumor_cpg_combined_bed),
    path(normal_cpg_combined_bed)
  )

  output:
  tuple (
    val(SampleName),
    path("${SampleName}.DMR.tsv") 
  )

  output:
  tuple (
    val(SampleName),
    path("${SampleName}.DMR.tsv"),
    emit: To_DMR_Annotation 
  )

  script:
  """
   gunzip -c ${tumor_cpg_combined_bed} | cut -f1,2,6,7 > ${SampleName}_tumor.tmp
   gunzip -c ${normal_cpg_combined_bed} | cut -f1,2,6,7 > ${SampleName}_normal.tmp

   Rscript --vanilla \
   /app/DSS_tumor_normal.R \
   ${SampleName}_tumor.tmp \
   ${SampleName}_normal.tmp \
   ${SampleName}.DMR.tsv \
   ${task.cpus}

   rm -f ${SampleName}_tumor.tmp ${SampleName}_normal.tmp
  """

}

process DMR_annotation {
   publishDir "${params.outdir}/${params.project}/DMR_annotation/${SampleName}", mode: 'copy', overwrite: true

   input:
   tuple (
    val(SampleName),
    path(DMR_tsv),
   )

   output:
   tuple (
    val(SampleName),
    path("${SampleName}*.tsv.gz") 
   )

   script:
   """
     Rscript --vanilla /app/annotatr_dmr.R \
     ${DMR_tsv} \
     ${SampleName} \
     ${task.cpus}
   """
}
