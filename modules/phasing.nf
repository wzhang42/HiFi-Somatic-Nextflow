
process HiFiPhase_Somatic {
publishDir "${params.outdir}/${params.project}/HiFiPhase_Somatic/${SampleName}", mode: 'copy', overwrite: true
  input:
  tuple(
    val(SampleName),
    path(tumor_aligned_bam),
    path(tumor_aligned_bam_bai),
    path(tumor_germline_vcf),
    path(tumor_germline_vcf_tbi),
    path(Somatic_SNV_vcf),
    path(Somatic_SNV_vcf_tbi)
  )
  val(Ref_Fa)

  output:
  tuple(
    val(SampleName),
    path("${SampleName}.hiphase.bam"),
    path("${SampleName}.hiphase.bam.bai"),
    path("${SampleName}_tumor_germline.hiphase.vcf.gz"),
    path("${SampleName}_tumor_germline.hiphase.vcf.gz.tbi"),
    path("${SampleName}_somatic_snv.hiphase.vcf.gz"),
    path("${SampleName}_somatic_snv.hiphase.vcf.gz.tbi"),
    path("${SampleName}.hiphase.stats"),
    path("${SampleName}.hiphase.summary.tsv")
  )

  tuple(
    val(SampleName),
    val("Tumor"),
    path("${SampleName}.hiphase.bam"),
    path("${SampleName}.hiphase.bam.bai"),
    emit: To_Tumor_pb_cpg
  )

  script:
  """
    hiphase \
      --bam ${tumor_aligned_bam} \
      -t ${task.cpus} \
      --output-bam ${SampleName}.hiphase.bam \
      --vcf ${tumor_germline_vcf} \
      --output-vcf ${SampleName}_tumor_germline.hiphase.vcf.gz \
      --vcf ${Somatic_SNV_vcf} \
      --output-vcf ${SampleName}_somatic_snv.hiphase.vcf.gz \
      -r ${Ref_Fa} \
      --stats-file ${SampleName}.hiphase.stats \
      --summary-file ${SampleName}.hiphase.summary.tsv \
      --ignore-read-groups
  """
}
 
process HiFiPhase {
publishDir "${params.outdir}/${params.project}/HiFiPhase/${SampleName}", mode: 'copy', overwrite: true
  input:
  tuple(
    val(SampleName),
    path(normal_aligned_bam),
    path(normal_aligned_bam_bai),
    path(normal_germline_vcf),
    path(normal_germline_vcf_tbi)
  )
  val(Ref_Fa)

  output:
  tuple(
    val(SampleName),
    path("${SampleName}.hiphase.bam"),
    path("${SampleName}.hiphase.bam.bai"),
    path("${SampleName}_germline.hiphase.vcf.gz"),
    path("${SampleName}.hiphase.stats"),
    path("${SampleName}.hiphase.summary.tsv")
  )

  tuple(
    val(SampleName),
    val("Normal"),
    path("${SampleName}.hiphase.bam"),
    path("${SampleName}.hiphase.bam.bai"),
    emit: To_Normal_pb_cpg
  )

  script:
  """
    hiphase \
      --bam ${normal_aligned_bam} \
      -t ${task.cpus} \
      --output-bam ${SampleName}.hiphase.bam \
      --vcf ${normal_germline_vcf} \
      --output-vcf ${SampleName}_germline.hiphase.vcf.gz \
      -r ${Ref_Fa} \
      --stats-file ${SampleName}.hiphase.stats \
      --summary-file ${SampleName}.hiphase.summary.tsv \
      --ignore-read-groups
  """
}
