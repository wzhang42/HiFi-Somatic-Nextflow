process BAM_QC_mosdepth {
  publishDir "${params.outdir}/${params.project}/BAM_QC_mosdepth/${SampleName}", mode: 'copy', overwrite: true
  input:
    tuple(
            val(SampleName),
            val(Normal_Tumor),
            path(Aligned_BAM),
            path(Aligned_BAM_bai)
     )

   output:
     tuple(
            val(SampleName),
            val(Normal_Tumor),
            path("${SampleName}.${Normal_Tumor}.mosdepth.global.dist.txt"),
            path("${SampleName}.${Normal_Tumor}.mosdepth.region.dist.txt"),
            path("${SampleName}.${Normal_Tumor}.mosdepth.summary.txt"),
            path("${SampleName}.${Normal_Tumor}.regions.bed.gz"),
            path("${SampleName}.${Normal_Tumor}.regions.bed.gz.csi")
      )  
  
  script:
  """
  mosdepth \
  -t ${task.cpus} \
  --by 500 \
  --no-per-base \
  --use-median \
  ${SampleName}.${Normal_Tumor} \
  ${Aligned_BAM}
  """
}

process BAM_QC_bam_coverage {
  publishDir "${params.outdir}/${params.project}/BAM_QC_bam_coverage/${SampleName}", mode: 'copy', overwrite: true
  input:
    tuple(
            val(SampleName),
            val(Normal_Tumor),
            path(Aligned_BAM),
            path(Aligned_BAM_bai)
     )

   output:
     tuple(
            val(SampleName),
            val(Normal_Tumor),
            path("${SampleName}.${Normal_Tumor}.bam_coverage.bw")
      )   

  //"/hpcf/authorized_apps/rhel7_apps/python/install/3.7.0/bin" 
  script:
  """
   ${params.BAM_QC_bam_coverage.bamCoverage_BIN_PATH}/bamCoverage \
   -b ${Aligned_BAM} \
   -o ${SampleName}.${Normal_Tumor}.bam_coverage.bw \
   --numberOfProcessors ${task.cpus} \
   -bs 1
  """
} 
