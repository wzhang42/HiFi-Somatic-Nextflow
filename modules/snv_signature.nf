process MutationalPattern {
  publishDir "${params.outdir}/${params.project}/MutationPattern/${SampleName}", mode: 'copy', overwrite: true
  input:
    tuple(
            val(SampleName),
            path(somatic_snv_vcf_gz),
            path(somatic_snv_vcf_gz_tbi)
     )

   output:
     tuple(
            val(SampleName),
            path("${SampleName}.mut_sigs.tsv"),
            path("${SampleName}.reconstructed_sigs.tsv"),
            path("${SampleName}.type_occurences.tsv"),
            path("${SampleName}.mut_sigs_bootstrapped.tsv"),
            path("${SampleName}.mutation_profile.pdf")
      )

  script:
  """
  Rscript --vanilla /app/mutational_pattern.R \
    ${somatic_snv_vcf_gz} \
    ${SampleName} \
    ${params.MutationalPattern.max_delta}
  """
}

