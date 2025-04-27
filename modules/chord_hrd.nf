process PredictHRD {
publishDir "${params.outdir}/${params.project}/PredictHRD/${SampleName}", mode: 'copy', overwrite: true
  input:
  tuple(
    val(SampleName),
    path(Somatic_SNV_vcf_gz),
    path(Somatic_SNV_vcf_gz_tbi),
    path(Somatic_SV_vcf)
  )  

  output:
  tuple(
    path("chord_hrd.log"),
    path("${SampleName}_chord_prediction.txt"),
    path("${SampleName}_chord_signatures.txt")
  )
  
  script:
  """
    gunzip -c ${Somatic_SNV_vcf_gz} > ${SampleName}_SNV.vcf      

    sed 's/gridss/manta/g' \
        /opt/chord/extractSigPredictHRD.R > ./extractSigPredictHRD.R
    chmod +x ./extractSigPredictHRD.R

    ./extractSigPredictHRD.R . \
    ${SampleName} \
    ${SampleName}_SNV.vcf \
    ${Somatic_SV_vcf} \
    38 2>&1 | tee chord_hrd.log

    rm -f ./extractSigPredictHRD.R
  """
}
