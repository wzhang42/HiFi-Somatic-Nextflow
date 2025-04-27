process SVPack{
   publishDir "${params.outdir}/${params.project}/SVPack/${SampleName}", mode: 'copy', overwrite: true
   input:
   tuple(
    val(SampleName),
    path(somatic_sv_vcf)
   )

   output:
   tuple(
    val(SampleName),
    path("${SampleName}.svpack.vcf.gz"),
    path("${SampleName}.svpack.vcf.gz.tbi")
   )
   
   tuple(
    val(SampleName),
    path(somatic_sv_vcf),
    path("${SampleName}.svpack.vcf.gz"),
    path("${SampleName}.svpack.vcf.gz.tbi"),
    emit: To_recover_mate_bnd 
   )
   
   tuple(
    val(SampleName),
    path("${SampleName}.svpack.vcf.gz"),
    path("${SampleName}.svpack.vcf.gz.tbi"),
    emit: To_AnnotSV
   )

   script:
   """
      svpack filter --pass-only --min-svlen ${params.SVPack.svlen} ${somatic_sv_vcf} | \
      svpack match -v - ${params.SVPack.control_vcf}  | \
      svpack consequence - ${params.SVPack.ref_gff} | \
      svpack tagzygosity - > ${SampleName}.svpack.vcf
     
      bgzip ${SampleName}.svpack.vcf

      tabix -p vcf ${SampleName}.svpack.vcf.gz
   """
}

process recover_mate_bnd {
   publishDir "${params.outdir}/${params.project}/recover_mate_bnd/${SampleName}", mode: 'copy', overwrite: true
   input:
   tuple(
     val(SampleName),
     path(somatic_sv_vcf),
     path(somatic_svpack_vcf_gz),
     path(somatic_svpack_vcf_tbi)
   )

   output:
   tuple(
     val(SampleName),
     path("${SampleName}.svpack.recover_mate_bnd.vcf.gz"),
     path("${SampleName}.svpack.recover_mate_bnd.vcf.gz.tbi")
   )
   
   tuple(
     val(SampleName),
     path("${SampleName}.svpack.recover_mate_bnd.vcf.gz"),
     path("${SampleName}.svpack.recover_mate_bnd.vcf.gz.tbi"),
     emit: To_AnnotSV 
   ) 
     
   /* comm -13 \
      <(bcftools query -f '%ID\t%MATE_ID\n' ${somatic_svpack_vcf_gz} | rg -v '\.' | cut -f1 | sort) \
      <(bcftools query -f '%ID\t%MATE_ID\n' ${somatic_svpack_vcf_gz} | rg -v '\.' | cut -f2 | sort) \
      > missing_mate.txt 

     bcftools query -f '%ID\t%MATE_ID\n' ${somatic_svpack_vcf_gz} | rg -v '\.' | cut -f1 | sort > query_sort1
     bcftools query -f '%ID\t%MATE_ID\n' ${somatic_svpack_vcf_gz} | rg -v '\.' | cut -f2 | sort > query_sort2
   */

   script:
   """
     bcftools query -f '%ID\t%MATE_ID\n' ${somatic_svpack_vcf_gz} | grep -v '.' | cut -f1 | sort > query_sort1
     bcftools query -f '%ID\t%MATE_ID\n' ${somatic_svpack_vcf_gz} | grep -v '.' | cut -f2 | sort > query_sort2
     comm -13 query_sort1 query_sort2 > missing_mate.txt
   
     bcftools view \
     -i 'ID=@missing_mate.txt' \
     ${somatic_sv_vcf} |\
     bcftools sort -Oz -o missing_mate.vcf.gz

     bcftools concat \
     ${somatic_svpack_vcf_gz } \
     missing_mate.vcf.gz |\
     bcftools sort -Oz -o tmp.vcf.gz

     mv tmp.vcf.gz ${SampleName}.svpack.recover_mate_bnd.vcf.gz
     tabix -p vcf ${SampleName}.svpack.recover_mate_bnd.vcf.gz
  
   """
}

process AnnotSV {
   publishDir "${params.outdir}/${params.project}/AnnotSV/${SampleName}", mode: 'copy', overwrite: true
   input:
   tuple(
     val(SampleName), 
     path(Somatics_SV_vcf_gz),
     path(Somatics_SV_vcf_tbi)
   )
   
   output:
   tuple(
    val(SampleName),
    path("${SampleName}.annotsv.tsv")
   )

   script:
   """
      mkdir -p annotsv_cache_dir 
      tar -xzvf ${params.AnnotSV.annotsv_cache} -C annotsv_cache_dir/

      AnnotSV \
            -annotationsDir annotsv_cache_dir/AnnotSV/ \
            -SVinputFile ${Somatics_SV_vcf_gz} \
            -outputDir . \
            -outputFile ${SampleName}.annotsv.tsv \
            -SVinputInfo 1 \
            -genomeBuild GRCh38

        # Delete cache after annotation
        rm -rf annotsv_cache_dir/    
   """
}
