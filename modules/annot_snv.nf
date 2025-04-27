process bcftool_norm_Somatic {
    publishDir "${params.outdir}/${params.project}/bcftool_norm_Somatic/${SampleName}", mode: 'copy', overwrite: true, pattern: "*.*" 

    input:
      tuple(
            val(SampleName),
            path(input_vcf_gz),
            path(input_vcf_gz_tbi)
      )
      val(Ref_Fa)  //path(Ref_Fa)
    
    output:
      tuple(
            val(SampleName),
            path("${SampleName}.norm.vcf.gz"),
            path("${SampleName}.norm.vcf.gz.tbi")
       )

      tuple(
            val(SampleName),
            path("${SampleName}.norm.vcf.gz"),
       	    path("${SampleName}.norm.vcf.gz.tbi"),
            emit: To_VEP_Annot_Somatic
      )
     
    // a=5; echo $((a+2)) // will display 7
    script:
    """
      export TMPDIR=${params.TMP_DIR}
      export _JAVA_OPTIONS=-Djava.io.tmpdir=${params.TMP_DIR}
      
      bcftools view ${input_vcf_gz} | sed -e 's/ID=AD,Number=./ID=AD,Number=R/' | \
      bcftools norm --threads ${task.cpus - 1} --multiallelics - --output-type b --fasta-ref ${Ref_Fa} | \
      bcftools sort -Oz -o ${SampleName}.norm.vcf.gz
            
      bcftools index -t ${SampleName}.norm.vcf.gz
    """
}

process VEP_Annot_Somatic {
   publishDir "${params.outdir}/${params.project}/VEP_Annot_Somatic/${SampleName}", mode: 'copy', overwrite: true, pattern: "*.*"

    input:
      tuple(
            val(SampleName),
            path(input_vcf_gz),
            path(input_vcf_gz_tbi)
      )
      val(Ref_Fa)  //path(Ref_Fa)

    output:
      tuple(
            val(SampleName),
            path("${SampleName}.vep.vcf.gz"),
            path("${SampleName}.vep.vcf.gz_summary.html")
       )

      tuple(
            val(SampleName),
            path("${SampleName}.vep.vcf.gz"),
            emit: To_Priorization
      )

    script:
    """
      export TMPDIR=${params.TMP_DIR}
      export _JAVA_OPTIONS=-Djava.io.tmpdir=${params.TMP_DIR}
      
      mkdir -p vep_data
      tar -xzvf ${params.VEP_Annot_Somatic.vep_cache} -C vep_data/
      vep \
            --cache \
            --offline \
            --dir vep_data/ \
            --fasta ${Ref_Fa} \
            --format vcf \
            --fork ${task.cpus} \
            --species homo_sapiens \
            --assembly GRCh38 \
            --symbol \
            --hgvs \
            --refseq \
            --check_existing \
            --vcf \
            --pick \
            --flag_pick_allele_gene \
            --everything \
            --compress_output bgzip \
            -i ${input_vcf_gz} \
            -o ${SampleName}.vep.vcf.gz

        # Delete cache after annotation
        rm -rf vep_data/
    """
}
