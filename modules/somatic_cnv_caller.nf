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
            path("${SampleName}_cnvkit/${SampleName}*aln.cns"),
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

process CNVKit_Merge_Germline_VCFs {
   publishDir "${params.outdir}/${params.project}/CNVKit_Merge_Germline_VCFs/${SampleName}", mode: 'copy', overwrite: true
   input:
	tuple(
            val(SampleName),
            path(normal_germline_vcf),
            path(normal_germline_vcf_tbi),
            path(tumor_germline_vcf),
            path(tumor_germline_vcf_tbi)
        )

    output:
        tuple(
           val(SampleName),
           path("${SampleName}.merged_germline.vcf.gz"),
           path("${SampleName}.merged_germline.vcf.gz.tbi"),
           path("${SampleName}.merged_germline_heterozygous.withAD.vcf.gz"),
           path("${SampleName}.merged_germline_heterozygous.withAD.vcf.gz.tbi") 
        )     
        
        tuple(
           val(SampleName),
           path("${SampleName}.merged_germline_heterozygous.withAD.vcf.gz"),
           path("${SampleName}.merged_germline_heterozygous.withAD.vcf.gz.tbi"),
           emit: To_CNVKit_Recall 
        )
        
    script:
      """
        echo -e "${SampleName}.normal">sample_names.txt
        echo -e "${SampleName}.tumor">>sample_names.txt
        
        bcftools merge --force-samples ${normal_germline_vcf} ${tumor_germline_vcf} | bcftools reheader -s sample_names.txt -o ${SampleName}.merged_germline.vcf
        
        bgzip ${SampleName}.merged_germline.vcf
        tabix -p vcf ${SampleName}.merged_germline.vcf.gz

        bcftools filter ${SampleName}.merged_germline.vcf.gz -i '(GT[0]=="0/1" || GT[0]=="1/0") && GT[1]!="./."' -Oz -o ${SampleName}.merged_germline_heterozygous.vcf.gz
        tabix -p vcf ${SampleName}.merged_germline_heterozygous.vcf.gz

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%AF\t]\n' ${SampleName}.merged_germline_heterozygous.vcf.gz | awk -v OFS='\t' '{ n_ref_1 = \$7 * (1 - \$8); n_alt_1 = \$7 * \$8; ad_1 = int(n_ref_1+0.5) "," int(n_alt_1+0.5); print \$1,\$2,ad_1}' | bgzip -c > ad_1.txt.gz
        tabix -s1 -b2 -e2 ad_1.txt.gz

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%AF\t]\n' ${SampleName}.merged_germline_heterozygous.vcf.gz | awk -v OFS='\t' '{ n_ref_2 = \$11 * (1 - \$12); n_alt_2 = \$11 * \$12; ad_2 = int(n_ref_2+0.5) "," int(n_alt_2+0.5); print \$1,\$2,ad_2}' | bgzip -c > ad_2.txt.gz
        tabix -s1 -b2 -e2 ad_2.txt.gz

        bcftools annotate -s ${SampleName}.normal -a ad_1.txt.gz -c CHROM,POS,FORMAT/AD ${SampleName}.merged_germline_heterozygous.vcf.gz | bcftools annotate -s ${SampleName}.tumor -a ad_2.txt.gz -c CHROM,POS,FORMAT/AD -Oz -o ${SampleName}.merged_germline_heterozygous.withAD.vcf.gz
  
        tabix -p vcf ${SampleName}.merged_germline_heterozygous.withAD.vcf.gz
      """
}

process CNVKit_Extract_ploidy_purity {
   input:
	tuple(
            val(SampleName),
            path(purity_ploidy)
        )

   output:
       tuple(
         val(SampleName),  
         env(ploidy),
         env(purity),
         emit: To_CNVKit_Recall
       )

    script:
      """
        ploidy=`cut -f5 $purity_ploidy | awk 'NR==2 {print int(\$1 + 0.5)}'`
        purity=`cut -f1 $purity_ploidy|awk 'NR==2 {print \$1}'`
      """
}

process CNVKit_Recall {
   publishDir "${params.outdir}/${params.project}/CNVKit_Recall/${SampleName}", mode: 'copy', overwrite: true   

   input:
	tuple(
            val(SampleName),
            path(cnvkit_cns),
            path(merged_germline_heterozygous_vcf),
            path(merged_germline_heterozygous_vcf_tbi),
            val(ploidy),
            val(purity)
        )

    output:
        tuple(
           val(SampleName),
           path("${SampleName}.CNVKit.with_major_minor_CN.cns")
         )
    
    script:
      """
        cnvkit.py call \
        ${cnvkit_cns} \
        -m clonal \
        --ploidy ${ploidy} \
        --purity ${purity} \
        -o ${SampleName}.CNVKit.with_major_minor_CN.cns \
        -v ${merged_germline_heterozygous_vcf} \
        -i ${SampleName}.tumor \
        -n ${SampleName}.normal \
        --min-variant-depth ${params.CNVKit.min_var_depth}
      """
}

