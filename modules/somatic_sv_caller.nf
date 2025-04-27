process Severus_SV {
    publishDir "${params.outdir}/${params.project}/Severus_SV/${SampleName}", mode: 'copy', overwrite: true

    input:
        tuple(
            val(SampleName),
            path(normal_bam),
            path(normal_bam_bai),
            path(tumor_bam),
            path(tumor_bam_bai)
        )

    output:
       
        /*
        tuple(
            val(SampleName),
            path("${SampleName}_breakpoints_double.csv"),
            path("${SampleName}_read_qual.txt"),
            path("${SampleName}_all_SVs.tar.gz"),
            path("${SampleName}_somatic_SVs.tar.gz")
        )
        */

        tuple(
            val(SampleName),
            path("${SampleName}_severus"),
         )

        tuple(
            val(SampleName),
            path("${SampleName}_severus/somatic_SVs/*.vcf"),
            emit: To_PredictHRD
         )
         
         tuple(
            val(SampleName),
            path("${SampleName}_severus/somatic_SVs/*.vcf"),
            emit: To_Annot_Somatic_SV
         )

    
    /* mv breakpoints_double.csv ${SampleName}_breakpoints_double.csv
       mv read_qual.txt ${SampleName}_read_qual.txt
       tar -zvcf ${SampleName}_all_SVs.tar.gz all_SVs
       tar -zvcf ${SampleName}_somatic_SVs.tar.gz  somatic_SVs

       severus \
       --target-bam ${tumor_bam} \
       --control-bam ${normal_bam} \
       --out-dir ${SampleName}_severus \
       -t ${task.cpus} \
       --vntr-bed ${params.Severus.trf_bed} \
       --min-support ${params.Severus.min_supp_reads}      

       mv ${SampleName}_severus/breakpoints_double.csv ${SampleName}_breakpoints_double.csv
       mv ${SampleName}_severus/read_qual.txt ${SampleName}_read_qual.txt
       tar -zvcf ${SampleName}_all_SVs.tar.gz ${SampleName}_severus/all_SVs
       tar -zvcf ${SampleName}_somatic_SVs.tar.gz ${SampleName}_severus/somatic_SVs
    */

    script:
    """
       severus \
       --target-bam ${tumor_bam} \
       --control-bam ${normal_bam} \
       --out-dir ${SampleName}_severus \
       -t ${task.cpus} \
       --vntr-bed ${params.Severus.trf_bed} \
       --min-support ${params.Severus.min_supp_reads}
       
    """
}

