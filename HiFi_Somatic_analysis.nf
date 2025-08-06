#!/usr/bin/env nextflow
/*******************************************************************************
 This Nextflow pipeline, is for Pacbio longread somatic analysis.  It take paired Somatic Hifi unaligned bam as input. 
  The downstream analysis include pbmm2 mapping, Somatic SNV/INDEL detection (DeepSomatic/Clairs), Structural Variant detection (Severus),
  tumor purity/ploidy inference (Purple/CNVkit), Mutation SNV signature and homologous recombination deficiency estimation.      

 Input: Longread paired somatic HiFi unaligned bam file list 
 Output: 
 Written by: Wenchao Zhang
 The Center for Applied Bioinformatics,St. Jude Children Research Hosptial
 Start Date: 11/26/2024, Final Modify Date: 07/29/2025
*********************************************************************************/
nextflow.enable.dsl=2   
include {Group_UnAligned_BAMs} from './modules/alignment.nf'
include {PBMM2_Align} from './modules/alignment.nf'
include {Pair_Normal_Tumor_Bams} from './modules/alignment.nf'
include {DeepSomatic_SNV} from './modules/somatic_snv_caller.nf'
include {Clair3_SNV} from './modules/somatic_snv_caller.nf'
include {Pair_Normal_Tumor_Vcfs} from './modules/somatic_snv_caller.nf'
include {Severus_SV} from './modules/somatic_sv_caller.nf'
include {BAM_QC_mosdepth} from './modules/bam_qc.nf'
include {BAM_QC_bam_coverage} from './modules/bam_qc.nf'
include {PredictHRD} from './modules/chord_hrd.nf'
include {AMBER} from './modules/clonality.nf'
include {COBALT} from './modules/clonality.nf'
include {PURPLE} from './modules/clonality.nf'
include {MutationalPattern} from './modules/snv_signature.nf'
include {HiFiPhase} from './modules/phasing.nf'
include {HiFiPhase_Somatic} from './modules/phasing.nf'
include {pb_cpg_score as Normal_pb_cpg_score} from './modules/pb_cpg.nf'
include {pb_cpg_score as Tumor_pb_cpg_score} from './modules/pb_cpg.nf'
include {DSS_DMR} from './modules/DMR.nf'
include {DMR_annotation} from './modules/DMR.nf'
include {bcftool_norm_Somatic} from './modules/annot_snv.nf'
include {VEP_Annot_Somatic} from './modules/annot_snv.nf'
include {SVPack} from './modules/annot_sv.nf'
// include {recover_mate_bnd} from './modules/annot_sv.nf'
include {AnnotSV} from './modules/annot_sv.nf'
include {CNVKit_Somatic} from './modules/somatic_cnv_caller.nf'
include {CNVKit_Merge_Germline_VCFs} from './modules/somatic_cnv_caller.nf'
include {CNVKit_Extract_ploidy_purity} from './modules/somatic_cnv_caller.nf'
include {CNVKit_Recall} from './modules/somatic_cnv_caller.nf'

REFERENCE=""
def parse_config_parameters() {
// Parse Nextflow Envrionment Paramters 

  if( params.genome_build == 'hg38' )
  {
     REFERENCE = params.Build_hg38.REFERENCE
     params.AMBER.REF_GENOME_VER = 'GRCh38'
     params.PURPLE.REF_GENOME_VER= '38'
  }
  else if( params.genome_build == 'hg19' )
  {
     REFERENCE = params.Build_hg19.REFERENCE
     params.AMBER.REF_GENOME_VER = 'hs37d5'
     params.PURPLE.REF_GENOME_VER= '19'
  }
  else if( params.genome_build == 'b37' )
  {
     REFERENCE = params.Build_b37.REFERENCE
     params.AMBER.REF_GENOME_VER = 'GRCh37-lite'
     params.PURPLE.REF_GENOME_VER= '37'

  }
  else if( params.genome_build == 'T2T' )
  {
     REFERENCE = params.Build_T2T.REFERENCE
     params.AMBER.REF_GENOME_VER = 'chm13v2'
     params.PURPLE.REF_GENOME_VER= 'T2T'
  }
  else
  {
     error "Invalid geome version: ${params.genome_build}"
     exit 1
  }
}

def DispConfig() {
 log.info """
  Welocme to run Nextflow Pipeline HiFi_Somatic_analysis.nf  
  Your configuration are the following:
  project                           : ${params.project}
  HiFi_Somatic_SampleFilelist       : ${params.HiFi_Somatic_SampleFilelist}
  Select_BAM_QC_mosdepth            : ${params.Select_BAM_QC_mosdepth}
  Select_BAM_QC_bam_coverage        : ${params.Select_BAM_QC_bam_coverage}
  Select_DeepSomatic                : ${params.Select_DeepSomatic}
  Select_Clair3_SNV                 : ${params.Select_Clair3_SNV}  
  Select_Severus_SV                 : ${params.Select_Severus_SV}
  Severus_trf_bed                   : ${params.Severus.trf_bed}
  Severus_min_supp_reads            : ${params.Severus.min_supp_reads}
  Select_PredictHRD                 : ${params.Select_PredictHRD}
  Select_MutationalPattern          : ${params.Select_MutationalPattern}
  MutationalPattern_max_delta       : ${params.MutationalPattern.max_delta} 
 
  ENSEMBLE_DATA_DIR                 : ${params.ENSEMBLE_DATA_DIR}
  Select_AMBER                      : ${params.Select_AMBER}
  Select_COBALT                     : ${params.Select_COBALT}
  Select_PURPLE                     : ${params.Select_PURPLE}

  Select_HiFiPhase                  : ${params.Select_HiFiPhase}
  Select_HiFiPhase_Somatic          : ${params.Select_HiFiPhase_Somatic}
  Select_Tumor_pb_cpg               : ${params.Select_Tumor_pb_cpg}
  Select_Normal_pb_cpg              : ${params.Select_Normal_pb_cpg}
  pb_cpg_min-mapq                   : ${params.pb_cpg_score.minmapq}
  pb_cpg_min-coverage               : ${params.pb_cpg_score.mincov}
  pb_cpg_mode                       : ${params.pb_cpg_score.model}

  Select_DSS_DMR                    : ${params.Select_DSS_DMR}
  Select_DMR_annotation             : ${params.Select_DMR_annotation}

  Select_Annot_SNV                  : ${params.Select_Annot_SNV}
  vep_cache                         : ${params.VEP_Annot_Somatic.vep_cache}

  Select_Annot_SV                   : ${params.Select_Annot_SV}
  SVPack_svlen                      : ${params.SVPack.svlen}
  SVPack_control_vcf                : ${params.SVPack.control_vcf}
  SVPack_ref_gff                    : ${params.SVPack.ref_gff}
  annotsv_cache                     : ${params.AnnotSV.annotsv_cache}

  Select_CNVKit_Somatic             : ${params.Select_CNVKit_Somatic}
  CNVKit_refFlat                    : ${params.CNVKit.refFlat}
  CNVKit_target_avg_size            : ${params.CNVKit.target_avg_size} 

  outdir                            : ${params.outdir} 
 """  
//  exit 0    
}

workflow {
   //Parse the input Parameters and configs. 
   parse_config_parameters() 
   // Display the configuration
   DispConfig() 
   
   // Set up the the external input channels for samples (fastq files) 
   UnalignBAM_Ch    = Channel
            .fromPath(params.HiFi_Somatic_SampleFilelist)
            .splitText()
            .splitCsv(sep: '\t')
  
  Group_UnAligned_BAM_Ch=Group_UnAligned_BAMs(UnalignBAM_Ch)
  PBMM2Align_Ch= PBMM2_Align(Group_UnAligned_BAM_Ch.To_GroupBAM.groupTuple(by:2), REFERENCE)
   
  BAM_QC_mosdepth_InCh = (params.Select_BAM_QC_mosdepth == "Y" )? PBMM2Align_Ch.To_BAM_QC_mosdepth : Channel.empty()  
  BAM_QC_mosdepth_Ch   = BAM_QC_mosdepth(BAM_QC_mosdepth_InCh) 

  BAM_QC_bam_coverage_InCh = (params.Select_BAM_QC_bam_coverage == "Y" )? PBMM2Align_Ch.To_BAM_QC_bam_coverage : Channel.empty()
  BAM_QC_bam_coverage_Ch   = BAM_QC_bam_coverage(BAM_QC_bam_coverage_InCh)
  
  Clair3_SNV_InCh = (params.Select_Clair3_SNV == "Y" )? PBMM2Align_Ch.To_Clair3_SNV : Channel.empty()
  Clair3_SNV_Ch   = Clair3_SNV(Clair3_SNV_InCh, REFERENCE)  
    
  Pair_Normal_Tumor_Bam_Ch =Pair_Normal_Tumor_Bams(PBMM2Align_Ch.To_Pair_Normal_Tumor.groupTuple())

  DeepSomatic_SNV_InCh = (params.Select_DeepSomatic == "Y" )? Pair_Normal_Tumor_Bam_Ch.To_DeepSomatic_SNV : Channel.empty()
  DeepSomatic_SNV_Ch = DeepSomatic_SNV(DeepSomatic_SNV_InCh, REFERENCE)
  
  Severus_SV_InCh = (params.Select_Severus_SV == "Y" )? Pair_Normal_Tumor_Bam_Ch.To_Severus_SV : Channel.empty()
  Severus_SV_Ch = Severus_SV(Severus_SV_InCh) 

  MutationalPattern_InCh = (params.Select_MutationalPattern == "Y" )? DeepSomatic_SNV_Ch.To_MutationalPattern : Channel.empty()
  MutationalPattern_Ch = MutationalPattern(MutationalPattern_InCh)

  SNV_To_PredictHRD = DeepSomatic_SNV_Ch.To_PredictHRD
  SV_To_PredictHRD  = Severus_SV_Ch.To_PredictHRD   
  PredictHRD_InCh = (params.Select_PredictHRD == "Y")? SNV_To_PredictHRD.combine(SV_To_PredictHRD, by:0 ): Channel.empty()
  PredictHRD_Ch = PredictHRD(PredictHRD_InCh)    

  AMBER_InCh = (params.Select_AMBER=="Y" )? Pair_Normal_Tumor_Bam_Ch.To_AMBER : Channel.empty()
  AMBER_Ch   = AMBER(AMBER_InCh, REFERENCE)

  COBALT_InCh= (params.Select_COBALT=="Y")? Pair_Normal_Tumor_Bam_Ch.To_COBALT : Channel.empty()  
  COBALT_Ch  = COBALT(COBALT_InCh, REFERENCE)

  AMBER_To_PURPLE= (params.Select_AMBER=="Y")? AMBER_Ch.To_Purple : Channel.empty()
  COLBALT_To_PURPLE= (params.Select_COBALT=="Y")? COBALT_Ch.To_Purple : Channel.empty() 
  
  PURPLE_InCh = (params.Select_PURPLE=="Y")? AMBER_To_PURPLE.combine(COLBALT_To_PURPLE, by:0) : Channel.empty()
  PURPLE_Ch = PURPLE(PURPLE_InCh, REFERENCE)

  Pair_Normal_Tumor_Vcf_Ch =Pair_Normal_Tumor_Vcfs(Clair3_SNV_Ch.To_Pair_Normal_Tumor.groupTuple())
  
  Normal_AlignedBam_Ch   = Pair_Normal_Tumor_Bam_Ch.To_HiPhasing_Normal_Bam
  Normal_Germline_Vcf_Ch = Pair_Normal_Tumor_Vcf_Ch.To_HiPhasing_Normal_Vcf
  HiFiPhase_Normal_InCh =(params.Select_HiFiPhase =="Y")? Normal_AlignedBam_Ch.combine(Normal_Germline_Vcf_Ch, by:0) : Channel.empty()
  HiFiPhase_Normal_Ch =HiFiPhase(HiFiPhase_Normal_InCh, REFERENCE)

  Tumor_AlignedBam_Ch   = Pair_Normal_Tumor_Bam_Ch.To_HiPhasing_Tumor_Bam
  Tumor_Germline_Vcf_Ch = Pair_Normal_Tumor_Vcf_Ch.To_HiPhasing_Tumor_Vcf
  Somatic_SNV_Vcf_Ch = DeepSomatic_SNV_Ch.To_Hiphasing_Somatic_Vcf
  HiFiPhase_Somatic_InCh =(params.Select_HiFiPhase_Somatic =="Y")? (Tumor_AlignedBam_Ch.combine(Tumor_Germline_Vcf_Ch, by:0)).combine(Somatic_SNV_Vcf_Ch, by:0) : Channel.empty()
  HiFiPhase_Somatic_Ch =HiFiPhase_Somatic(HiFiPhase_Somatic_InCh, REFERENCE)
  
  Tumor_pb_cpg_InCh = (params.Select_Tumor_pb_cpg =="Y")? HiFiPhase_Somatic_Ch.To_Tumor_pb_cpg : Channel.empty()
  Tumor_pb_cpg_Ch = Tumor_pb_cpg_score(Tumor_pb_cpg_InCh, REFERENCE)  

  Normal_pb_cpg_InCh = (params.Select_Normal_pb_cpg =="Y")? HiFiPhase_Normal_Ch.To_Normal_pb_cpg : Channel.empty()
  Normal_pb_cpg_Ch    = Normal_pb_cpg_score(Normal_pb_cpg_InCh, REFERENCE)
    
  DSS_DMR_InCh = (params.Select_DSS_DMR =="Y")? (Tumor_pb_cpg_Ch.To_DMR).combine(Normal_pb_cpg_Ch.To_DMR, by:0) : Channel.empty()
  DSS_DMR_Ch = DSS_DMR(DSS_DMR_InCh)

  DMR_annotation_InCh = (params.Select_DMR_annotation =="Y")? DSS_DMR_Ch.To_DMR_Annotation : Channel.empty()  
  DMR_annotation_Ch = DMR_annotation(DMR_annotation_InCh)
  
  bcftool_norm_Somatic_InCh= (params.Select_Annot_SNV =="Y")? DeepSomatic_SNV_Ch.To_Annot_Somatic : Channel.empty()
  bcftool_norm_Somatic_Ch = bcftool_norm_Somatic (bcftool_norm_Somatic_InCh, REFERENCE ) 
  
  VEP_Annot_Somatic_InCh = (params.Select_Annot_SNV =="Y")? bcftool_norm_Somatic_Ch.To_VEP_Annot_Somatic : Channel.empty()
  VEP_Annot_Somatic_Ch = VEP_Annot_Somatic( VEP_Annot_Somatic_InCh, REFERENCE)

  SVPack_InCh =(params.Select_Annot_SV =="Y")? Severus_SV_Ch.To_Annot_Somatic_SV : Channel.empty()
  SVPack_Ch = SVPack(SVPack_InCh)

  AnnotSV_InCh = (params.Select_Annot_SV =="Y")? SVPack_Ch.To_AnnotSV : Channel.empty()
  AnnotSV_Ch  =AnnotSV(AnnotSV_InCh)
  
  CNVKit_Somatic_InCh = (params.Select_CNVKit_Somatic == "Y" )? Pair_Normal_Tumor_Bam_Ch.To_CNVKit_Somatic : Channel.empty()
  CNVKit_Somatic_Ch = CNVKit_Somatic(CNVKit_Somatic_InCh, REFERENCE)
 
  CNVKit_Merge_Germline_VCFs_InCh = (params.Select_CNVKit_Somatic == "Y" )? Pair_Normal_Tumor_Vcf_Ch.To_CNVKit_Merge_Germline_VCFs : Channel.empty()
  CNVKit_Merge_Germline_VCFs_Ch = CNVKit_Merge_Germline_VCFs(CNVKit_Merge_Germline_VCFs_InCh)

  CNVKit_Extract_ploidy_purity_Ch = CNVKit_Extract_ploidy_purity(PURPLE_Ch.To_CNVKit_Extract_ploidy_purity)

  CNVKit_Recall_InCh = (params.Select_CNVKit_Somatic == "Y" )? (CNVKit_Somatic_Ch.To_CNVKit_Recall.combine(CNVKit_Merge_Germline_VCFs_Ch.To_CNVKit_Recall, by:0)).combine(CNVKit_Extract_ploidy_purity_Ch.To_CNVKit_Recall,by:0) : Channel.empty() 
  CNVKit_Recall_Ch = CNVKit_Recall(CNVKit_Recall_InCh)   
}
