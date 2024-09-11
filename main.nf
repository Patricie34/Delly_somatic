include { test_bamToFastq_bwa; SAMTOOLS_BWA; MARK_DUPLICATES; DELLY_CALL; DELLY_PREFILTER; DELLY_GEN;
DELLY_POSTFILTER; BCFTOOLS; ALFRED; MERGE_TSV; DELLY_CNV; CNV_UNZIP; CNV_PROFILES; CNV_CALLS; CNV_CALL_PLT  } from "${params.projectDirectory}/modules"

workflow {
    tsv_ch = Channel.fromPath("/storage2/delly/Delly_somatic/samples_all_WGS.csv")
    . splitCsv( header: true)
    . map { row -> 
    def path_bam = file("/storage/delly_hg38_allSamples_BAM/${row.sample_id}.recal.bam")
    def path_bai = file("/storage/delly_hg38_allSamples_BAM/${row.sample_id}.recal.bai")
    [[patient_id: row.patient_id, sample_id: row.sample_id, type: row.type, cont: row.cont], [bam: path_bam, bai: path_bai]]
    } 
    // .view()
        
    fastq_ch = tsv_ch.map{ it -> [it[0].patient_id, it[0], it[1].bam, it[1].bai]}
    // .view() 
            
    bwa_ch = SAMTOOLS_BWA(fastq_ch)
    mark_dup_ch = MARK_DUPLICATES(bwa_ch)//.view()

    // test_ch = Channel.fromFilePairs(['/home/user/Delly_somatic/test_bams/*.bam', '/home/user/Delly_somatic/test_bams/*.bai'])//.view()
    // .map{it -> [it[0], it[1][1], it[1][0]]}
    // .view()
    // test_bamToFastq_bwa(test_ch)
 
    tumor_bam_bai = mark_dup_ch
    .filter{ it[1].type == 'tumor' }
    .map{ it -> [it[0], it[1], it[2], it[3]]} 
    //.view()

    normal_bam_bai = mark_dup_ch
    .filter{ it[1].type == 'control' }
    .map{ it -> [it[0], it[1], it[2], it[3]]} 
    .view()

    sample_ch = normal_bam_bai.cross(tumor_bam_bai)//.view()
    .map{it -> [it[0][0], it[0][1], it[0][2], it[0][3], it[1][1], it[1][2], it[1][3]]}
    // .view()

    call_ch = DELLY_CALL(sample_ch)
    // .view()

    
    pre_filter_ch = DELLY_PREFILTER(call_ch)
    // .view()

    control_bams = normal_bam_bai.map { it[2] }.collect(sort: true).toList()
    // .view()
    control_bais = normal_bam_bai.map { it[3] }.collect(sort: true).toList()
    all_controls = control_bams.combine(control_bais)
    // .view()
     .map { it ->
            def control_bams_str = it[0].join(' ')
            def control_bais_str = it[1].join(' ')
            [control_bams_str, control_bais_str]
        }
    // .view()

    tumor_ch = tumor_bam_bai
    .map { it -> [it[0], it[1].sample_id, it[1], it[2], it[3]] }
    // .view()

    gen_ch = pre_filter_ch
    .map{it -> [it[0], it[2].sample_id, it[1], it[2], it[3], it[4]] } 
    .join(tumor_ch, by: [0, 1])
    // .view()
    .combine(all_controls)
    .map{it -> [it[0], it[1], it[2], it[3], it[4], it[5], it[7], it[8], it[9], it[10]]} 
    // .view()

    delly_gen_ch = DELLY_GEN(gen_ch)
    // .view()
    post_filter_ch = DELLY_POSTFILTER(delly_gen_ch)
    // .view()
    // bcf_ch = BCFTOOLS(post_filter_ch)
    // alf_ch = ALFRED(bcf_ch)
    // var_ch = MERGE_TSV(alf_ch)

    // cnv_ch = mark_dup_ch.map{ it -> [it[0], it[1], it[2], it[4]]}.view()
    // delly_cnv_ch = DELLY_CNV(cnv_ch)
    // .filter { it[0] == "BRNO0501" }
    // .view()

    // // this runs only ones
    // unzip_ch = CNV_UNZIP(delly_cnv_ch)

    // plot_ch = CNV_PROFILES(delly_cnv_ch) //(delly_cnv_ch.first())

    // // optional
    // bcfcnv_ch = CNV_CALLS(delly_cnv_ch)
    // cnvPlt_ch = CNV_CALL_PLT (bcfcnv_ch)

}