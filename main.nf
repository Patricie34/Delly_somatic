nextflow.enable.dsl=2

// TEST bam_to_fastq
process test_bamToFastq_bwa{
    publishDir "${params.pubdir}/test/", mode: 'copy'
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:23592e4ad15ca2acfca18facab87a1ce22c49da1-0"
    memory "7 GB"
    debug true
    

    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.fastq"), path("${sample_id}.sam"), path("${sample_id}.bam")

    script:
    """
    samtools bam2fq ${bam} > ${sample_id}.fastq
    bwa mem ${params.ref_bwa} ${sample_id}.fastq > ${sample_id}.sam
    samtools view -bS ${sample_id}.sam > ${sample_id}.bam
    """

}


// Convert bams to fastq, fastq to bam, sort
process samtools_bwa{
    publishDir "${params.pubdir}/${patient_id}/bwa/", mode: 'copy'
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:23592e4ad15ca2acfca18facab87a1ce22c49da1-0"
    memory "32 GB"
    cpus 8
    debug true
    tag "${meta_files.sample_id}"

    input:
    tuple val(patient_id), val(meta_files), path(bam), path(bai)
    
    output:
    tuple val(patient_id), 
    path("${meta_files.sample_id}.fastq"), 
    path("${meta_files.sample_id}.sam"), 
    path("${meta_files.sample_id}.bam"), 
    path("${meta_files.sample_id}_sorted.bam")

    script:
    """
    echo "Samtools bam2fq is running on sample ${meta_files.sample_id}"
    samtools bam2fq ${bam} > ${meta_files.sample_id}.fastq

    echo "BWA is running on sample ${meta_files.sample_id}"
    bwa mem -t ${task.cpus} ${params.ref_bwa} ${meta_files.sample_id}.fastq > ${meta_files.sample_id}.sam

    echo "Samtools view is running on sample ${meta_files.sample_id}"
    samtools view -@ ${task.cpus} -bS ${meta_files.sample_id}.sam > ${meta_files.sample_id}.bam

    echo "Samtools sort is running on sample ${meta_files.sample_id}"
    samtools sort -@ ${task.cpus} ${meta_files.sample_id}.bam -o ${meta_files.sample_id}_sorted.bam
    """
}

// 2nd option bam reheader

// process reheader{
//     publishDir "${params.pubdir}/fastq_bwa/", mode: 'copy'
//     container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:23592e4ad15ca2acfca18facab87a1ce22c49da1-0"
//     memory "10 GB"
//     cpus 16
//     debug true

//     input:
//     tuple val(patient_id), val(meta_files), path(bam), path(bai)
    
//     output:
//     tuple val(patient_id), 
//     path("${meta_files.sample_id}_rehead.bam")

//     script:
//     """
//     samtools reheader -c 'perl -pe "s/^(@SQ.*)(\tSN:)Chr/\$1\$2/"' ${bam}

//     """

// }



// Mark duplicates

process mark_duplicates{
    publishDir "${params.pubdir}/${patient_id}/mark_dup/", mode: 'copy'
    container "broadinstitute/gatk:latest"
    memory "16 GB"
    cpus 8
    debug true
    tag "${meta_files.sample_id}"

    input:
    tuple val(patient_id), val(meta_files), path(sort_mapped_bam)
    
    output:
    tuple val(patient_id), path("${meta_files.sample_id}_sorted_mapped_MD.bam"), path("${meta_files.sample_id}_marked_dup_metrics.txt")

    script:
    """
    echo Mark duplicates running on sample ${meta_files.sample_id}
    gatk --java-options "-Xmx12g -XX:ParallelGCThreads=${task.cpus}" MarkDuplicates \
        -I "${sort_mapped_bam}" \
        -O "${meta_files.sample_id}_sorted_mapped_MD.bam" \
        -M "${meta_files.sample_id}_marked_dup_metrics.txt" \
        --CREATE_INDEX
    """
}


// somatic SV calling

process delly_call {

    publishDir "${params.pubdir}/results/", mode: 'copy'
    container "dellytools/delly:latest"

    input:
    tuple val(patient_id), val(normal_meta), val(normal_files), val(tumor_meta), val(tumor_files)
    
    output:
    path '*'

    script:
    """

    delly call -x ${params.reg} -o ${patient_id}_SV_delly.bcf -g ${params.ref} $tumor_files.bam $normal_files.bam 
    
    """

}

// delly call -x hg38.excl -o t1.bcf -g hg38.fa tumor1.bam control1.bam


// somatic pre-filtering

process delly_pre_filter {

    publishDir "${params.pubdir}/results/SV_delly/filter", mode: 'copy'
    container "dellytools/delly:latest"

    input:
    tuple val(sample), path(bcf), path(bcf_csi)
    
    output:
    tuple val(sample), path("${sample}_SV_delly.pre.bcf")

    script:
    """
    delly filter -f somatic -o ${sample}_SV_delly.pre.bcf -s samples.tsv ${sample}_SV_delly.bcf 
    """

}

// genotype pre-filtered somatic sites

process delly_gen {

    publishDir "${params.pubdir}/results/SV_delly", mode: 'copy'
    container "dellytools/delly:latest"

    input:
    tuple val(sample), path(pre_bcf), path(t_bam), path(t_bai), path(n_bam), path(n_bai)
    
    output:
    tuple val(sample), path("${sample}_SV_delly.geno.bcf")

    script:
    """
    delly call -g hg38.fa -v ${pre_bcf} -o ${sample}_SV_delly.geno.bcf -x hg38.excl ${t_bam} ${n_bam}
    """

}

// post-filtering for somatic SVs

process delly_post_fiter {

    publishDir "${params.pubdir}/results/SV_delly", mode: 'copy'
    container "dellytools/delly:latest"

    input:
    tuple val(sample), path(geno_bcf)
        
    
    output:
    tuple val(sample), path("${sample}_SV_delly.somatic.bcf")

    script:
    """
    delly filter -f somatic -o ${sample}_SV_delly.somatic.bcf -s samples.tsv ${geno.bcf}
    """

}


workflow {
    tsv_ch = Channel.fromPath("/home/user/Delly_somatic/samples_first.csv")
    . splitCsv( header: true)
    . map { row -> 
    def path_bam = file("/storage/delly_hg38_allSamples_BAM/${row.sample_id}.recal.bam")
    def path_bai = file("/storage/delly_hg38_allSamples_BAM/${row.sample_id}.recal.bai")
    [[patient_id: row.patient_id, sample_id: row.sample_id, type: row.type], [bam: path_bam, bai: path_bai]]
    // [row.sample_id, t_bam, t_bai, n_bam, n_bai]

    } 
    // .view()

    fastq_ch = tsv_ch.map{ it -> [it[0].patient_id, it[0], it[1].bam, it[1].bai]}
    // .view() 
    
    // samtools_bwa(fastq_ch)
    
    bwa_ch = samtools_bwa(fastq_ch)
    mark_dup_ch = mark_duplicates(bwa_ch)


    tumor_bam_bai = tsv_ch
    .filter{ it[0].type == 'tumor' }
    .map{ it -> [it[0].patient_id, it[0], it[1]]} 
    // .view()


    normal_bam_bai = tsv_ch
    .filter{ it[0].type == 'normal' }
    .map{ it -> [it[0].patient_id, it[0], it[1]]} 
    //.view()


    sample_ch = normal_bam_bai.cross(tumor_bam_bai).map{it -> [it[0][0], it[0][1], it[0][2], it[1][1], it[1][2]]}//.view()
    // .view()
    call_ch = delly_call(sample_ch.first())

}



