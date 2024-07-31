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
    memory "40 GB"
    cpus 8
    debug true
    tag "${meta_files.sample_id}"

    input:
    tuple val(patient_id), val(meta_files), path(bam), path(bai)
    
    output:
    tuple val(patient_id), 
    val(meta_files), 
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

// Mark duplicates

process mark_duplicates{
    publishDir "${params.pubdir}/${patient_id}/mark_dup/", mode: 'copy'
    container "broadinstitute/gatk:latest"
    memory "8 GB"
    cpus 2
    debug true
    tag "${meta_files.sample_id}"

    input:
    tuple val(patient_id), val(meta_files), path(fastq), path(sam), path(bam), path(sort_mapped_bam)
    
    output:
    tuple val(patient_id), 
    val(meta_files), 
    path("${meta_files.sample_id}_sorted_mapped_MD.bam"), 
    path("${meta_files.sample_id}_marked_dup_metrics.txt"), 
    path("${meta_files.sample_id}_sorted_mapped_MD.bai")

    script:
    """
    echo Mark duplicates running on sample ${meta_files.sample_id}
    gatk --java-options "-Xmx8g -XX:ParallelGCThreads=${task.cpus}" MarkDuplicates \
        -I "${sort_mapped_bam}" \
        -O "${meta_files.sample_id}_sorted_mapped_MD.bam" \
        -M "${meta_files.sample_id}_marked_dup_metrics.txt" \
        --CREATE_INDEX true
    """
}


// somatic SV calling

process delly_call {

    publishDir "${params.pubdir}/${patient_id}/delly_call/", mode: 'copy'
    container "dellytools/delly:latest"
    tag "delly_call_${patient_id}"
    debug true
    memory "8 GB"
    cpus 2

    input:
    tuple val(patient_id), val(normal_meta), path(normal_bam), path(normal_bai), val(tumor_meta), path(tumor_bam), path(tumor_bai)
    
    output:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path("${tumor_meta.sample_id}_vs_${normal_meta.sample_id}.bcf"), 
    path("${tumor_meta.sample_id}_vs_${normal_meta.sample_id}.bcf.csi")

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    delly call -x ${params.reg} -o ${tumor_meta.sample_id}_vs_${normal_meta.sample_id}.bcf -g ${params.ref} $tumor_bam $normal_bam
    """

}

// somatic pre-filtering

process delly_pre_filter {

    publishDir "${params.pubdir}/${patient_id}/delly_prefil", mode: 'copy'
    container "dellytools/delly:latest"
    tag "delly_prefilter_${patient_id}"
    debug true
    memory "8 GB"
    cpus 2
    
    input:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path(bcf), path(bcf_csi)
    
    
    output:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path("${patient_id}.pre.bcf"), path("${patient_id}.pre.bcf.csi")

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}

    echo -e "${normal_meta.sample_id}_sorted_mapped_MD\\t${normal_meta.type}" >> sample.tsv
    echo -e "${tumor_meta.sample_id}_sorted_mapped_MD\\t${tumor_meta.type}" >> sample.tsv

    delly filter -f somatic -o ${patient_id}.pre.bcf -s sample.tsv $bcf 
    """

}   

// genotype pre-filtered somatic sites

process delly_gen {

    publishDir "${params.pubdir}/${patient_id}/delly_gen", mode: 'copy'
    container "dellytools/delly:latest"
    tag "delly_gen_${patient_id}"
    debug true
    memory "8 GB"
    cpus 2

    input:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path(pre_bcf), path(pre_bcf_csi), path(tumor_bam), path(tumor_bai), val(normal_bams), val(normal_bais)
    
    output:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path("${patient_id}.geno.bcf"), path("${patient_id}.geno.bcf.csi")

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    delly call -g ${params.ref} -v ${pre_bcf} -o ${patient_id}.geno.bcf -x ${params.reg} $tumor_bam $normal_bams
    """

}
  

// post-filtering for somatic SVs

process delly_post_fiter {

    publishDir "${params.pubdir}/delly_postfil", mode: 'copy'
    container "quay.io/biocontainers/mulled-v2-aaf75b349de6d380ac6d4a206d51c2d696678b2a:f3c6faf275e70708a7635731117f172a7eafdd14-0"
    tag "delly_post_filter_${patient_id}"
    memory "4 GB"
    
    input:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path(geno_bcf), path(geno_bcf_csi)
        
    output:
    tuple val(patient_id), path("${patient_id}.bcf"), path("${patient_id}.bcf.csi"), path("samples.tsv")

    script:
    """
    sample_names=\$(bcftools query -l "${geno_bcf}")

    echo -e "${tumor_meta.sample_id}_sorted_mapped_MD\\ttumor" > samples.tsv
    
    for sample_name in \${sample_names}; do
        if [ "\${sample_name}" != "${tumor_meta.sample_id}_sorted_mapped_MD" ]; then
            echo -e "\${sample_name}\\tcontrol" >> samples.tsv
        fi
    done
    echo "Content of samples.tsv:"
    cat samples.tsv

    delly filter -f somatic -o ${patient_id}.bcf -s samples.tsv ${geno_bcf}
    """
}

process bcftools{
    publishDir "${params.pubdir}/bcftools", mode: 'copy'
    container "staphb/bcftools:latest"
    tag "bcftools_${patient_id}"
    debug true
    memory '20 MB'

    input:
    tuple val(patient_id), path(bcf), path(bcf_sci)

    output:
    tuple val(patient_id), path("${patient_id}.sv.tsv"), path("${patient_id}.input.tsv")

    script:
    """
    bcftools query -f "%CHROM\t%POS\t%CHROM\t%INFO/END\t%ID\n" ${bcf} | grep -v "BND" > ${patient_id}.sv.tsv 
	bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\n" ${bcf} | grep "BND" >> ${patient_id}.sv.tsv

    cat ${patient_id}.sv.tsv | cut -f 1,2,5 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Left";}' > ${patient_id}.input.tsv
	cat ${patient_id}.sv.tsv | cut -f 3,4,5 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Right";}' >> ${patient_id}.input.tsv
    """
}

process alfred{
    publishDir "${params.pubdir}/alfred"
    memory '70 MB'

    input:
    tuple val(patient_id), path("${patient_id}.sv.tsv"), path("${patient_id}.input.tsv")

    output:
    tuple val(patient_id), path("${patient_id}.sv.gene.tsv")

    container "trausch/alfred:latest"

    script:
    """
    alfred annotate -d 3000 -g /home/user/delly_k8s/Homo_sapiens.GRCh38.107.gtf.gz -o ${patient_id}.sv.gene.tsv ${patient_id}.input.tsv
	rm ${patient_id}.sv.tsv ${patient_id}.input.tsv gene.bed
    """
}

process combine_variants {
    publishDir("/storage/01.NanoBreak/data/samples/${patient_id}/delly_hg38/", mode: 'copy')
    //publishDir "${params.pubdir}/final_tab", mode: 'copy'
    container "quay.io/biocontainers/mulled-v2-76a0b41e09773ed659596514b804bc832021772e:d50820554e481bef847e84d0590124e985594f5d-0"

    input:
    tuple val(patient_id), path("${patient_id}.sv.gene.tsv")

    output:
    path("${patient_id}.combined.sv.tsv")

    script:
    """
    python $params.script --input_file ${patient_id}.sv.gene.tsv --output_file ${patient_id}.combined.sv.tsv
    """
}

process delly_cnv {
    publishDir("${params.pubdir}/CNV", mode: 'copy')
    container "dellytools/delly:latest"
    tag "delly_cnv_${patient_id}"
    debug true
    memory "8 GB"

    input:
    tuple val(patient_id), val(meta), path(bam), path(bai)

    output:
    tuple val(patient_id), path("${patient_id}.cnv.bcf"), path("${patient_id}.cov.gz")

    script:
    """
    delly cnv -a -g ${params.ref} -m ${params.map} -c ${patient_id}.cov.gz -o ${patient_id}.cnv.bcf $bam
    """
}


workflow {
    tsv_ch = Channel.fromPath("/home/user/Delly_somatic/samplesheet_ten.csv")
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
    mark_dup_ch = mark_duplicates(bwa_ch)//.view()

    // test_ch = Channel.fromFilePairs(['/home/user/Delly_somatic/test_bams/*.bam', '/home/user/Delly_somatic/test_bams/*.bai'])//.view()
    // .map{it -> [it[0], it[1][1], it[1][0]]}
    // .view()
    // test_bamToFastq_bwa(test_ch)
 
    tumor_bam_bai = mark_dup_ch
    .filter{ it[1].type == 'tumor' }
    .map{ it -> [it[0], it[1], it[2], it[4]]} 
    // .view()

    normal_bam_bai = mark_dup_ch
    .filter{ it[1].type == 'control' }
    .map{ it -> [it[0], it[1], it[2], it[4]]} 
    // .view()

    sample_ch = normal_bam_bai.cross(tumor_bam_bai)
    .map{it -> [it[0][0], it[0][1], it[0][2], it[0][3], it[1][1], it[1][2], it[1][3]]}
    // .view()

    call_ch = delly_call(sample_ch)
    // .view()

    pre_filter_ch = delly_pre_filter(call_ch)
    //.view()

    control_bams = normal_bam_bai.map { it[2] }.collect(sort: true).toList()//.view()
    control_bais = normal_bam_bai.map { it[3] }.collect(sort: true).toList()
    all_controls = control_bams.combine(control_bais)//.view()
     .map { it ->
            def control_bams_str = it[0].join(' ')
            def control_bais_str = it[1].join(' ')
            [control_bams_str, control_bais_str]
        }
    // .view()

    gen_ch = pre_filter_ch
    .join(tumor_bam_bai)
    .map{it -> [it[0], it[1], it[2], it[3], it[4], it[6], it[7]]}
    .combine(all_controls)
    // .view()
    
    // delly_gen_ch = delly_gen(gen_ch)
    // post_filter_ch = delly_post_fiter(delly_gen_ch)
    //bcf_ch = bcftools(post_filter_ch)
    //alf_ch = alfred(bcf_ch)
    //var_ch = combine_variants(alf_ch)

    cnv_ch = mark_dup_ch.map{ it -> [it[0], it[1], it[2], it[4]]}.view()
    delly_cnv_ch = delly_cnv(cnv_ch.first())
}



