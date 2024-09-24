nextflow.enable.dsl=2

// TEST bam_to_fastq
process test_bamToFastq_bwa{
    publishDir "${params.pubdir}/test/", mode: 'copy'
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:23592e4ad15ca2acfca18facab87a1ce22c49da1-0"
    label "m_mem"
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

// Convert bams to fastq, fastq to bam & sort
process SAMTOOLS_BWA {
    // publishDir "${params.pubdir}/${patient_id}/bwa/", mode: 'copy'
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:23592e4ad15ca2acfca18facab87a1ce22c49da1-0"
    tag "samtools_bwa on ${patient_id}_${meta.sample_id}"
    debug true
    label "l_mem"
    
    input:
    tuple val(patient_id), val(meta), path(bam), path(bai)
    
    output:
    tuple val(patient_id), val(meta), path("${meta.sample_id}_sorted.bam")

    script:
    """
    echo SAMTOOLS_BWA ${patient_id}_${meta.sample_id}
    samtools bam2fq ${bam} > ${meta.sample_id}.fastq
    bwa mem -t ${task.cpus} ${params.ref_bwa} ${meta.sample_id}.fastq > ${meta.sample_id}.sam
    samtools view -@ ${task.cpus} -bS ${meta.sample_id}.sam > ${meta.sample_id}.bam
    samtools sort -@ ${task.cpus} ${meta.sample_id}.bam -o ${meta.sample_id}_sorted.bam
    """
}

// Mark duplicates & index
process MARK_DUPLICATES {
    publishDir "${params.pubdir2}/${patient_id}/delly_hg38/", mode: 'copy'
    container "broadinstitute/gatk:latest"
    tag "MARK_DUPLICATES on ${patient_id}_${meta.sample_id}"
    debug true
    label "m_mem"
        
    input:
    tuple val(patient_id), val(meta), path(sort_mapped_bam)
    
    output:
    tuple val(patient_id), 
    val(meta), 
    path("${meta.sample_id}_sorted_mapped_MD.bam"), 
    path("${meta.sample_id}_sorted_mapped_MD.bai")

    script:
    """
    echo MARK_DUPLICATES ${patient_id}_${meta.sample_id}
    gatk --java-options "-Xmx8g -XX:ParallelGCThreads=${task.cpus}" MarkDuplicates \
        -I "${sort_mapped_bam}" \
        -O "${meta.sample_id}_sorted_mapped_MD.bam" \
        -M "${meta.sample_id}_marked_dup_metrics.txt" \
        --CREATE_INDEX true
    """
}

// Somatic SV calling
process DELLY_CALL {
    // publishDir "${params.pubdir}/${patient_id}/delly_call/", mode: 'copy'
    container "dellytools/delly:latest"
    tag "DELLY_CALL on ${patient_id}"
    debug true
    label "m_mem"
    
    input:
    tuple val(patient_id), val(normal_meta), path(normal_bam), path(normal_bai), val(tumor_meta), path(tumor_bam), path(tumor_bai)
    
    output:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path("${tumor_meta.sample_id}_vs_${normal_meta.sample_id}.bcf"), 
    path("${tumor_meta.sample_id}_vs_${normal_meta.sample_id}.bcf.csi")

    script:
    """
    echo DELLY_CALL ${patient_id}
    export OMP_NUM_THREADS=${task.cpus}
    delly call -x ${params.reg} -o ${tumor_meta.sample_id}_vs_${normal_meta.sample_id}.bcf -g ${params.ref} $tumor_bam $normal_bam
    """
}

// Somatic pre-filtering
process DELLY_PREFILTER {
    publishDir "${params.pubdir}/${patient_id}/delly_prefil", mode: 'copy'
    container "dellytools/delly:latest"
    tag "DELLY_PREFILTER on ${patient_id}_${tumor_meta.sample_id}"
    debug true
    label "xxs_mem"
        
    input:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path(bcf), path(bcf_csi)
        
    output:
    tuple val(patient_id), val(normal_meta), val(tumor_meta), path("${patient_id}_${tumor_meta.sample_id}.pre.bcf"), path("${patient_id}_${tumor_meta.sample_id}.pre.bcf.csi")//, path("sample.tsv")

    script:
    """
    echo DELLY_PREFILTER ${patient_id}
    export OMP_NUM_THREADS=${task.cpus}
    echo -e "${normal_meta.sample_id}_sorted_mapped_MD\\t${normal_meta.type}\\t${normal_meta.cont}" >> sample.tsv
    echo -e "${tumor_meta.sample_id}_sorted_mapped_MD\\t${tumor_meta.type}\\t${normal_meta.cont}" >> sample.tsv
    if [[ "${normal_meta.type}" == "control" ]]; then
        cont_value=${normal_meta.cont}
    elif [[ "${tumor_meta.type}" == "tumor" ]]; then
        cont_value=${tumor_meta.cont}
    fi
    delly filter -f somatic -c \$cont_value -o ${patient_id}_${tumor_meta.sample_id}.pre.bcf -s sample.tsv $bcf 
    """
}   

// Genotype pre-filtered somatic sites
process DELLY_GEN {
    publishDir "${params.pubdir}/${patient_id}/delly_gen", mode: 'copy'
    container "dellytools/delly:latest"
    tag "DELLY_GEN ON ${patient_id}_${sample_id}"
    debug true
    label "m_mem"
    
    input:
    tuple val(patient_id), val(sample_id), val(normal_meta), val(tumor_meta), path(pre_bcf), path(pre_bcf_csi), path(tumor_bam), path(tumor_bai), val(normal_bams), val(normal_bais)
    
    output:
    tuple val(patient_id), val(sample_id), val(normal_meta), val(tumor_meta), path("${patient_id}_${sample_id}.geno.bcf"), path("${patient_id}_${sample_id}.geno.bcf.csi")

    script:
    """
    echo DELLY_GEN ${patient_id}
    export OMP_NUM_THREADS=${task.cpus}
    delly call -g ${params.ref} -v ${pre_bcf} -o ${patient_id}_${sample_id}.geno.bcf -x ${params.reg} $tumor_bam $normal_bams
    """
}
  
// Post-filtering for somatic SVs
process DELLY_POSTFILTER {
    publishDir "${params.pubdir2}/${patient_id}/delly_hg38/", mode: 'copy'
    publishDir "${params.pubdir}/${patient_id}/delly_postfil", mode: 'copy'
    container "quay.io/biocontainers/mulled-v2-aaf75b349de6d380ac6d4a206d51c2d696678b2a:f3c6faf275e70708a7635731117f172a7eafdd14-0"
    tag "DELLY_POSTFILTER ON ${patient_id}_${sample_id}"
    label "s_mem"
        
    input:
    tuple val(patient_id), val(sample_id), val(normal_meta), val(tumor_meta), path(geno_bcf), path(geno_bcf_csi)
        
    output:
    tuple val(patient_id), val(sample_id), path("${patient_id}_${sample_id}.bcf"), path("${patient_id}_${sample_id}.bcf.csi")//, path("samples.tsv")

    script:
    """
    echo DELLY_POSTFILTER ${patient_id}
    sample_names=\$(bcftools query -l "${geno_bcf}")
    echo -e "${sample_id}_sorted_mapped_MD\\ttumor" > samples.tsv
    
    for sample_name in \${sample_names}; do
        if [ "\${sample_name}" != "${sample_id}_sorted_mapped_MD" ]; then
            echo -e "\${sample_name}\\tcontrol" >> samples.tsv
        fi
    done

    delly filter -a 0.02 -f somatic -o ${patient_id}_${sample_id}.bcf -s samples.tsv ${geno_bcf}
    """
}

// Parsing bcf files
process BCFTOOLS {
    publishDir "${params.pubdir}/${patient_id}/bcftools", mode: 'copy'
    container "staphb/bcftools:latest"
    tag "bcftools_${patient_id}_${sample_id}"
    debug true
    errorStrategy "ignore"
    label "xxs_mem"

    input:
    tuple val(patient_id), val(sample_id), path(bcf), path(bcf_sci)//, path(tsv)

    output:
    tuple val(patient_id), val(sample_id), path("${patient_id}_${sample_id}.sv.tsv"), path("${patient_id}_${sample_id}.input.tsv")

    script:
    """
    bcftools query -f "%CHROM\t%POS\t%CHROM\t%INFO/END\t%ID\t%INFO/SR\t%INFO/CT[\t%RC][\t%RCL][\t%RCR]\n" ${bcf} | grep -v "BND" > ${patient_id}_${sample_id}.sv.tsv 
    if bcftools view -H ${bcf} | grep "BND" > /dev/null; then
        bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\t%INFO/SR\t%INFO/CT[\t%RC][\t%RCL][\t%RCR]\n" ${bcf} | grep "BND" >> ${patient_id}_${sample_id}.sv.tsv
    else
        echo "No BND variants found, skipping BND query step."
    fi
    cat ${patient_id}_${sample_id}.sv.tsv | cut -f 1,2,5,6,7,8,19,30 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Left\\t"\$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"\$8;}' > ${patient_id}_${sample_id}.input.tsv
    cat ${patient_id}_${sample_id}.sv.tsv | cut -f 3,4,5,6,7,8,19,30 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Right\\t"\$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"\$8;}' >> ${patient_id}_${sample_id}.input.tsv
    """
}

// Feature annotation
process ALFRED {
    publishDir "${params.pubdir}/${patient_id}/alfred"
    container "trausch/alfred:latest"
    tag "ALFRED on ${patient_id}_${sample_id}"
    debug true
    errorStrategy "ignore"
    label "xxs_mem"
    
    input:
    tuple val(patient_id), val(sample_id), path("${patient_id}_${sample_id}.sv.tsv"), path("${patient_id}_${sample_id}.input.tsv")

    output:
    tuple val(patient_id), val(sample_id), path("${patient_id}_${sample_id}.sv.gene.tsv"), path("${patient_id}_${sample_id}.input.tsv")
    
    script:
    """
    alfred annotate -d 3000 -g ${params.gtf} -o ${patient_id}_${sample_id}.sv.gene.tsv ${patient_id}_${sample_id}.input.tsv
	rm ${patient_id}_${sample_id}.sv.tsv gene.bed
    """
}

// Modify final tables
process MERGE_TSV {
    publishDir "${params.pubdir2}/${patient_id}/delly_hg38/", mode: 'copy'
    publishDir "${params.pubdir}/${patient_id}/mergeTsv"
    container "quay.io/biocontainers/mulled-v2-76a0b41e09773ed659596514b804bc832021772e:d50820554e481bef847e84d0590124e985594f5d-0"
    tag "MERGE_TSV on ${patient_id}_${sample_id}"
    debug true
    errorStrategy "ignore"
    label "xxs_mem"
    
    input:
    tuple val(patient_id), val(sample_id), path("${patient_id}_${sample_id}.sv.gene.tsv"), path("${patient_id}_${sample_id}.input.tsv")

    output:
    path("${patient_id}_${sample_id}.combined.sv.tsv")

    script:
    """
    python ${params.script} --input_file1 ${patient_id}_${sample_id}.sv.gene.tsv --input_file2 ${patient_id}_${sample_id}.input.tsv --intermediate_file ${patient_id}_${sample_id}.intermediate.tsv --output_file ${patient_id}_${sample_id}.combined.sv.tsv
    """
}

// CNV calling & read-depth profiling
process DELLY_CNV {
    publishDir("${params.pubdir}/${patient_id}/CNV", mode: 'copy')
    container "dellytools/delly:latest"
    tag "DELLY_CNV on ${patient_id}_${meta.sample_id}"
    debug true
    label "s_mem"
    
    input:
    tuple val(patient_id), val(sample_id), val(meta), path(bam), path(bai)

    output:
    tuple val(patient_id), val(sample_id), val(meta), path("${meta.sample_id}.cnv.bcf"), path("${meta.sample_id}.cov.gz")

    script:
    """
    echo DELLY_CNV ${patient_id}_${meta.sample_id}
    delly cnv -a -g ${params.ref} -m ${params.map} -c ${meta.sample_id}.cov.gz -o ${meta.sample_id}.cnv.bcf $bam
    """
}

// Save files with read-depth coverages
process CNV_UNZIP {
    publishDir "${params.pubdir2}/${patient_id}/delly_hg38/CNV_cov/", mode: 'copy', pattern: "*.cov"
    container "ubuntu:22.04"
    tag "CNV_UNZIP on ${patient_id}_${sample_id}"
    debug true
    label "xxs_mem"
    
    input:
    tuple val(patient_id), val(sample_id), val(meta), path(bcf), path(cov)

    output:
    tuple val(patient_id), val(sample_id), val(meta), path("${sample_id}.cnv.bcf"), path("${sample_id}.cov.gz"), path("${sample_id}.cov")

    script:
    """
    gunzip -c ${sample_id}.cov.gz > ${sample_id}.cov
    """
}

// Generation of plots using CNV profiles & segmentation
process CNV_PROFILES {
    publishDir "${params.pubdir2}/${patient_id}/delly_hg38/CNV_plots/", mode: 'copy'
    publishDir("${params.pubdir}/${patient_id}/Rplots", mode: 'copy')
    container "patricie/my-r-dnacopy-image:latest"
    tag "CNV_PROFILES on ${patient_id}_${meta.sample_id}"
    debug true
    label "xs_mem"
    
    input:
    tuple val(patient_id), val(sample_id), val(meta), path(bcf), path(cov)

    output:
    path("*")

    script:
    """
    Rscript ${params.Rscript} $cov ${meta.sample_id}
    """
}

// Optional
process CNV_CALLS {
    // publishDir("${params.pubdir}/${patient_id}/cnv_calls", mode: 'copy')
    container "staphb/bcftools:latest"
    tag "CNV_CALLS on ${patient_id}_${meta.sample_id}"
    label "xxs_mem"
    
    input:
    tuple val(patient_id), val(meta), path(bcf), path(cov)

    output:
    tuple val(patient_id), val(meta), path("${meta.sample_id}.bed"), path(cov)

    script:
    """
    bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" ${meta.sample_id}.cnv.bcf > ${meta.sample_id}.bed
    """
}

// Optional - visualize the CNV calls
process CNV_CALL_PLT {
    // publishDir("${params.pubdir}/${patient_id}/cnv_calls", mode: 'copy')
    container "patricie/my-r-dnacopy-image:latest"
    tag "CNV_CALL_PLT on ${patient_id}_${meta.sample_id}"
    debug true
    label "xs_mem"
    
    input:
    tuple val(patient_id), val(meta), path(bed), path(cov)

    output:
    path("*")

    script:
    """
    Rscript ${params.Rscript} $cov ${meta.sample_id} $bed
    """
}

// Calculate coverage, average read_lenght, SV
process COUNT_STATS{
    publishDir("${params.pubdir}/${patient_id}/stats", mode: 'copy')
    container "patricie/python-sam-bcf-image:latest"
    tag "COUNT_STATS on ${patient_id}_${meta.sample_id}"
    debug true
    label "s_mem"
    
    input:
    tuple val(patient_id), val(sample_id), val(meta), path(bam), path(bcf), path(csi)

    output:
    path("aggregate_metrics.txt")

    script:
    """
    samtools stats ${bam} > ${patient_id}_stats.txt
    grep '^COV' ${patient_id}_${sample_id}_stats.txt | cut -f 2- >> ${patient_id}_${sample_id}_stats.txt
    python ${params.script2} --input_file ${patient_id}_${sample_id}_stats.txt
    echo -n "Average length: " >> ${patient_id}_${sample_id}_stats.txt
    grep '^SN	average length:' ${patient_id}_${sample_id}_stats.txt | cut -f 2- >> ${patient_id}_${sample_id}_stats.txt
    echo -n "SV count: " >> ${patient_id}_${sample_id}_stats.txt
    bcftools view -H ${bcf} | wc -l >> ${patient_id}_${sample_id}_stats.txt
    cat ${patient_id}_${sample_id}_stats.txt >> aggregate_metrics.txt
    """

}
