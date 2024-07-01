nextflow.enable.dsl=2

// somatic SV calling

process delly_call {

    publishDir "${params.pubdir}/results/SV_delly/call", mode: 'copy'
    container "dellytools/delly:latest"

    input:
    tuple val(sample), path(t_bam), path(t_bai), path(n_bam), path(n_bai)
    
    output:
    path '*'

    shell:
    """
    delly call -o ${sample}_SV_delly.bcf -g ${params.ref} ${t_bam} ${n_bam}
    """

}

// somatic pre-filtering

process delly_pre_filter {

    publishDir "${params.pubdir}/results/SV_delly/filter", mode: 'copy'
    container "dellytools/delly:latest"

    input:
    tuple val(sample), path(bcf), path(bcf_csi)
    
    output:
    tuple val(sample), path("${sample}_SV_delly.pre.bcf")

    shell:
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

    shell:
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

    shell:
    """
    delly filter -f somatic -o ${sample}_SV_delly.somatic.bcf -s samples.tsv ${geno.bcf}
    """

}



workflow {
    // create channel
    tsv_ch = Channel.fromPath("Delly_somatic/samplesheet.tsv")
    | splitCsv( header: true )
    | map { row -> [[id: row.id, repeat: row.repeat, type: row.type], [file(row.bam), file(row.bai)]] }
    | branch { meta, reads ->
        tumor: meta.type == "tumor"
        normal: meta.type == "normal"
    }
    | set { samples }

    samples.tumor.view()
    samples.normal.view()
}



