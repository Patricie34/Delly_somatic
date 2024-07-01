workflow {
    input_bams = file("path/to/bam")
    tsv_ch = Channel.fromPath("samplesheet.tsv")
    | splitCsv( header: true )
    | map { row -> 
    def bam_path = file("${input_bams}/${row.sample_id}.bam")
    def bai_path = file("${input_bams}/${row.sample_id}.bai")
    [[sample_id: row.sample_id, name: row.name, type: row.type], [bam_path, bai_path]
    }
    | branch { meta, reads ->
        tumor: meta.type == "tumor"
        normal: meta.type == "normal"
    }
    | set { samples }

    samples.tumor.view()
    samples.normal.view()
}
