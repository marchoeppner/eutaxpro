include { VSEARCH_FASTQMERGE }                          from './../../modules/vsearch/fastqmerge'
include { VSEARCH_DEREPFULL }                           from './../../modules/vsearch/derep'
include { VSEARCH_DEREPFULL as VSEARCH_DEREPFULL_ALL }  from './../../modules/vsearch/derep'
include { VSEARCH_FASTQFILTER }                         from './../../modules/vsearch/fastqfilter'
include { VSEARCH_CLUSTER_SIZE }                        from './../../modules/vsearch/cluster_size'
include { VSEARCH_CLUSTER_UNOISE }                      from './../../modules/vsearch/unoise'
include { VSEARCH_UCHIME3_DENOVO }                      from './../../modules/vsearch/uchime3/denovo'
include { VSEARCH_USEARCH_GLOBAL }                      from './../../modules/vsearch/usearch_global'
include { VSEARCH_SINTAX }                              from './../../modules/vsearch/sintax'
include { SINTAX_OTU2TAB }                              from './../../modules/helper/sintax_otu2tab'

ch_versions = Channel.from([])
ch_reports  = Channel.from([])
ch_qc_files = Channel.from([])

workflow VSEARCH_WORKFLOW {
    take:
    reads
    sintax_db

    main:

    reads.branch { m, r ->
        paired: !m.single_end
        unpaired: m.single_end
    }.set { ch_trimmed_reads }

    // Merge PE files and attach sample names
    VSEARCH_FASTQMERGE(
        ch_trimmed_reads.paired.map { m, r -> [m, r[0], r[1]] }
    )
    ch_versions = ch_versions.mix(VSEARCH_FASTQMERGE.out.versions)

    // paired and unpaired reads after optional merging and read name tagging
    // we now have [ meta, fastq ]
    ch_merged_reads = VSEARCH_FASTQMERGE.out.fastq

    // Files merged reads using static parameters
    // This is not ideal and could be improved!
    VSEARCH_FASTQFILTER(
        VSEARCH_FASTQMERGE.out.fastq
    )
    ch_versions = ch_versions.mix(VSEARCH_FASTQFILTER.out.versions)

    VSEARCH_FASTQFILTER.out.fasta.map { m, f -> f }.collectFile(name: 'all.fasta').map { fasta ->
        [ [sample_id: params.run_name ], fasta ]
    }.set { all_seqs }

    // Dereplicate the concatenated set of sequences
    VSEARCH_DEREPFULL_ALL(
        all_seqs
    )
    ch_versions = ch_versions.mix(VSEARCH_DEREPFULL_ALL.out.versions)

    // The initial set of ASUs
    VSEARCH_CLUSTER_UNOISE(
        VSEARCH_DEREPFULL_ALL.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_UNOISE.out.versions)

    // the denoised and chimera-filtered ASUs
    VSEARCH_UCHIME3_DENOVO(
        VSEARCH_CLUSTER_UNOISE.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_UCHIME3_DENOVO.out.versions)

    // We now make OTUs because that is good enough for our purpose
    VSEARCH_CLUSTER_SIZE(
        VSEARCH_UCHIME3_DENOVO.out.fasta,
        params.vsearch_cluster_id

    )
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_SIZE.out.versions)

    // We taxonomically map the OTUS
    VSEARCH_SINTAX(
        VSEARCH_CLUSTER_SIZE.out.fasta,
        sintax_db
    )
    ch_versions = ch_versions.mix(VSEARCH_SINTAX.out.versions)

    // We generate the OTU Table with sample IDs
    VSEARCH_USEARCH_GLOBAL(
        VSEARCH_CLUSTER_SIZE.out.fasta,
        ch_merged_reads.map { m, f -> f }.collectFile(name: 'all.merged.fastq')
    )
    ch_versions = ch_versions.mix(VSEARCH_USEARCH_GLOBAL.out.versions)

    SINTAX_OTU2TAB(
        VSEARCH_SINTAX.out.tsv.join(VSEARCH_USEARCH_GLOBAL.out.tab)
    )

    emit:
    tsv = SINTAX_OTU2TAB.out.tsv
    versions = ch_versions
    fasta = VSEARCH_CLUSTER_SIZE.out.fasta
    qc = ch_qc_files
}
