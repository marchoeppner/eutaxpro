/*
Include Modules
*/
include { VSEARCH_FASTQMERGE }                          from './../../modules/vsearch/fastqmerge'
include { VSEARCH_DEREPFULL }                           from './../../modules/vsearch/derep'
include { VSEARCH_DEREPFULL as VSEARCH_DEREPFULL_ALL }  from './../../modules/vsearch/derep'
include { VSEARCH_FASTQFILTER }                         from './../../modules/vsearch/fastqfilter'
include { VSEARCH_CLUSTER_SIZE }                        from './../../modules/vsearch/cluster_size'
include { VSEARCH_CLUSTER_UNOISE }                      from './../../modules/vsearch/unoise'
include { VSEARCH_UCHIME3_DENOVO }                      from './../../modules/vsearch/uchime3/denovo'
include { VSEARCH_USEARCH_GLOBAL }                      from './../../modules/vsearch/usearch_global'
include { VSEARCH_SINTAX }                              from './../../modules/vsearch/sintax'
include { SINTAX_OTU2JSON }                             from './../../modules/helper/sintax_otu2json'

/* 
Set default channels
*/
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

    /*
    Merge PE files
    */
    VSEARCH_FASTQMERGE(
        ch_trimmed_reads.paired.map { m, r -> [m, r[0], r[1]] }
    )
    ch_versions = ch_versions.mix(VSEARCH_FASTQMERGE.out.versions)

    // paired and unpaired reads after optional merging and read name tagging
    // we now have [ meta, fastq ]
    ch_merged_reads = VSEARCH_FASTQMERGE.out.fastq

    /*
    Filter merged reads using static parameters
    This is not ideal and could be improved!
    */
    VSEARCH_FASTQFILTER(
        VSEARCH_FASTQMERGE.out.fastq
    )
    ch_versions = ch_versions.mix(VSEARCH_FASTQFILTER.out.versions)

    /*
    Dereplicate the individual samples
    */
    VSEARCH_DEREPFULL(
        VSEARCH_FASTQFILTER.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_DEREPFULL.out.versions)

    /*
    Denoise invidual samples
    */
    VSEARCH_CLUSTER_UNOISE(
        VSEARCH_DEREPFULL.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_UNOISE.out.versions)

    /*
    We concatenate all the filtered reads for joint de-replication
    */
    VSEARCH_CLUSTER_UNOISE.out.fasta.map { m, f -> f }.collectFile(name: 'all.fasta').map { fasta ->
        [ [sample_id: params.run_name ], fasta ]
    }.set { all_seqs }

    /*
    Dereplicate the concatenated set of sequences
    */
    VSEARCH_DEREPFULL_ALL(
        all_seqs
    )
    ch_versions = ch_versions.mix(VSEARCH_DEREPFULL_ALL.out.versions)

    /*
    Remove chimera from the set of ASUs
    */
    VSEARCH_UCHIME3_DENOVO(
        VSEARCH_DEREPFULL_ALL.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_UCHIME3_DENOVO.out.versions)

    /*
    Perform clustering of the chimera-filteredf ASUs into OTUs
    */
    VSEARCH_CLUSTER_SIZE(
        VSEARCH_UCHIME3_DENOVO.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_SIZE.out.versions)

    // We taxonomically map the OTUS
    VSEARCH_SINTAX(
        VSEARCH_CLUSTER_SIZE.out.fasta,
        sintax_db
    )
    ch_versions = ch_versions.mix(VSEARCH_SINTAX.out.versions)

    /*
    We generate the OTU Table with sample IDs
    To this end all reads get sample IDs attached and are jointly mapped
    and counted against our final set of OTUs
    */
    VSEARCH_USEARCH_GLOBAL(
        VSEARCH_CLUSTER_SIZE.out.fasta,
        ch_merged_reads.map { m, f -> f }.collectFile(name: 'all.merged.fastq')
    )
    ch_versions = ch_versions.mix(VSEARCH_USEARCH_GLOBAL.out.versions)

    /*
    Convert the sintax results into a count table in JSON format
    */
    SINTAX_OTU2JSON(
        VSEARCH_SINTAX.out.tsv.join(VSEARCH_USEARCH_GLOBAL.out.tab)
    )

    emit:
    json = SINTAX_OTU2JSON.out.json
    versions = ch_versions
    fasta = VSEARCH_CLUSTER_SIZE.out.fasta
    qc = ch_qc_files
}
