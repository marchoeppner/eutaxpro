/*
Include Modules
*/
include { VSEARCH_FASTQMERGE }                              from './../../modules/vsearch/fastqmerge'
include { VSEARCH_DEREPFULL }                               from './../../modules/vsearch/derep'
include { VSEARCH_SORTBYSIZE }                              from './../../modules/vsearch/sortbysize'
include { VSEARCH_FASTQFILTER }                             from './../../modules/vsearch/fastqfilter'
include { VSEARCH_CLUSTER_SIZE as VSEARCH_CLUSTER_SIZE_SINGLE } from './../../modules/vsearch/cluster_size'
include { VSEARCH_CLUSTER_UNOISE }                          from './../../modules/vsearch/unoise'
include { VSEARCH_UCHIME_DENOVO as VSEARCH_UCHIME_DENOVO_SINGLE } from './../../modules/vsearch/uchime/denovo'
include { VSEARCH_USEARCH_GLOBAL }                          from './../../modules/vsearch/usearch_global'
include { VSEARCH_SINTAX }                                  from './../../modules/vsearch/sintax'
include { SINTAX_OTU2JSON }                                 from './../../modules/helper/sintax_otu2json'

/* 
Set default channels
*/
ch_versions = Channel.from([])
ch_reports  = Channel.from([])
ch_qc_files = Channel.from([])

workflow VSEARCH_SINGLE {
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

    /*
    paired and unpaired reads after optional merging and read name tagging
    we now have [ meta, fastq ]
    */
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
    Dereplicate the filtered reads
    */
    VSEARCH_DEREPFULL(
        VSEARCH_FASTQFILTER.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_DEREPFULL.out.versions)

    /*
    Cluster dereplicated sequences
    */
    VSEARCH_CLUSTER_SIZE_SINGLE(
        VSEARCH_DEREPFULL.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_SIZE_SINGLE.out.versions)

    /*
    Detect chimeras denovo and remove from OTU set
    */
    VSEARCH_UCHIME_DENOVO_SINGLE(
        VSEARCH_CLUSTER_SIZE_SINGLE.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_UCHIME_DENOVO_SINGLE.out.versions)

    /*
    Map our reads back to our OTUs, including filtered reads
    to get an accurate count 
    */
    VSEARCH_USEARCH_GLOBAL(
        VSEARCH_UCHIME_DENOVO_SINGLE.out.fasta.join(
            VSEARCH_FASTQMERGE.out.fastq
        )
    )
    ch_versions = ch_versions.mix(VSEARCH_USEARCH_GLOBAL.out.versions)

    /*
    Determine best-guess taxon for our OTUs
    */
    VSEARCH_SINTAX(
        VSEARCH_UCHIME_DENOVO_SINGLE.out.fasta,
        sintax_db
    )
    ch_versions = ch_versions.mix(VSEARCH_SINTAX.out.versions)

    /*
    Convert the sintax results into a count table in JSON format
    */
    SINTAX_OTU2JSON(
        VSEARCH_SINTAX.out.tsv.join(
            VSEARCH_USEARCH_GLOBAL.out.tab
        )
    )

    emit:
    json = SINTAX_OTU2JSON.out.json
    versions = ch_versions
    fasta = VSEARCH_CLUSTER_SIZE_SINGLE.out.fasta
    qc = ch_qc_files

}