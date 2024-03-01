include { VSEARCH_FASTQMERGE }                          from './../../modules/vsearch/fastqmerge'
include { VSEARCH_DEREPFULL }                           from './../../modules/vsearch/derep'
include { VSEARCH_DEREPFULL as VSEARCH_DEREPFULL_ALL }  from './../../modules/vsearch/derep'
include { VSEARCH_FASTQFILTER }                         from './../../modules/vsearch/fastqfilter'
include { VSEARCH_CLUSTER_SIZE }                        from './../../modules/vsearch/cluster_size'
include { VSEARCH_CLUSTER_SIZE as VSEARCH_CLUSTER_ALL } from './../../modules/vsearch/cluster_size'
include { VSEARCH_CLUSTER_UNOISE }                      from './../../modules/vsearch/unoise'
include { VSEARCH_UCHIME3_DENOVO }                      from './../../modules/vsearch/uchime3/denovo'
include { VSEARCH_USEARCH_GLOBAL }                      from './../../modules/vsearch/usearch_global'

include { VSEARCH_SINTAX }                              from './../../modules/vsearch/sintax'
include { VSEARCH_EXTRACT_NONCHIMERIC }                 from './../../modules/vsearch/extract_nonchimeric'
include { VSEARCH_EXTRACT_NONCHIMERIC as VSEARCH_EXTRACT_NONCHIMERIC_PER_SAMPLE } from './../../modules/vsearch/extract_nonchimeric'

ch_versions = Channel.from([])
ch_reports  = Channel.from([])
ch_qc_files = Channel.from([])

workflow VSEARCH_WORKFLOW {
    take:
    reads
    sintax_db

    main:

    // Merge PE files
    VSEARCH_FASTQMERGE(
        reads.map { m, r -> [m, r[0], r[1]] }
    )
    ch_versions = ch_versions.mix(VSEARCH_FASTQMERGE.out.versions)

    // Files merged reads using static parameters
    // This is not ideal and could be improved!
    VSEARCH_FASTQFILTER(
        VSEARCH_FASTQMERGE.out.fastq
    )
    ch_versions = ch_versions.mix(VSEARCH_FASTQFILTER.out.versions)

    VSEARCH_FASTQFILTER.out.fasta.map { m, f -> f }.collectFile(name: 'all.fasta').map { fasta ->
        [ [sample_id: 'all' ], fasta ]
    }.set { all_seqs }

    // Dereplicate the concatenated set of sequences
    VSEARCH_DEREPFULL_ALL(
        all_seqs
    )

    // The initial set of ASUs
    VSEARCH_CLUSTER_UNOISE(
        VSEARCH_DEREPFULL_ALL.out.fasta
    )

    // the denoised and chimera-filtered ASUs
    VSEARCH_UCHIME3_DENOVO(
        VSEARCH_CLUSTER_UNOISE.out.fasta
    )

    // We now make OTUs because that is good enough for our purpose
    VSEARCH_CLUSTER_ALL(
        VSEARCH_UCHIME3_DENOVO.out.fasta,
        params.vsearch_cluster_id

    )
    // We taxonomically map the OTUS
    VSEARCH_SINTAX(
        VSEARCH_CLUSTER_ALL.out.fasta,
        sintax_db
    )

    
    emit:
    tsv = VSEARCH_SINTAX.out.tsv
    versions = ch_versions
    fasta = VSEARCH_CLUSTER_ALL.out.fasta
    qc = ch_qc_files
}
