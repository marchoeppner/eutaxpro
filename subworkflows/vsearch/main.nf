include { VSEARCH_FASTQMERGE }                          from './../../modules/vsearch/fastqmerge'
include { VSEARCH_DEREPFULL }                           from './../../modules/vsearch/derep'
include { VSEARCH_DEREPFULL as VSEARCH_DEREPFULL_ALL}   from './../../modules/vsearch/derep'
include { VSEARCH_FASTQFILTER }                         from './../../modules/vsearch/fastqfilter'
include { VSEARCH_CLUSTER }                             from './../../modules/vsearch/cluster'
include { VSEARCH_CLUSTER as VSEARCH_CLUSTER_ALL }      from './../../modules/vsearch/cluster'
include { VSEARCH_UCHIME_DENOVO }                       from './../../modules/vsearch/uchime/denovo'
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
        reads.map {m,r -> [m, r[0],r[1]]}
    )
    ch_versions = ch_versions.mix(VSEARCH_FASTQMERGE.out.versions)

    // Files merged reads using static parameters
    // This is not ideal and could be improved!
    VSEARCH_FASTQFILTER(
        VSEARCH_FASTQMERGE.out.fastq
    )
    ch_versions = ch_versions.mix(VSEARCH_FASTQFILTER.out.versions)

    // Reduce reads into unique sequences
    VSEARCH_DEREPFULL(
        VSEARCH_FASTQFILTER.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_DEREPFULL.out.versions)

    VSEARCH_DEREPFULL.out.fasta.map {m,f -> f}.collectFile(name: 'all.fasta').map { fasta ->
        [ [sample_id: "all" ], fasta ]
    }.set { all_seqs }

    VSEARCH_DEREPFULL_ALL(
        all_seqs
    )

    VSEARCH_CLUSTER(
        VSEARCH_DEREPFULL_ALL.out.fasta,
        params.vsearch_cluster_id
    )

    VSEARCH_UCHIME_DENOVO(
        VSEARCH_CLUSTER.out.fasta
    )

    VSEARCH_EXTRACT_NONCHIMERIC(
        VSEARCH_DEREPFULL_ALL.out.fasta.collect(),
        VSEARCH_CLUSTER.out.uc.map {m,u -> u}.collect(),
        VSEARCH_UCHIME_DENOVO.out.fasta.map{ m,f -> f }.collect(),
        params.vsearch_min_cov
    )

    VSEARCH_EXTRACT_NONCHIMERIC_PER_SAMPLE(
        all_seqs.collect(),
        VSEARCH_DEREPFULL_ALL.out.uc.map{m,u ->u}.collect(),
        VSEARCH_EXTRACT_NONCHIMERIC.out.fasta.map{m,f -> f}.collect(),
        params.vsearch_min_cov
    )

    VSEARCH_SINTAX(
        VSEARCH_EXTRACT_NONCHIMERIC_PER_SAMPLE.out.fasta,
        sintax_db
    )

    emit:
    tsv = VSEARCH_SINTAX.out.tsv
    versions = ch_versions
    fasta = VSEARCH_EXTRACT_NONCHIMERIC_PER_SAMPLE.out.fasta
    qc = ch_qc_files
    
}
