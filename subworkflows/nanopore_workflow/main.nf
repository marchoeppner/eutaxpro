/*
include modules
*/
include { PORECHOP_ABI }            from './../../modules/porechop/abi'
include { CHOPPER }                 from './../../modules/chopper'
include { VSEARCH_CLUSTER_SIZE }    from './../../modules/vsearch/cluster_size'
include { VSEARCH_ORIENT }          from './../../modules/vsearch/orient'
include { VSEARCH_USEARCH_GLOBAL as VSEARCH_USEARCH_GLOBAL_ONT }  from './../../modules/vsearch/usearch_global'
include { VSEARCH_SINTAX }          from './../../modules/vsearch/sintax'
include { SINTAX_OTU2TAB }          from './../../modules/helper/sintax_otu2tab'
include { SINTAX_OTU2JSON }         from './../../modules/helper/sintax_otu2json'
include { SAMTOOLS_CONSENSUS }      from './../../modules/samtools/consensus'
include { MINIMAP2 }                from './../../modules/minimap2/align'

/*
Include sub-workflows
*/
include { REMOVE_PCR_PRIMERS }      from './../remoce_pcr_primers'

ch_versions = Channel.from([])

workflow NANOPORE_WORKFLOW {

    take:
    reads
    sintax_db

    main:

    PORECHOP_ABI(
        reads
    )
    ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions)

    // Chop reads and remove reads outside the expected size range
    CHOPPER(
        PORECHOP_ABI.out.reads,
        params.ont_min_length,
        params.ont_max_length
    )
    ch_versions = ch_versions.mix(CHOPPER.out.versions)

    REMOVE_PCR_PRIMERS(
        CHOPPER.out.reads
    )

    // Dummy call to rename fastq entries for downstream quantification
    VSEARCH_ORIENT(
        REMOVE_PCR_PRIMERS.out.reads,
        sintax_db
    )

    // Make per-sample OTU, because well...
    VSEARCH_CLUSTER_SIZE(
        NANOFILT.out.reads,
        params.ont_cluster_id
    )
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_SIZE.out.versions)

    MINIMAP2(
        VSEARCH_CLUSTER_SIZE.out.fasta,
        reads
    )

    SAMTOOLS_CONSENSUS(
        MINIMAP2.out.bam
    )

   // We taxonomically map the OTUS
    VSEARCH_SINTAX(
        VSEARCH_CLUSTER_SIZE.out.fasta,
        sintax_db
    )
    ch_versions = ch_versions.mix(VSEARCH_SINTAX.out.versions)

    // We generate the OTU Table with sample IDs
    VSEARCH_USEARCH_GLOBAL_ONT(
        VSEARCH_CLUSTER_SIZE.out.fasta,
        VSEARCH_ORIENT.out.reads.map{ m,f -> f }
    )
    ch_versions = ch_versions.mix(VSEARCH_USEARCH_GLOBAL.out.versions)

    SINTAX_OTU2TAB(
        VSEARCH_SINTAX.out.tsv.join(VSEARCH_USEARCH_GLOBAL.out.tab)
    )

    SINTAX_OTU2JSON(
        VSEARCH_SINTAX.out.tsv.join(VSEARCH_USEARCH_GLOBAL.out.tab)
    )

    emit:
    fasta = VSEARCH_CLUSTER_SIZE.out.fasta
    versions = ch_versions

}