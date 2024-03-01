
include { DADA2_PREPROCESSING }         from './../preprocessing'
include { DADA2_ERROR }                 from './../../../modules/dada2/error'
include { DADA2_STATS }                 from './../../../modules/dada2/stats'
include { DADA2_MERGE }                 from './../../../modules/dada2/merge'
include { DADA2_DENOISING }             from './../../../modules/dada2/denoising'
include { DADA2_RMCHIMERA }             from './../../../modules/dada2/rmchimera'

ch_versions = Channel.from([])

workflow DADA2_WORKFLOW {
    take:
    reads
    single_end
    find_trancation_values
    trunclenf
    trunclenr

    main:

    // Clean up the reads so DADA2 can work with them
    DADA2_PREPROCESSING(
        reads,
        single_end,
        find_truncation_values,
        trunclenf,
        trunclenr
    )
    ch_filt_reads   = DADA2_PREPROCESSING.out.reads
    ch_versions     = ch_versions.mix(DADA2_PREPROCESSING.out.versions)

    // Model read errors per sample
    DADA2_ERROR(
        ch_filt_reads
    )
    ch_versions     = ch_versions.mix(DADA2_ERROR.out.versions)
    ch_errormodel   = DADA2_ERROR.out.errormodel

    ch_filt_reads
        .join(ch_errormodel)
        .set { ch_derep_errormodel }

    // Perform denoising of reads
    DADA2_DENOISING(ch_derep_errormodel.dump(tag: 'into_denoising'))
    ch_versions = ch_versions.mix(DADA2_DENOISING.out.versions.first())

    // Remove chimeras, if any
    DADA2_RMCHIMERA(DADA2_DENOISING.out.seqtab)

    //group by sequencing run & group by meta
    DADA2_PREPROCESSING.out.logs
        .join(DADA2_DENOISING.out.denoised)
        .join(DADA2_DENOISING.out.mergers)
        .join(DADA2_RMCHIMERA.out.rds)
    .set { ch_track_numbers }
    DADA2_STATS(ch_track_numbers)

    //merge if several runs, otherwise just publish
    DADA2_MERGE(
        DADA2_STATS.out.stats.map { meta, stats -> stats }.collect(),
        DADA2_RMCHIMERA.out.rds.map { meta, rds -> rds }.collect()
    )

    emit:
    fasta = DADA2_MERGE.out.fasta
    versions = ch_versions
}
