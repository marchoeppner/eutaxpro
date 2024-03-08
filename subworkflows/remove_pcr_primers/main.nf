include { CUTADAPT }                    from './../../modules/cutadapt'
include { FASTX_REVERSE_COMPLEMENT }    from './../../modules/fastx_toolkit/fastx_reverse_complement'
include { PTRIMMER }                    from './../../modules/ptrimmer'

ch_versions = Channel.from([])

workflow REMOVE_PCR_PRIMERS {
    take:
    reads
    ch_ptrimmer_config
    ch_primers
    ch_primers_rc

    main:
    // Allow use of cutadapt if need be
    if (params.cutadapt) {
        if (params.cutadapt_trim_3p) {
            FASTX_REVERSE_COMPLEMENT(
                ch_primers
            )
            ch_primers_rc = FASTX_REVERSE_COMPLEMENT.out.fasta
            ch_versions = ch_versions.mix(FASTX_REVERSE_COMPLEMENT.out.versions)
        }
        CUTADAPT(
            reads,
            ch_primers,
            ch_primers_rc
        )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)
        ch_reads_no_primers = CUTADAPT.out.reads
    } else {
        PTRIMMER(
            reads,
            ch_ptrimmer_config
        )
        ch_versions = ch_versions.mix(PTRIMMER.out.versions)
        ch_reads_no_primers = PTRIMMER.out.reads
    }

    emit:
    reads = ch_reads_no_primers
    versions = ch_versions
}
