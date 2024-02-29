include { VSEARCH_WORKFLOW }            from './../vsearch'
include { DADA2_PREPROCESSING }         from './../dada2'
include { DADA2_ERROR }                 from './../../modules/dada2/error'
include { DADA2_STATS }                 from './../../modules/dada2/stats'
include { DADA2_MERGE }                 from './../../modules/dada2/merge'
include { DADA2_DENOISING }             from './../../modules/dada2/denoising'
include { DADA2_RMCHIMERA }             from './../../modules/dada2/rmchimera'
include { FASTP }                       from './../../modules/fastp'
include { PTRIMMER }                    from './../../modules/ptrimmer'
include { CAT_FASTQ }                   from './../../modules/cat_fastq'

ch_versions = Channel.from([])
multiqc_files = Channel.from([])
ch_assembled_fasta = Channel.from([])

workflow ILLUMINA_WORKFLOW {

    take:
    tools
    reads
    ptrimmer_config
    find_truncation_values
    trunclenf
    trunclenr

    main:
    // Trim Illumina reads
    FASTP(
        reads
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)
    multiqc_files = multiqc_files.mix(FASTP.out.json)

    // Split trimmed reads by sample to find multi-lane data set
    FASTP.out.reads.groupTuple().branch { meta,reads ->
        single: reads.size() == 1
            return [ meta, reads.flatten()]
        multi: reads.size() > 1
            return [ meta, reads.flatten()]
    }.set {ch_reads_illumina }

    // Concatenate samples with multiple PE files
    CAT_FASTQ(
        ch_reads_illumina.multi
    )

    // Merge singleton reads and concatenated reads
    ch_reads_illumina_merged = ch_reads_illumina.single.mix(CAT_FASTQ.out.reads)

    if ('vsearch' in tools) {

        // Remove PCR primers
        PTRIMMER(
            ch_reads_illumina_merged,
            ptrimmer_config
        )

        // Process reads with VSEARCH
        VSEARCH_WORKFLOW(
            PTRIMMER.out.reads.map { m,reads ->
                [ m, reads[0],reads[1] ]
            }
        )

        ch_assembled_fasta = ch_assembled_fasta.mix(VSEARCH_WORKFLOW.out.fasta)
    }

    if ('dada2' in tools ) {

        // PRocess reads with Dada2
        DADA2_PREPROCESSING(
            ch_reads_illumina_merged,
            params.single_end,
            find_truncation_values,
            trunclenf,
            trunclenr 
        )
        ch_filt_reads   = DADA2_PREPROCESSING.out.reads
        ch_versions     = ch_versions.mix(DADA2_PREPROCESSING.out.versions)

        DADA2_ERROR(
            ch_filt_reads
        )
        ch_versions     = ch_versions.mix(DADA2_ERROR.out.versions)
        ch_errormodel   = DADA2_ERROR.out.errormodel

        ch_filt_reads
            .join( ch_errormodel )
            .set { ch_derep_errormodel }
        DADA2_DENOISING ( ch_derep_errormodel.dump(tag: 'into_denoising')  )
        ch_versions = ch_versions.mix(DADA2_DENOISING.out.versions.first())

        DADA2_RMCHIMERA ( DADA2_DENOISING.out.seqtab )

        //group by sequencing run & group by meta
        DADA2_PREPROCESSING.out.logs
            .join( DADA2_DENOISING.out.denoised )
            .join( DADA2_DENOISING.out.mergers )
            .join( DADA2_RMCHIMERA.out.rds )
            .set { ch_track_numbers }
        DADA2_STATS ( ch_track_numbers )

        //merge if several runs, otherwise just publish
        DADA2_MERGE (
            DADA2_STATS.out.stats.map { meta, stats -> stats }.collect(),
            DADA2_RMCHIMERA.out.rds.map { meta, rds -> rds }.collect() 
        )
        ch_assembled_fasta = ch_assembled_fasta.mix(DADA2_MERGE.out.fasta)
    }

    emit:
    versions = ch_versions
    fasta = ch_assembled_fasta
    qc = multiqc_files

}