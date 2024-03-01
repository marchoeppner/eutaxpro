
include { INPUT_CHECK }                 from './../modules/input_check'
include { MULTIQC }                     from './../modules/multiqc/main'
include { VSEARCH_WORKFLOW }            from './../subworkflows/vsearch'
include { DADA2_WORKFLOW }              from './../subworkflows/dada2/dada2'
include { PTRIMMER }                    from './../modules/ptrimmer'
include { FASTP }                       from './../modules/fastp'
include { CAT_FASTQ }                   from './../modules/cat_fastq'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'

samplesheet             = params.input      ? Channel.fromPath(params.input)                         : Channel.value([])

// Each primer set has a specific ptrimmer config file
ptrimmer_config         = file(params.references.primers[params.primer_set].ptrimmer_config)

gene                    = params.references.primers[params.primer_set].gene

ch_db_dada2             = Channel.fromPath(params.references.genes[gene].dada2_db).collect()
ch_db_sintax            = Channel.fromPath(params.references.genes[gene].sintax_db).collect()

single_end              = params.single_end

trunc_def               = 0
trunclenf               = params.trunclenf ?: trunc_def
trunclenr               = params.trunclenr ?: trunc_def

tools = params.tools ? params.tools.split(',').collect { tool -> clean_tool(tool) } : []

if (('dada2' in tools) && !single_end && (params.trunclenf == null || params.trunclenr == null)) {
    find_truncation_values = true
    log.warn 'No DADA2 cutoffs were specified (`--trunclenf` & `--trunclenr`), will try to truncate automatically...'
}

ch_versions             = Channel.from([])
multiqc_files           = Channel.from([])

ch_reads_for_vsearch    = Channel.from([])
ch_reads_for_dada2      = Channel.from([])

ch_assembled_fasta      = Channel.from([])

workflow EUTAXPRO {
    main:

    INPUT_CHECK(samplesheet)

    // Branch input reads by sequencing technology
    INPUT_CHECK.out.reads.branch { m, r ->
        illumina: m.platform == 'ILLUMINA'
        nanopore: m.platform == 'NANOPORE'
        pacbio: m.platform == 'PACBIO'
    }.set { ch_reads_by_platform }
    // channel: [[ sample_id: xxx, platform: xxx ], [ reads ] ]

    FASTP(
        ch_reads_by_platform.illumina
    )

    // Split trimmed reads by sample to find multi-lane data set
    FASTP.out.reads.groupTuple().branch { meta, reads ->
        single: reads.size() == 1
            return [ meta, reads.flatten()]
        multi: reads.size() > 1
            return [ meta, reads.flatten()]
    }.set { ch_reads_illumina }

    // Concatenate samples with multiple PE files
    CAT_FASTQ(
        ch_reads_illumina.multi
    )

    ch_illumina_trimmed = ch_reads_illumina.single.mix(CAT_FASTQ.out.reads)

    ch_reads_for_vsearch = ch_reads_for_vsearch.mix(ch_illumina_trimmed)
    ch_reads_for_dada2 = ch_reads_for_dada2.mix(ch_illumina_trimmed)

    if ('vsearch' in tools) {
        // Remove PCR primers the right way
        // Dada2 wants hard-clipped reads instead...
        PTRIMMER(
            ch_reads_for_vsearch,
            ptrimmer_config
        )
        ch_versions = ch_versions.mix(PTRIMMER.out.versions)

        VSEARCH_WORKFLOW(
            PTRIMMER.out.reads,
            ch_db_sintax

        )
        multiqc_files = multiqc_files.mix(VSEARCH_WORKFLOW.out.qc)
        ch_versions = ch_versions.mix(VSEARCH_WORKFLOW.out.versions)
        ch_assembled_fasta = ch_assembled_fasta.mix(VSEARCH_WORKFLOW.out.fasta)
    }

    if ('dada2' in tools) {
        DADA2_WORKFLOW(
            ch_reads_for_dada2,
            single_end,
            find_trancation_values,
            trunclenf,
            trunclenr
        )
        ch_assembled_fasta = ch_assembled_fasta.mix(DADA2_MERGE.out.fasta)
        ch_versions = ch_versions.mix(DADA2_WORKFLOW.out.versions)
    }

    // Create list of software packages used
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    multiqc_files = multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    MULTIQC(
        multiqc_files.collect()
    )

    emit:
    qc = MULTIQC.out.html
    }

def clean_tool(String tool) {
    return tool.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')
}
