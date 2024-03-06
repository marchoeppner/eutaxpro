
include { INPUT_CHECK }                 from './../modules/input_check'
include { MULTIQC }                     from './../modules/multiqc/main'
include { VSEARCH_WORKFLOW }            from './../subworkflows/vsearch'
include { DADA2_WORKFLOW }              from './../subworkflows/dada2/dada2'
include { FASTP }                       from './../modules/fastp'
include { CUTADAPT }                    from './../modules/cutadapt'
include { CAT_FASTQ }                   from './../modules/cat_fastq'
include { PORECHOP_PORECHOP }           from './../modules/porechop/porechop'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'

// The input sample sheet
samplesheet             = params.input ? Channel.fromPath(file(params.input, checkIfExists:true)) : Channel.value([])

// Primer sequences are either pre-configured or can be supplied by user
if (params.primer_set) {
    ch_primers              = Channel.fromPath(file(params.references.primers[params.primer_set].fasta, checkIfExits: true)).collect()
} else if (params.primers) {
    ch_primers              = Channel.fromPath(file(params.primers, checkIfExists: true)).collect()
} else {
    ch_primers              = Channel.from([])
} 

// The gene configured for the primer set
gene                    ? params.gene : params.references.primers[params.primer_set].gene

// The taxonomy databases for this gene
if (params.reference_base) {
    ch_db_sintax            = Channel.fromPath(params.references.genes[gene].sintax_db).collect()
} else {
    ch_db_sintax            = Channel.fromPath(file(params.references.genes[gene].sintax_url)).collect()
}

// PLACE HOLDER: data is single-end?
single_end              = params.single_end

// List of tool chains to run
tools = params.tools ? params.tools.split(',').collect { tool -> clean_tool(tool) } : []

// DADA2 options to determine need for automatic truncation
trunc_def               = 0
trunclenf               = params.trunclenf ?: trunc_def
trunclenr               = params.trunclenr ?: trunc_def
find_truncation_values  = false

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

    // trim nanopore reads
    PORECHOP_PORECHOP(
        ch_reads_by_platform.nanopore
    )
    ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions)
    multiqc_files = multiqc_files.mix(PORECHOP_PORECHOP.out.log)
    ch_reads_for_vsearch = ch_reads_for_vsearch.mix(PORECHOP_PORECHOP.out.reads)

    // trim illumina reads
    FASTP(
        ch_reads_by_platform.illumina
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)
    multiqc_files = multiqc_files.mix(FASTP.out.json)

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

    // Trim PCR primers with CUTADAPT
    CUTADAPT(
        ch_illumina_trimmed,
        ch_primers
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    ch_reads_for_vsearch = ch_reads_for_vsearch.mix(CUTADAPT.out.reads)
    ch_reads_for_dada2 = ch_reads_for_dada2.mix(CUTADAPT.out.reads)

    if ('vsearch' in tools) {
        
        VSEARCH_WORKFLOW(
            ch_reads_for_vsearch,
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
            find_truncation_values,
            trunclenf,
            trunclenr
        )
        ch_assembled_fasta = ch_assembled_fasta.mix(DADA2_WORKFLOW.out.fasta)
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
