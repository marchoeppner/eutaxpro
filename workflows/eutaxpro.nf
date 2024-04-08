
/*
Import modules
*/
include { INPUT_CHECK }                 from './../modules/input_check'
include { MULTIQC }                     from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'

/*
Import sub workflows
*/
include { ILLUMINA_WORKFLOW }           from './../subworkflows/illumina_workflow'
include { REPORT }                      from './../subworkflows/report'

/*
Set default channels
*/
samplesheet             = params.input ? Channel.fromPath(file(params.input, checkIfExists:true)) : Channel.value([])

gene = null

/* 
Primer sets are either pre-configured or can be supplied by user,
preferably as Ptrimmer config, or as fasta for cutadapt.
*/
if (params.primer_set) {
    ch_ptrimmer_config      = Channel.fromPath(file(params.references.primers[params.primer_set].ptrimmer_config, checkIfExits: true)).collect()
    gene                    = params.references.primers[params.primer_set].gene
    ch_primers              = Channel.fromPath(file(params.references.primers[params.primer_set].fasta, checkIfExits: true)).collect()
    ch_primers_rc           = Channel.fromPath(file(params.references.primers[params.primer_set].fasta, checkIfExits: true)).collectFile(name: 'primers_rc.fasta')
} else if (params.primers_txt) {
    ch_ptrimmer_config      = Channel.fromPath(file(params.primers_txt, checkIfExists: true)).collect()
    gene                    = params.gene.toLowerCase()
    ch_primers              = Channel.from([])
    ch_primers_rc           = Channel.from([])
} else if (params.primers_fa) {
    ch_ptrimmer_config      = Channel.from([])
    ch_primers              = Channel.fromPath(file(params.primers_fa, checkIfExists: true)).collect()
    ch_primers_rc           = Channel.fromPath(file(params.primers_fa, checkIfExists: true)).collectFile(name: 'primers_rc.fasta')
    gene                    = params.gene.toLowerCase()
}

// The taxonomy database for this gene
if (params.reference_base && gene) {
    ch_db_sintax            = Channel.fromPath(params.references.genes[gene].sintax_db, checkIfExists: true).collect()
} else if (gene) {
    ch_db_sintax            = Channel.fromPath(file(params.references.genes[gene].sintax_url)).collect()
}

ch_versions             = Channel.from([])
multiqc_files           = Channel.from([])

workflow EUTAXPRO {
    main:

    INPUT_CHECK(samplesheet)

    // Branch input reads by sequencing technology
    INPUT_CHECK.out.reads.branch { m, r ->
        illumina: m.platform == 'ILLUMINA'
        torrent: m.platform = 'TORRENT'
        nanopore: m.platform == 'NANOPORE'
        pacbio: m.platform == 'PACBIO'
    }.set { ch_reads_by_platform }
    // channel: [[ sample_id: xxx, platform: xxx ], [ reads ] ]

    /*
    Processing of Illumina reads
    */
    ILLUMINA_WORKFLOW(
        ch_reads_by_platform.illumina,
        ch_ptrimmer_config,
        ch_primers,
        ch_primers_rc,
        ch_db_sintax
    )
    ch_versions = ch_versions.mix(ILLUMINA_WORKFLOW.out.versions)
    multiqc_files = multiqc_files.mix(ILLUMINA_WORKFLOW.out.qc)
    ch_json_report = ILLUMINA_WORKFLOW.out.json

    /*
    Create human-readable report(s)
    */
    REPORT(
        ch_json_report
    )
    multiqc_files = multiqc_files.mix(REPORT.out.mqc_json)

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
