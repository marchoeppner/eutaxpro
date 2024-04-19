include { JSON2TSV }           from './../../modules/helper/json2tsv'

ch_qc = Channel.from([])
ch_versions = Channel.from([])
ch_mqc = Channel.from([])
ch_tsv = Channel.from([])

workflow REPORT {

    take:
    json

    main:

    JSON2TSV(
        json
    )
    ch_tsv = JSON2TSV.out.tsv

    emit:
    tsv = ch_tsv
    mqc_json = ch_mqc
    versions = ch_versions
    qc = ch_qc
}