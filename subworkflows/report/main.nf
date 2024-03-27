include { JSON2TSV }           from './../../modules/helper/json2tsv'
include { JSON2MQC }           from './../../modules/helper/json2mqc'

ch_qc = Channel.from([])
ch_versions = Channel.from([])

workflow REPORT {

    take:
    json

    main:

    JSON2TSV(
        json
    )

    JSON2MQC(
        json
    )

    emit:
    tsv = JSON2TSV.out.tsv
    mqc_json = JSON2MQC.out.mqc_json
    versions = ch_versions
    qc = ch_qc
}