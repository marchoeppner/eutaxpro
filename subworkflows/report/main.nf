include { JSON2TSV }           from './../../modules/helper/json2tsv'

workflow REPORT {

    take:
    json

    main:

    JSON2TSV(
        json
    )

    emit:
    tsv = JSON2TSV.out.tsv

}