include { GUNZIP } from './../modules/gunzip'
include { UNZIP } from './../modules/unzip'

genes = params.references.genes.keySet()

midori_files = []

// For all genes of interest, recover supported tools and the corresponding database link
genes.each { gene ->
    tools = params.references.genes[gene].keySet()
    tools.each { tool ->
        meta = [:]
        meta.target = gene
        meta.tool = tool
        midori_files << [ meta ,file(params.references.genes[gene][tool].url, checkIfExists: true) ]
    }
}

ch_files = Channel.fromList(midori_files)

workflow BUILD_REFERENCES {

    main:

    ch_files.branch {
        zipped: it[1].toString().contains(".zip")
        gzipped: it[1].toString().contains(".gz")
        uncompressed: !it[1].toString().contains(".zip") && !it[1].toString().contains(".gz")
    }.set { ch_branched_files }

    UNZIP(
        ch_branched_files.zipped
    )
    
    GUNZIP(
        ch_branched_files.gzipped
    )

}