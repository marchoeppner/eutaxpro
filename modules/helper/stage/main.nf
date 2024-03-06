process STAGE {
    publishDir "${params.outdir}/${meta.target}/${meta.tool}", mode: 'copy'

    tag "${meta.target}|${meta.tool}"

    input:
    tuple val(meta), path(thisFile)

    output:
    tuple val(meta), path(thisFile), emit: staged

    script:

    '''

    '''
}
