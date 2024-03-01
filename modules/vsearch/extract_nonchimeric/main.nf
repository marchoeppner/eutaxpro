process VSEARCH_EXTRACT_NONCHIMERIC {
    tag "${meta.sample_id}"

    label 'short_serial'

    conda 'sed=4.7'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta1)
    path(uc)
    path(fasta2)
    val(mincov)

    output:
    
    tuple val(meta), path(nonchim), emit: fasta

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: fasta1.getBaseName()

    nonchim = prefix + '.nonchimeric.fasta'

    """
    vsearch_extract_nonchimeric.pl $fasta1 $uc $fasta2 $mincov > $nonchim

    """
}
