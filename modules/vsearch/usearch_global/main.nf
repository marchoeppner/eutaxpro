process VSEARCH_USEARCH_GLOBAL {
    
    tag "${meta.sample_id}"

    label 'short_serial'

    conda 'bioconda::vsearch=2.27.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.27.0--h6a68c12_0' :
        'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0' }"

    input:
    tuple val(meta),path(fasta)
    path(db)

    output:
    tuple val(meta), path(tabbed), emit: fasta
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    tabbed = prefix + '.usearch_global.tsv'

    """
    vsearch --usearch_global $fasta \
    -db $db \
    -otutabout $tabbed $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1 | sed -e "s/vsearch //g" -e "s/,.*//")
    END_VERSIONS
    """
}
