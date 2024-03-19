process CHOPPER {
    tag "$meta.sample_id"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.7.0--hdcf5f25_0' :
        'quay.io/biocontainers/chopper:0.7.0--hdcf5f25_0' }"

    input:
    tuple val(meta), path(reads)
    val(min_len)
    val(max_len)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path 'versions.yml'                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: reads.getBaseName() + '.chopped'
    
    """
    gunzip -c $reads | \\
    chopper \\
        --threads $task.cpus \\
        $args \\
        --minlength $min_len \\
        --maxlength $max_len \\
        -q 5 |
    gzip > ${prefix}.fastq.gz \\
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$( chopper --version | sed 's/chopper //' )
    END_VERSIONS
    """
}
