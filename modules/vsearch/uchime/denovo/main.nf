process VSEARCH_UCHIME_DENOVO {
    tag "${meta.sample_id}"

    label 'short_serial'

    conda 'bioconda::vsearch=2.27.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.27.0--h6a68c12_0' :
        'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0' }"

    input:
    tuple val(meta), path(fa)

    output:
    tuple val(meta), path(nonchimera), emit: fasta
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    
    nonchimera = prefix + '.uchime_denovo.fasta'
    derep_uc = prefix + '.uchime_denovo.uc'

    """
    vsearch --uchime_denovo $fa \
    --sizein \
    --sizeout \
    --nonchimera $nonchimera

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1 | sed -e "s/vsearch //g" -e "s/,.*//")
    END_VERSIONS
    """
}
