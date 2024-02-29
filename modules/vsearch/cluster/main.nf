process VSEARCH_CLUSTER {
    
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${meta.sample_id}/VSEARCH", mode: 'copy'

    label 'short_serial'

    conda 'bioconda::vsearch=2.27.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.27.0--h6a68c12_0' :
        'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0' }"

    input:
    tuple val(meta), path(fa)
    val(cluster_id)

    output:
    tuple val(meta), path(cluster), emit: fasta
    tuple val(meta),path(uc), emit: uc
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix?: "${meta.sample_id}"
    cluster = meta.sample_id + '.precluster.fasta'
    uc = meta.sample_id + '.precluster.uc'

    """
    vsearch --cluster_size $fa \
    --threads ${task.cpus} \
    --id $cluster_id \
    --uc $uc \
    --centroids $cluster $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1 | sed -e "s/vsearch //g" -e "s/,.*//")
    END_VERSIONS
    """
}
