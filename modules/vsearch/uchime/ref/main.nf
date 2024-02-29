process VSEARCH_DEREPFULL {
    
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${meta.sample_id}/VSEARCH", mode: 'copy'

    label 'short_serial'

    conda 'bioconda::vsearch=2.27.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.27.0--h6a68c12_0' :
        'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0' }"

    input:
    tuple val(meta), path(fa)

    output:
    tuple val(meta), path(derep), emit: fasta
    path("versions.yml"), emit: versions

    script:
    derep = meta.sample_id + '.derep.fasta'
    derep_uc = meta.sample_id + '.uc'
    options = ""
    if (!meta.sample_id == "all") {
        options = "--relabel ${meta.sample_id}_Derep"
    }
    """
    vsearch --derep_fulllength $fa \
    --strand plus \
    --sizeout \
    --uc $derep_uc \
    --output $derep $options

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1 | sed -e "s/vsearch //g" -e "s/,.*//")
    END_VERSIONS
    """
}