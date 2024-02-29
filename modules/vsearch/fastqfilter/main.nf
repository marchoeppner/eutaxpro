process VSEARCH_FASTQFILTER {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${meta.sample_id}/VSEARCH", mode: 'copy'

    label 'short_serial'

    conda 'bioconda::vsearch=2.27.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.27.0--h6a68c12_0' :
        'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0' }"

    input:
    tuple val(meta), path(fq)

    output:
    tuple val(meta), path(filtered), emit: fasta
    path("versions.yml"), emit: versions

    script:
    filtered = fq.getBaseName() + '.filtered.fasta'

    """
    vsearch -fastq_filter $fq -fastq_maxee 0.5 -relabel Filtered -fastaout $filtered

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1  | sed -e "s/vsearch //g" -e "s/,.*//")
    END_VERSIONS
    """
}
