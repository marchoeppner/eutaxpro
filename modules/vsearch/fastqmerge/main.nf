process VSEARCH_FASTQMERGE {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${meta.sample_id}/VSEARCH", mode: 'copy'

    label 'short_serial'

    conda 'bioconda::vsearch=2.27.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.27.0--h6a68c12_0' :
        'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0' }"

    input:
    tuple val(meta), path(fwd), path(rev)

    output:
    tuple val(meta), path(merged), emit: fastq
    path("versions.yml"), emit: versions

    script:
    merged = meta.sample_id + '.merged.fastq'

    """
    vsearch --fastq_merge $fwd --reverse $rev \
    --fastqout $merged \
    --fastq_eeout

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1  | sed -e "s/vsearch //g" -e "s/,.*//")
    END_VERSIONS
    """
}
