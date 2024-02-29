process BIOBLOOMTOOLS_CATEGORIZER {
    publishDir "${params.outdir}/Processing/Bloomfilter", mode: 'copy'

    label 'short_parallel'

    tag "${meta.sample_id}|${meta.library_id}|${meta.readgroup_id}"

    conda 'bioconda::biobloomtools=2.3.5'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biobloomtools:2.3.5--h4056dc3_2' :
        'quay.io/biocontainers/biobloomtools:2.3.5--h4056dc3_2' }"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path(r1_trim), path(r2_trim), emit: reads
    path('versions.yml'), emit: versions
    path("*summary.tsv"), emit: results

    script:
    filtered = meta.sample_id + '_' + meta.library_id + '_' + meta.readgroup_id
    r1_trim = filtered + '_noMatch_1.fq.gz'
    r2_trim = filtered + '_noMatch_2.fq.gz'

    """
    biobloomcategorizer -p $filtered -t ${task.cpus} -n --fq --gz_out -i -e -f "${params.bloomfilter}" $r1 $r2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BioBloomtools: \$(biobloomcategorizer -version 2>&1 | head -n1 | sed -e "s/.*) //g")
    END_VERSIONS

    """
}
