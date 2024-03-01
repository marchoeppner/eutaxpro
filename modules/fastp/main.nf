process FASTP {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${meta.sample_id}/FASTP", mode: 'copy'

    label 'short_parallel'

    conda 'bioconda::fastp=0.23.4'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'quay.io/biocontainers/fastp:0.23.4--hadf994f_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*trimmed.fastq.gz'), emit: reads
    path("*.json"), emit: json
    path('versions.yml'), emit: versions

    script:
    prefix = meta.sample_id
    suffix = '.trimmed.fastq.gz'

    json = prefix + '.fastp.json'
    html = prefix + '.fastp.html'

    if (meta.single_end) {
    } else {
        r1 = reads[0]
        r2 = reads[1]
        r1_trim = prefix + '_1' + suffix
        r2_trim = prefix + '_2' + suffix
        """
        fastp -c --in1 $r1 --in2 $r2 \
        --out1 $r1_trim \
        --out2 $r2_trim \
        --detect_adapter_for_pe \
        -w ${task.cpus} \
        -j $json \
        -h $html \
        -n 0 \
        --length_required 35

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }
}
