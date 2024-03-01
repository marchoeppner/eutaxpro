process PTRIMMER {
    publishDir "${params.outdir}/${meta.sample_id}/PTRIMMER", mode: 'copy'

    label 'short_serial'

    tag "${meta.sample_id}"

    conda 'bioconda::ptrimmer=1.3.3.'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ptrimmer:1.3.3--h50ea8bc_5' :
        'quay.io/biocontainers/ptrimmer:1.3.3--h50ea8bc_5' }"

    input:
    tuple val(meta), path(reads)
    path(amplicon_txt)

    output:
    tuple val(meta), path('*ptrimmed.fastq.gz'), emit: reads
    path('versions.yml'), emit: versions

    script:

    if (meta.single_end) {
    } else {
        r1 = reads[0]
        r2 = reads[1]
        r1_trimmed = meta.sample_id + '_1.ptrimmed.fastq'
        r2_trimmed = meta.sample_id + '_2.ptrimmed.fastq'
        r1_trimmed_gz = r1_trimmed + '.gz'
        r2_trimmed_gz = r2_trimmed + '.gz'

        """
        ptrimmer -t pair -a $amplicon_txt -f $r1 -d $r1_trimmed -r $r2 -e $r2_trimmed
        gzip $r1_trimmed
        gzip $r2_trimmed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            Ptrimmer: \$(ptrimmer --help 2>&1 | grep Version | sed -e "s/Version: //g")
        END_VERSIONS

        """
    }
}
