process SINTAX_OTU2TAB {
    tag "$meta  "
    label 'process_low'

    conda 'conda-forge::perl=5.32'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta),path(sintax),path(otu_tab)

    output:
    tuple val(meta), path(result), emit: tsv
    path 'versions.yml'    , emit: versions

    script:
    def args = task.ext.args ?: ''
    result = meta.sample_id + ".taxonomy_by_sample.tsv"
    
    """
    sintax_otu2tab.pl --sintax $sintax --otu $otu_tab --outfile $result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version  | head -n2 | tail -n1 | sed -e "s/.*(//" -e "s/).*//")
    END_VERSIONS
    """
}
