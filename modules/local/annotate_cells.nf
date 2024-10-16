process ANNOTATE_CELLS {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/celltypist:1.6.3--pyhdfd78af_0' :
        'biocontainers/celltypist:1.6.3--pyhdfd78af_0' }"

    input:
    path h5ad_filtered
    path model_file

    output:
    path "annotation"
    path "annotation/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if(params.ctmodel)
    """
    annotate_cells.py \\
        --h5ad ${h5ad_filtered} \\
        --model_file ${model_file} \\
        --outdir annotation \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS        
    """
    
    else
    """
    annotate_cells.py \\
        --h5ad ${h5ad_filtered} \\
        --outdir annotation \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS        
    """    



    // stub:
    // """
    // touch combined_matrix.h5ad
    // """
}
