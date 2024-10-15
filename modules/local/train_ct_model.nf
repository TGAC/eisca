process TRAIN_CT_MODEL {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/celltypist:1.6.3--pyhdfd78af_0' :
        'teichlab/celltypist:latest' }"

    input:
    path h5ad_filtered

    output:
    path "annotation"
    path "annotation/models/*.pkl",  emit: model
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    train_ct_model.py \\
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
