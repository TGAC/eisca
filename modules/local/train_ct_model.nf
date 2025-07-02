process TRAIN_CT_MODEL {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/celltypist:1.6.3--pyhdfd78af_0' :
        'docker.io/myeihub/celltypist_ps:1.6.3' }"

    input:
    path h5ad_filtered

    output:
    path "ctmodel"
    path "ctmodel/*.pkl",  emit: model
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    train_ct_model.py \\
        --h5ad ${h5ad_filtered} \\
        --outdir ctmodel \\
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
