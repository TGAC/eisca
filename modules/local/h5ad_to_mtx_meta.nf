process H5AD_TO_MTX_META {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy-scripts:1.9.301--pyhdfd78af_0' :
        'biocontainers/scanpy-scripts:1.9.301--pyhdfd78af_0' }"

    input:
    path(h5ad_annotated)
    // path samplesheet

    output:
    path "h5ad_mtx_meta"
    path "h5ad_mtx_meta/counts.mtx",  emit: h5ad_mtx
    path "h5ad_mtx_meta/metadata.csv",  emit: h5ad_meta
    path "h5ad_mtx_meta/genes.csv",  emit: h5ad_genes
    path "h5ad_mtx_meta/cells.csv",  emit: h5ad_cells
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    """
    h5ad_to_mtx_meta.py \\
        --h5ad ${h5ad_annotated} \\
        --outdir h5ad_mtx_meta \\


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
