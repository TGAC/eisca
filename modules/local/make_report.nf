
process MAKE_REPORT {
    // tag "$meta.id"
    label 'process_low'

    conda "bioconda::pysam=0.19.0 bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0' :
        'biocontainers/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0' }"

    input:
    path "versions.txt"
    path "params.json"

    output:
    path "eisca_report.html"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    report.py \\
        --report-fname eisca_report.html \\
        --results ${params.outdir} \\
        --versions versions.txt \\
        --params params.json \\
        --wf-version ${workflow.manifest.version} \\


    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     python: \$(python --version | sed 's/Python //g')
    // END_VERSIONS
    """
}
