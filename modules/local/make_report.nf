import groovy.json.JsonBuilder

process MAKE_REPORT {
    // tag "$meta.id"
    label 'process_low'
    // debug true

    // conda "bioconda::pysam=0.19.0 bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/ontresearch/wf-common:sha338caea0a2532dc0ea8f46638ccc322bb8f9af48' :
        'docker.io/ontresearch/wf-common:sha338caea0a2532dc0ea8f46638ccc322bb8f9af48' }"

    input:
    // val  ready
    // path "versions.txt"
    path results
    path multiqc_report // just for executing in the end of pipeline
    // path "params.json"

    output:
    path "eisca_report.html"

    when:
    task.ext.when == null || task.ext.when

    script:
    report_name = "eisca_report.html"
    String paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    echo '$paramsJSON' > params.json
    
cat << EOF > versions.txt
Python,3.9.18
pandas,2.2.0
scanpy,1.9.3 
anndata,0.10.5.post1
EOF

    report.py \\
        $report_name \\
        --results ${results} \\
        --versions versions.txt \\
        --params params.json \\
        --wf-version ${workflow.manifest.version} \\
        --logo ${workflow.projectDir}/bin/images/EI_logo.png \\


    """
}
