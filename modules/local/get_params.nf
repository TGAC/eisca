import groovy.json.JsonBuilder

process GET_PARAMS {
    label "eisca"
    label 'process_low'

    output:
        path "params.json",  emit: json
    
    script:
    String paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}
