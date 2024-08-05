/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_eisca_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EISCA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_genome_fasta = params.fasta ? file(params.fasta, checkIfExists: true) : ( params.genome ? file( getGenomeAttribute('fasta'), checkIfExists: true ) : [] )
    ch_gtf          = params.gtf   ? file(params.gtf  , checkIfExists: true) : ( params.genome ? file( getGenomeAttribute('gtf')  , checkIfExists: true ) : [] )

    // general input and params
    ch_transcript_fasta = params.transcript_fasta ? file(params.transcript_fasta): []
    ch_motifs = params.motifs ? file(params.motifs) : []
    ch_cellrangerarc_config = params.cellrangerarc_config ? file(params.cellrangerarc_config) : []
    ch_txp2gene = params.txp2gene ? file(params.txp2gene) : []
    ch_multiqc_files = Channel.empty()
    if (params.barcode_whitelist) {
        ch_barcode_whitelist = file(params.barcode_whitelist)
    } else if (protocol_config.containsKey("whitelist")) {
        ch_barcode_whitelist = file("$projectDir/${protocol_config['whitelist']}")
    } else {
        ch_barcode_whitelist = []
    }

    //kallisto params
    ch_kallisto_index = params.kallisto_index ? file(params.kallisto_index) : []
    kb_workflow = params.kb_workflow
    kb_t1c = params.kb_t1c ? file(params.kb_t1c) : []
    kb_t2c = params.kb_t2c ? file(params.kb_t2c) : []

    // samplesheet - this is passed to the MTX conversion functions to add metadata to the
    // AnnData objects.
    ch_input = file(params.input)

    //kallisto params
    ch_kallisto_index = params.kallisto_index ? file(params.kallisto_index) : []
    kb_workflow = params.kb_workflow

    //salmon params
    ch_salmon_index = params.salmon_index ? file(params.salmon_index) : []

    //star params
    ch_star_index = params.star_index ? file(params.star_index) : []
    star_feature = params.star_feature


    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


















    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
