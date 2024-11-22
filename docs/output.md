# nf-core/eisca: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Primary analysis](#primary-analysis)
  - [FastQC](#fastqc) - Raw read QC
  - [Kallisto & Bustools](#kallisto--bustools) - Mapping & quantification by Kallisto & Bustools
  - [Salmon Alevin](#salmon-alevin--alevinqc) - Mapping & quantification by Salmon Alevin
  - [STARsolo](#starsolo) - Mapping & quantification by STAR
- [Secondary analysis](#secondary-analysis)
  - [QC & cell filtering](#qc--cell-filtering) - Cell filtering and QC on raw data and filtered data
  - [Clustering analysis](#clustering-analysis) - Single-cell clustering analysis
- Tertiary analysis
  - [Cell-type annotation analysis](#annotation-analysis) - Single-cell cell-type annotation analysis
  - [Differential expression analysis](#dea-analysis) - Single-cell differential expression analysis
- [Pipeline reporting](#pipeline-reporting)
  - [Analysis report](#analysis-report) - Single-ell Analysis Report
  - [MultiQC](#multiqc) - Aggregate report describing results and QC for tools registered in nf-core
  - [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Primary analysis

### <u>FastQC</u>

<!-- <details markdown="1"> -->
<!-- <summary>Output files</summary> -->

**Output directory: `results/fastqc`**
- `*_fastqc.html`: FastQC report containing quality metrics.
- `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

<!-- </details> -->

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

<!-- ![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png) -->

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.
:::


### <u>Kallisto & Bustools</u>

**Output directory: `results/kallisto`**
- `<sample_name>.count/`: Contains count matrix and mapping information.
  - `counts_unfiltered/`
    - `cells_x_genes.mtx`: unfiltered count matrix.
  - `run_info.json`: json file showing mapping statistics.  
  - `*.bus`: Contains the same BUS formatted data, sorted and corrected with the supplied barcode whitelist.
- `mtx_conversions/`
  - `<sample_name>/`
    - `sample_name_*_matrix.h5ad`: AnnData object file for this sample.
    - `sample_name_*_matrix.rds`: Seurat object file for this sample.
  - `combined_raw_matrix.h5ad`: AnnData object file for combined samples.

**Output directory: `results/reference_genome`**
- `kallisto_index/`: Contains the index of the supplied (genome/transcriptome) fasta file.

The kallisto & Bustools workflow can analyze data from single cell rnaseq experiments and generates a set of folders with respective outputs from various steps of the analysis. For a detailed summary what the pipeline does specifically, please follow the [excellent tutorial](https://www.kallistobus.tools/getting_started.html) that also describes specific steps for downstream analysis of the generated matrices.

For details on how to load these into R and perform further downstream analysis, please refer to the [BusTools HowTo](https://github.com/BUStools/getting_started/blob/master/getting_started.ipynb). See [Kallisto](https://pachterlab.github.io/kallisto/about) for details about Kallisto and [Bustools](https://bustools.github.io/) for more information on BusTools.


### <u>Salmon Alevin & AlevinQC</u>

**Output directory: `results/alevin`**
- `<sample_name>_alevin_results/`: Contains intermediate results.
  - `af_quant/`
    - 'quants_mat.mtx`: unfiltered count matrix.
- `mtx_conversions/`
  - `<sample_name>/`
    - `sample_name_*_matrix.h5ad`: AnnData object file for this sample.
    - `sample_name_*_matrix.rds`: Seurat object file for this sample.
  - `combined_raw_matrix.h5ad`: AnnData object file for combined samples.

**Output directory: `results/alevinqc`**
- `alevin_report_<sample_name>.html`: QC report for the Salmon Alevin output data.

**Output directory: `results/reference_genome`**
- `salmon_index/`: Contains the indexed reference transcriptome for Salmon Alevin.
- `alevin/txp2gene.tsv`: The transcriptome to gene mapping TSV file utilized by Salmon Alevin.


### <u>STARsolo</u>

**Output directory: `results/star`**
- Files will be organized in one directory per sample.
- Contains the mapped BAM files and output metrics created by STARsolo.

**Output directory: `results/reference_genome`**
- `star_index/`: Contains the index of the supplied genome fasta file.


## Secondary analysis

### <u>QC & cell filtering</u>

**Output directory: `results/qc_cell_filter`**
- `sample_summary.csv`: overall summary of the single-cell count matrix
- `adata_filtered_normalized.h5ad`: AnnData object file after cell filtering and normalization
- `raw_counts/sample_*/`
  - `scatter_total_counts_genes.png`: scatter plot shows the relationship between total read counts and the number of genes.
  - `violin*.png`: violin plots display the distribution of cells based on the number of genes, total counts, and the percentage of counts in mitochondrial genes.  
- `cell_filtering/`
  - `highly_variable_genes.png`: plot of mean expressions against dispersions of genes for highly variable genes.
  - `umap_samples.png`: UMAP plot between samples.
  - `sample_summary_filtered.csv`: overall summary of the single-cell count matrix after cell filtering
  - `sample_*/`
    - `umap_total_counts_genes_mt.png`: UMAP plots for the number of genes, total counts, and the percentage of counts in mitochondrial genes.
    - `violin*.png`: violin plots display the distribution of cells based on the number of genes, total counts, and the percentage of counts in mitochondrial genes after cell filtering.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.
    

### <u>Clustering analysis</u>

**Output directory: `results/clustering`**
- `adata_clustering.h5ad`: AnnData object file after clustering analysis.
- `sample_*/` or `group_*/`
  - `umap_leiden_res_*.png`: UMAP plots showing clustering results with differnt resoultuion settings.
- `resolution_*/`
  - `prop_leiden_res_*.png`: plot showing a stacked bar chart that presents the proportions of clusters across samples/groups.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


### <u>Annotation analysis</u>

**Output directory: `results/annotation`**
- `adata_annotation.h5ad`: AnnData object file after cell-type annotation analysis.
- `sample_*/` or `group_*/`
  - `umap_cell_type.png`: UMAP plots showing predicted cell-type clusters.
  - `umap_conf_score.png`: UMAP plots showing mapped confidence scores of the cells.
- `prop_majority_voting.png`: plot showing a stacked bar chart that presents the proportions of cell-type clusters across samples/groups.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


### <u>DEA analysis</u>

**Output directory: `results/dea`**
- `adata_annotation.h5ad`: AnnData object file after cell-type annotation analysis.
- `sample_*/` or `group_*/` or `celltype_*/` (no subfolder for DEA betweeen groups)
  - `plot_genes_*.png`: plots showing top number of DE genes across groups.
  - `dotplot_genes_*.png`: dot plot showing top number of DE genes across groups.
  - `dea_*.csv`: a csv table file showing DEA results for all genes, e.g. log fold change, p-values.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


## Pipeline reporting

### <u>Analysis report</u>

**Output directory: `results/report`**
- `eisca_report.html`: this is HTML report file showing all major analysis results.


### <u>MultiQC</u>

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### <u>Pipeline information</u>

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
