# nf-core/eisca: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.6 - [2024-11-05]
- add following tertiary analysis: Differential expression analysis
- fix plot size issue in report


## v1.4 - [2024-10-22]
- add following tertiary analysis: Cell type annotation
- add 'analysis parameters' for analysis sections in report 

## v1.1 - [2024-10-04]
- add module 'Merging/integration of samples' for secondary analysis

## v1.0 - [2024-09-19]

Initial release of nf-core/eisca, created with the [nf-core](https://nf-co.re/) template.
The pipeline including following analyses:
1. Primary analysis:
   - Quality control of raw reads
   - Mapping reads and quantification
   - Convert count matrix into Anndata and Securat objects
2. Secondary analysis
   - Single-cell quality control
   - Cell filtering
   - Clustering analysis

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
