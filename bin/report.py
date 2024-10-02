#!/usr/bin/env python

import dominate.tags as html
import ezcharts as ezc
from ezcharts.components.reports.labs import LabsReport, LabsNavigation, ILabsNavigationClasses
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.components.theme import LAB_head_resources
from ezcharts.components.ezchart import EZChart

import argparse
from pathlib import Path
from report_util import *
import util
import json
import pandas as pd
import sys


logger = util.get_named_logger('Report')

report_title = 'EI Single Cell Analysis Report'
workflow_name = 'eisca'


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Quality control before and after cell filtering",
        # epilog="python count_reads_from_bam.py --bam file.bam --bed file.bed --json output.json",
    )
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--samplesheet",
        metavar="SAMPLESHEET",
        type=Path,
        help="Input samplesheet file.",
        required=True,
    )
    parser.add_argument(
        "--results",
        metavar="RESULTS_DIR",
        type=Path,
        help="Results directory.",
        required=True,
    )
    # parser.add_argument(
    #     "--stats", nargs='+',
    #     help="Fastcat per-read stats, ordered as per entries in --metadata.")
    # parser.add_argument(
    #     "--images", nargs='+',
    #     help="Sample directories containing various images to put in report")
    # parser.add_argument(
    #     "--survival",
    #     help="Read survival data in TSV format")
    parser.add_argument(
        "--params",
        metavar="FILE_PARAMS",
        type=Path,
        help="Workflow params json file",
        required=True,
    )
    parser.add_argument(
        "--versions",
        metavar="FILE_VERSIONS",
        type=Path,
        help="Workflow versions file",
        required=True,
    )
    parser.add_argument(
        "--logo",
        metavar="FILE_LOGO",
        type=Path,
        help="Logo image file",
        required=True,
    )    
    # parser.add_argument(
    #     "--umap_dirs", nargs='+',
    #     help="Sample directories containing umap and gene expression files")
    # parser.add_argument(
    #     "--umap_genes", help="File containing list of genes to annnotate UMAPs")
    # parser.add_argument(
    #     "--metadata", default='metadata.json', required=True,
    #     help="sample metadata")
    parser.add_argument(
        "--wf-version", default='unknown',
        help="version of the executed workflow")               
    return parser.parse_args(argv)


def main(argv=None):
    logger.info('Building report')
    args = parse_args(argv)


    if not args.params.is_file():
        logger.error(f"The given input file {args.params} was not found!")
        sys.exit(2)

    if not args.versions.is_file():
        logger.error(f"The given input file {args.versions} was not found!")
        sys.exit(2)

    report = LabsReport(
        report_title, workflow_name,
        args.params, args.versions, args.wf_version,
        head_resources=[*LAB_head_resources],
    )
    
    report.nav.getElementsByTagName('a')[0].clear()
    navdiv = report.nav.getElementsByTagName('div')[0]
    link = html.a(href='https://www.earlham.ac.uk/', cls=ILabsNavigationClasses().logo)
    with link:
        # with open(Path('images/EI_logo.png'), 'rb') as f:
        with open(Path(args.logo), 'rb') as f:
            b64img = base64.b64encode(f.read()).decode()
            html.img(src=f'data:image/png;base64,{b64img}', width=120)
    navdiv[0] = link

    report.banner.clear()
    report.footer.clear()
    report.intro_content.add(EIBanner(report_title, workflow_name))
    report.footer.add(EILabsAddendum(workflow_name, args.wf_version))

    samplesheet = pd.read_csv(args.samplesheet)
    batch = 'group' if hasattr(samplesheet, 'group') else 'sample'
    Nbatch = len(samplesheet[batch].unique())
    # samples = samplesheet['sample'].unique()

    path_quant_qc = Path(args.results, 'qc_cell_filter')
    path_quant_qc_raw = Path(path_quant_qc, 'raw_counts')
    # path_quant_qc_scatter = Path(path_quant_qc, 'scatter')
    # path_quant_qc_violin = Path(path_quant_qc, 'violin')
    path_cell_filtering = Path(path_quant_qc, 'cell_filtering')
    # path_cell_filtering_dist = Path(path_cell_filtering, 'distribution')
    path_clustering = Path(args.results, 'clustering')

    # print(path_quant_qc) #tst

    if path_quant_qc.exists():
        summary = pd.read_csv(Path(path_quant_qc, 'sample_summary.csv')).set_index('Sample ID')
        Nsample = summary.shape[0]
        with report.add_section('Single cell summary', 'Summary'):
            html.p("""This section gives an overall summary of the single-cell count matrix for 
                   each sample. The statistics include the total number of cells with at least 
                   one gene expressed, the total number of genes expressed in at least one cell, 
                   the median number of genes per cell, and the median percentage of counts in 
                   mitochondrial genes (pct-mt).""")
            DataTable.from_pandas(summary)   

        with report.add_section('Quantification QC', 'QC'):
            html.p("""This section presents the QC plots of the raw count matrix generated during 
                   the quantification step. These plots provide insight into the quality of the 
                   experiments and guide the filtering of low-quality cells.""")
            html.p("""The following scatter plot shows the relationship between total 
                   read counts and the number of genes, with the percentage of counts in 
                   mitochondrial genes indicated by color.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', widths=['800'], =['scatter*.png'])
            html.p("""The following violin plots display the distribution of cells based on the number of 
                   genes, total counts, and the percentage of counts in mitochondrial genes.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', ncol=3, suffix=['violin*.png'])
    else:
        logger.info('Skipping Quantification QC')

    if path_cell_filtering.exists():
        summary = pd.read_csv(Path(path_cell_filtering, 'sample_summary_filtered.csv')).set_index('Sample ID')
        with report.add_section('Cell filtering', 'Cell filtering'):
            html.p("""This section presents the statistics and QC plots after cell filtering process. 
                   The filtering process includes hard thresholds for minimum genes, minimum counts, 
                   minimum cells and the percentage of counts in mitochondrial genes. Additionally, 
                   users can set quantile limits on the number of genes. Then doublet detection is 
                   performed using Scrublet.""")
            html.p("""The following table shows summary statistics, with percentages in brackets 
                   indicating the comparison to the raw counts.""")
            DataTable.from_pandas(summary)
            html.p("""The following violin plots display the distribution of cells based on the number of 
                   genes, total counts, and the percentage of counts in mitochondrial genes after filtering.""")
            plots_from_image_files(path_cell_filtering, meta='sample', ncol=3, suffix=['violin*.png'])            
            # html.p("""The following plots show the distribution KDE curves before and after filtering 
            #        for the number of genes, total counts, and the percentage of counts in mitochondrial genes.""")            
            # plots_from_image_files(path_cell_filtering, meta='sample', ncol=2, suffix=['dist*.png'])
            html.p("""The following plots show the UMAP plots
                   for the number of genes, total counts, and the percentage of counts in mitochondrial genes.""")                        
            plots_from_image_files(path_cell_filtering, meta='sample', suffix=['umap_total*.png'])
            if util.check_file(f"{path_cell_filtering}/sample_*", 'umap_doublet.png'):
                html.p("""The following plots show the UMAP plots for the predicted doublets and doublet scores.""")                        
                plots_from_image_files(path_cell_filtering, meta='sample', suffix=['umap_doublet.png'])            
    else:
        logger.info('Skipping Cell filtering')


    if path_clustering.exists():
        with report.add_section('Clustering analysis', 'Clustering'):
            html.p("""This section shows clustering UMAP plots for each {batch}. The clustering 
                   was performed using Leiden graph-clustering method. The resolution parameter 
                   was set for different values to get different number of clusters which 
                   could match to biologically-meaningful cell types.""")
            plots_from_image_files(path_clustering, meta=batch, ncol=2)
            html.p("""The following plot shows a stacked bar chart that presents the proportions of clusters 
                   across {batch}s, calculated for each resolution value. The plot illustrates the distribution 
                   profiles of predicted clusters between {batch}s.""")                        
            plots_from_image_files(path_clustering, meta='resolution', widths=[str(min(Nbatch*280, 1200))])                    
    else:
        logger.info('Skipping Cell filtering')        


    report.write(args.report)
    logger.info('Report writing finished')
    

if __name__ == "__main__":
    sys.exit(main())
