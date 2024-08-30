#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import scrublet as scr
from scipy import io
import anndata
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
import util
# from .util import get_named_logger,

logger = util.get_named_logger('QC_CellFitering')



def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Quality control before and after cell filtering",
    )
    parser.add_argument(
        "--h5ad",
        metavar="FILE_H5AD",
        type=Path,
        help="Input anndata data file.",
        required=True,
    )
    parser.add_argument(
        "--outdir",
        metavar="OUT_DIR",
        type=Path,
        help="Output directory.",
        required=True,
    )
    parser.add_argument(
        "--min_genes",
        type=int,
        help="Filter cells by minimum number of genes.",
        default=1,
    )
    parser.add_argument(
        "--min_counts",
        type=int,
        help="Filter cells by minimum number of counts.",
        default=750,
    )    
    parser.add_argument(
        "--min_cells",
        type=int,
        help="Filter genes by upper quantile of number of genes.",
        default=3,
    )
    parser.add_argument(
        "--pct_mt",
        type=int,
        help="Filter genes by the maximum percentage of mitochondrial counts.",
        default=35,
    )
    parser.add_argument(
        "--doublet_rate",
        type=float,
        help="The expected fraction of transcriptomes that are doublets.",
        default=0.1,
    )        
    parser.add_argument(
        "--quantile_upper",
        type=float,
        help="Filter genes by upper limit of quantile on number of genes.",
        default=1, # 0.98
    )
    parser.add_argument(
        "--quantile_lower",
        type=float,
        help="Filter genes by lower limit of quantile on number of genes.",
        default=0, # 0.02
    )            
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.h5ad.is_file():
        logger.error(f"The given input file {args.h5ad} was not found!")
        sys.exit(2)

    util.check_and_create_folder(args.outdir)
    # path_quant_qc = Path(args.outdir, 'quant_qc')
    path_quant_qc = Path(args.outdir)
    path_quant_qc_scatter = Path(path_quant_qc, 'scatter')
    path_quant_qc_violin = Path(path_quant_qc, 'violin')
    
    path_cell_filtering = Path(args.outdir, 'cell_filtering')
    path_cell_filtering_dist = Path(path_cell_filtering, 'distribution')
    util.check_and_create_folder(path_quant_qc)
    util.check_and_create_folder(path_quant_qc_scatter)
    util.check_and_create_folder(path_quant_qc_violin)
    util.check_and_create_folder(path_cell_filtering)

    adata = anndata.read_h5ad(args.h5ad)
    summary = []
    summary_filtered = []

    # import pdb; pdb.set_trace() #debug code
    # QC on raw counts
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    obs_raw = adata.obs

    # create summary csv file for all samples
    for sid in adata.obs['sample'].unique():
        adata_s = adata[adata.obs['sample']==sid]
        summary += [{
            'Sample ID': sid,
            'Number of cells': adata_s.obs[adata_s.obs['n_genes_by_counts']>0].shape[0],
            'Number of genes': adata_s.var[adata_s.var['n_cells_by_counts']>0].shape[0],
            'Median genes per cell': np.median(adata_s.obs['n_genes_by_counts']),
            'Median of pct-mt':  np.median(adata_s.obs['pct_counts_mt']),
        }]
    summary = pd.DataFrame(summary, index=adata.obs['sample'].unique())
    summary.to_csv(Path(path_quant_qc, 'sample_summary.csv'), index=False)
        

    for sid in adata.obs['sample'].unique():
        adata_s = adata[adata.obs['sample']==sid]
        path_sample = Path(path_quant_qc_scatter, f"sample_{sid}")
        util.check_and_create_folder(path_sample)

        with plt.rc_context():
            sc.pl.scatter(adata_s, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
            plt.savefig(Path(path_sample, 'scatter_total_counts_genes.png'), bbox_inches="tight")

    for sid in adata.obs['sample'].unique():
        adata_s = adata[adata.obs['sample']==sid]
        path_sample = Path(path_quant_qc_violin, f"sample_{sid}")
        util.check_and_create_folder(path_sample)
        with plt.rc_context():
            sc.pl.violin(
                adata_s,
                ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                jitter=0.4,
                multi_panel=True,
                show=False
            )
            plt.savefig(Path(path_sample, 'violin_total_counts_genes_mt.png'), bbox_inches="tight")


    # Cell filtering
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    if args.quantile_upper < 1:
        upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, args.quantile_upper)
        adata = adata[adata.obs.n_genes_by_counts.values < upper_lim]
    if args.quantile_lower > 0:
        lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, args.quantile_lower)
        adata = adata[adata.obs.n_genes_by_counts.values > lower_lim]
    
    adata = adata[adata.obs.pct_counts_mt < args.pct_mt]


    # Doublet detection
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=args.doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets


    # create summary csv file for all samples after filtering
    for sid in adata.obs['sample'].unique():
        adata_s = adata[adata.obs['sample']==sid]
        obs_s = adata_s.obs[adata_s.obs['n_genes_by_counts']>0]
        if 'predicted_doublet' in adata_s.obs.columns:
            obs_s = obs_s[~obs_s['predicted_doublet']]
        n_cells = obs_s.shape[0]
        n_genes = adata_s.var[adata_s.var['n_cells_by_counts']>0].shape[0]
        summary_filtered += [{
            'Sample ID': sid,
            'Number of cells': f"{n_cells} ({n_cells/summary.loc[sid, 'Number of cells']:.0%})",
            'Number of genes': f"{n_genes} ({n_genes/summary.loc[sid, 'Number of genes']:.0%})",
            'Median genes per cell': np.median(obs_s['n_genes_by_counts']),
            'Median of pct-mt':  np.median(obs_s['pct_counts_mt']),
        }]
    summary_filtered = pd.DataFrame(summary_filtered, index=adata.obs['sample'].unique()) 
    summary_filtered.to_csv(Path(path_cell_filtering, 'sample_summary_filtered.csv'), index=False)


    # distributions of n_genes_by_counts, total_counts and pct_counts_mt before and after filtering
    for sid in adata.obs['sample'].unique():
        obs1 = obs_raw[obs_raw['sample']==sid]    
        obs2 = adata.obs[adata.obs['sample']==sid] 
        if 'predicted_doublet' in obs2.columns:
            obs2 = obs2[~obs2['predicted_doublet']]
        path_sample = Path(path_cell_filtering_dist, f"sample_{sid}")
        util.check_and_create_folder(path_sample)

        hist_df = pd.concat([obs1[['n_genes_by_counts']], obs2[['n_genes_by_counts']]], axis=1)
        hist_df.columns = ['Before filtering', 'After filtering']
        plot = sns.displot(data=hist_df, kind="kde",  fill=True, palette=sns.color_palette('bright'), aspect=1.2)
        plot.set(xlabel='n_genes_by_counts')
        plot.add_legend(loc="upper right")
        plot.savefig(Path(path_sample, 'dist1_n_genes_by_counts.png'))

        hist_df = pd.concat([obs1[['total_counts']], obs2[['total_counts']]], axis=1)
        hist_df.columns = ['Before filtering', 'After filtering']
        plot = sns.displot(data=hist_df, kind="kde",  fill=True, palette=sns.color_palette('bright'), aspect=1.2)
        plot.set(xlabel='total_counts')
        plot.add_legend(loc="upper right")
        plot.savefig(Path(path_sample, 'dist2_total_counts.png'))

        hist_df = pd.concat([obs1[['pct_counts_mt']], obs2[['pct_counts_mt']]], axis=1)
        hist_df.columns = ['Before filtering', 'After filtering']
        plot = sns.displot(data=hist_df, kind="kde",  fill=True, palette=sns.color_palette('bright'), aspect=1.2)
        plot.set(xlabel='pct_counts_mt')
        plot.add_legend(loc="upper right")
        plot.savefig(Path(path_sample, 'dist3_pct_counts_mt.png'))


    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # save a filtered and normalized h5ad file
    adata.write_h5ad(Path(path_quant_qc, 'adata_filtered_normalized.h5ad'))

    # Feature selection and dimensionality reduction
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    with plt.rc_context():
        sc.pl.highly_variable_genes(adata, show=False)
        plt.savefig(Path(path_cell_filtering, 'highly_variable_genes.png'), bbox_inches="tight")
    sc.tl.pca(adata)

    # QC after cell filtering
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    with plt.rc_context():
        sc.pl.umap(
            adata,
            color="sample",
            size=2,
            show=False
        )
        plt.savefig(Path(path_cell_filtering, 'umap_samples.png'), bbox_inches="tight")


    for sid in adata.obs['sample'].unique():
        adata_s = adata[adata.obs['sample']==sid]   
        path_cell_filtering_s = Path(path_cell_filtering, f"sample_{sid}")
        util.check_and_create_folder(path_cell_filtering_s)
        with plt.rc_context():
            sc.pl.umap(
                adata_s,
                color=["predicted_doublet", "doublet_score"],
                wspace=0.3, # increase horizontal space between panels
                size=3,
                show=False
            )
            plt.savefig(Path(path_cell_filtering_s, 'umap_doublet.png'), bbox_inches="tight")

        with plt.rc_context():
            sc.pl.umap(
                adata_s,
                color=["log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
                wspace=0.3,
                ncols=2,
                show=False
            )
            plt.savefig(Path(path_cell_filtering_s, 'umap_total_counts_genes_mt.png'), bbox_inches="tight")    


if __name__ == "__main__":
    sys.exit(main())
