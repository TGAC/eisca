#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import pandas as pd
import numpy as np
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
        default=100,
    )
    parser.add_argument(
        "--min_cells",
        type=int,
        help="Filter genes by upper quantile of number of genes.",
        default=3,
    )
    # parser.add_argument(
    #     "--quantile_upper",
    #     type=float,
    #     help="Filter genes by upper limit of quantile on number of genes.",
    #     default=0.98,
    # )
    # parser.add_argument(
    #     "--quantile_lower",
    #     type=float,
    #     help="Filter genes by lower limit of quantile on number of genes.",
    #     default=0.02,
    # )
    parser.add_argument(
        "--pct_mt",
        type=int,
        help="Filter genes by the maximum percentage of mitochondrial counts.",
        default=20,
    )
    # parser.add_argument(
    #     "--resolution",
    #     type=float,
    #     help="Resulution is used to control number of clusters.",
    #     default=0.05,
    # )                 
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
    util.check_and_create_folder(path_quant_qc)
    util.check_and_create_folder(path_quant_qc_scatter)
    util.check_and_create_folder(path_quant_qc_violin)
    util.check_and_create_folder(path_cell_filtering)

    adata = anndata.read_h5ad(args.h5ad)
    summary = []

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

    # create summary csv file for all samples
    for sid in adata.obs['sample'].unique():
        adata_s = adata[adata.obs['sample']==sid]
        summary += [{
            'sample_id': sid,
            'n_cells': adata_s.obs[adata_s.obs['n_genes_by_counts']>0].shape[0],
            'n_genes': adata_s.var[adata_s.var['n_cells_by_counts']>0].shape[0],
            'gpc_median': np.median(adata_s.obs['n_genes_by_counts']),
            'pctmt_median':  np.median(adata_s.obs['pct_counts_mt']),
        }]
    summary = pd.DataFrame(summary) 
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
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    # upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, args.quantile_upper)
    # lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, args.quantile_lower)
    # adata = adata[(adata.obs.n_genes_by_counts.values < upper_lim) & (adata.obs.n_genes_by_counts.values > lower_lim)]
    
    adata = adata[adata.obs.pct_counts_mt < args.pct_mt]

    # Doublet detection
    # sc.pp.scrublet(adata, batch_key="sample")

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # save a filtered and normalized h5ad file
    adata.write_h5ad(Path(path_quant_qc, 'cell_filtered_normalized.h5ad'))

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
        # with plt.rc_context():
        #     sc.pl.umap(
        #         adata_s,
        #         color=["predicted_doublet", "doublet_score"],
        #         wspace=0.5, # increase horizontal space between panels
        #         size=3,
        #         show=False
        #     )
        #     plt.savefig(Path(path_cell_filtering_s, 'umap_doublet.png'), bbox_inches="tight")

        with plt.rc_context():
            sc.pl.umap(
                adata_s,
                color=["log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
                wspace=0.5,
                ncols=2,
                show=False
            )
            plt.savefig(Path(path_cell_filtering_s, 'umap_total_counts_genes_mt.png'), bbox_inches="tight")    


if __name__ == "__main__":
    sys.exit(main())
