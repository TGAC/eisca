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
import sys
from pathlib import Path
import util

logger = util.get_named_logger('QC_CellFitering')


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Perform cell clustering and plot UMAPs of clustering",
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
        "--normalize",
        help="Whether to apply normalization to the data.",
        # action='store_true',
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--regress",
        help="Whether to regress out the variations from the total counts and the percentage of mitochondrial genes expressed.",
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--scale",
        help="Whether to scale the expression to have zero mean and unit variance.",
        action=argparse.BooleanOptionalAction,
    )        
    parser.add_argument(
        "--resolutions",
        type=util.floatlist,
        help="Resolution is used to control number of clusters.",
        default=[0.02, 0.05, 0.1, 0.5],
    )                 
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.h5ad.is_file():
        logger.error(f"The given input file {args.h5ad} was not found!")
        sys.exit(2)

    util.check_and_create_folder(args.outdir)
    path_clustering = Path(args.outdir)
    util.check_and_create_folder(path_clustering)

    adata = anndata.read_h5ad(args.h5ad)

    # save the original counts 
    adata.layers["counts"] = adata.X.copy()

    # Normalization
    if args.normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)


    # Feature selection and dimensionality reduction
    sc.pp.highly_variable_genes(
        adata,
        flavor='seurat', 
        min_mean=0.0125, 
        max_mean=1000000000, 
        min_disp=0.5, 
        max_disp=50, 
        n_bins=20, 
        batch_key="sample"
    )
    #sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

    # save the raw counts
    #s adata.raw = adata
    #s adata = adata[:, adata.var.highly_variable]

    # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
    if args.regress:
        if not(hasattr(adata.obs, 'total_counts') and hasattr(adata.obs, 'pct_counts_mt')):
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)
            # sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)    
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

    # scale the expression to have zero mean and unit variance
    if args.scale:
        sc.pp.scale(adata, max_value=10)

    # Dimensionality reduction
    sc.tl.pca(
        adata, 
        n_comps=50, 
        chunked=False,
        zero_center=False, 
        svd_solver='arpack'
    )

    # find nearest neighbor graph constuction
    sc.pp.neighbors(
        adata, 
        n_neighbors=15, 
        n_pcs=50,
        knn=True, 
        method='umap', 
        metric='euclidean'
    )
    sc.tl.umap(adata)

    # perform clustering using Leiden graph-clustering method
    # sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=args.resolution)
    for sid in adata.obs['sample'].unique():
        adata_s = adata[adata.obs['sample']==sid]   
        path_clustering_s = Path(path_clustering, f"sample_{sid}")
        util.check_and_create_folder(path_clustering_s)

        for res in args.resolutions:
            sc.tl.leiden(
                adata_s, n_iterations=2, 
                key_added=f"leiden_res_{res:4.2f}", resolution=res
            )
            with plt.rc_context():
                sc.pl.umap(
                    adata_s,
                    color=[f"leiden_res_{res:4.2f}"],
                    legend_loc="on data",
                    show=False
                )
                plt.savefig(Path(path_clustering_s, f"umap_leiden_res_{res:4.2f}.png"), bbox_inches="tight")

 
    # save a filtered and normalized h5ad file
    adata.write_h5ad(Path(path_clustering, 'adata_clustering.h5ad'))

if __name__ == "__main__":
    sys.exit(main())
