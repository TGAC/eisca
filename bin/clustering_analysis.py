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
        action='store_true',
        # action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--keep_doublets",
        help="Whether to filter out the cells called as doublets.",
        action='store_true',
    )       
    parser.add_argument(
        "--regress",
        help="Whether to regress out the variations from the total counts and the percentage of mitochondrial genes expressed.",
        action='store_true',
    )
    parser.add_argument(
        "--scale",
        help="Whether to scale the expression to have zero mean and unit variance.",
        action='store_true',
    )        
    parser.add_argument(
        "--resolutions",
        type=util.floatlist,
        help="Resolution is used to control number of clusters.",
        default=[0.02, 0.05, 0.1, 0.5],
    )
    parser.add_argument(
        "--integrate",
        default=None,
        choices=['bbknn', 'harmony', 'scanorama'],
        help="Choose a method for data integration",
    )
    parser.add_argument(
        "--meta",
        default='auto',
        choices=['auto', 'sample', 'group', 'plate'],
        help="Choose a metadata column as the batch for clustering",
    )
    parser.add_argument(
        "--fontsize",
        type=int,
        help="Set font size for plots.",
        default=12,
    )                     
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.h5ad.is_file():
        logger.error(f"The given input file {args.h5ad} was not found!")
        sys.exit(2)

    plt.rcParams.update({
        "font.size": args.fontsize,
        # "axes.titlesize": 'medium',
        # "axes.labelsize": 'small',
        # "xtick.labelsize": 'small',
        # "ytick.labelsize": 'small',
        # "legend.fontsize": 'small',
    })

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
    if not adata.uns.get('log1p'): # to fix issue in scanpy function
        adata.uns['log1p'] = {'base': None}

    # remove doublets before clustering
    if not args.keep_doublets:
        if hasattr(adata.obs, 'predicted_doublet'):
            adata = adata[~adata.obs['predicted_doublet']]

    batch_key = 'plate' if hasattr(adata.obs, 'plate') else 'sample' # correct on plates for smart-seq data

    # Feature selection and dimensionality reduction
    sc.pp.highly_variable_genes(
        adata,
        flavor='seurat', 
        min_mean=0.0125, 
        max_mean=1000000000, 
        min_disp=0.5, 
        max_disp=50, 
        n_bins=20, 
        batch_key=batch_key
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
    n_comps = min((min(adata.X.shape)-1), 50)
    sc.tl.pca(
        adata, 
        n_comps=n_comps, 
        chunked=False,
        zero_center=False, 
        svd_solver='arpack'
    )

    # perform data integration
    if args.meta == 'auto':
        # batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
        batch = 'sample'
        if hasattr(adata.obs, 'group'):
            batch = 'group'
        elif hasattr(adata.obs, 'plate'):
            batch = 'plate' 
    else:
        batch = args.meta
    
    n_pcs = adata.obsm['X_pca'].shape[1]
    if args.integrate == 'bbknn':
        sc.external.pp.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs)
    elif args.integrate == 'harmony':
        sc.external.pp.harmony_integrate(adata, batch_key)
    elif args.integrate == 'scanorama':
        sc.external.pp.scanorama_integrate(adata, batch_key)
        
    # find nearest neighbor graph constuction
    pca_rep = 'X_pca'
    if args.integrate == 'harmony':
        pca_rep = 'X_pca_harmony'
    elif args.integrate == 'scanorama':
        pca_rep = 'X_scanorama'
    if args.integrate != 'bbknn':
        sc.pp.neighbors(
            adata, 
            n_neighbors=15, 
            n_pcs=n_pcs,
            knn=True, 
            method='umap', 
            metric='euclidean',
            use_rep=pca_rep
        )
    sc.tl.umap(adata)

    # perform clustering using Leiden graph-clustering method
    for res in args.resolutions:
        sc.tl.leiden(
            adata, n_iterations=2, 
            key_added=f"leiden_res_{res:4.2f}", resolution=res
        )

    # save the AnnData into a h5ad file 
    adata.write_h5ad(Path(path_clustering, 'adata_clustering.h5ad'))

    
    # sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=args.resolution)
    for sid in sorted(adata.obs[batch].unique()):
        adata_s = adata[adata.obs[batch]==sid]   
        path_clustering_s = Path(path_clustering, f"{batch}_{sid}")
        util.check_and_create_folder(path_clustering_s)
        for res in args.resolutions:
            with plt.rc_context():
                sc.pl.umap(
                    adata_s,
                    color=[f"leiden_res_{res:4.2f}"],
                    # legend_loc="on data",
                    show=False
                )
                plt.savefig(Path(path_clustering_s, f"umap_leiden_res_{res:4.2f}.png"), bbox_inches="tight")

    
    # stacked proportion bar plot showing all samples for each resolution      
    for res in args.resolutions:
        path_res = Path(path_clustering, f"resolution_{res:4.2f}")
        util.check_and_create_folder(path_res)
        n_cluster = len(adata.obs[f'leiden_res_{res:4.2f}'].unique())+1
        ncol = min((n_cluster//20 + min(n_cluster%20, 1)), 3)
        with plt.rc_context():
            prop = pd.crosstab(adata.obs[f'leiden_res_{res:4.2f}'],adata.obs[batch], normalize='columns').T.plot(kind='bar', stacked=True)
            prop.legend(bbox_to_anchor=(1+(args.fontsize-10)/50+ncol*0.17, 1.02), loc='upper right', ncol=ncol)
            plt.savefig(Path(path_res, f"prop_leiden_res_{res:4.2f}.png"), bbox_inches="tight")


    # save analysis parameters into a json file
    with open(Path(path_clustering, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})        
        params.update({"--resolutions": args.resolutions})        
        if args.integrate: params.update({"--integrate": args.integrate})        
        params.update({"--meta": args.meta})        
        if args.normalize: params.update({"--normalize": args.normalize})
        if args.regress: params.update({"--regress": args.regress})
        if args.scale: params.update({"--scale": args.scale})
        if args.keep_doublets: params.update({"--keep_doublets": args.keep_doublets})
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())
