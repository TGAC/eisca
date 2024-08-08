#!/usr/bin/env python

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import io
import anndata
from matplotlib import pyplot as plt
import argparse
import logging
import sys
import json
import gzip
from pathlib import Path
import util
# from .util import get_named_logger,

logger = util.get_named_logger('QC_CellFitering')



def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Quality control before and after cell filtering",
        epilog="python count_reads_from_bam.py --bam file.bam --bed file.bed --json output.json",
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

    adata = anndata.read_h5ad(args.h5ad)
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

    with plt.rc_context():
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False
        )
        plt.savefig(Path(args.outdir, 'violin_total_counts_genes_mt.png'), bbox_inches="tight")

    with plt.rc_context():
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        plt.savefig(Path(args.outdir, 'scatter_total_counts_genes.png'), bbox_inches="tight")

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

    # Feature selection and dimensionality reduction
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    with plt.rc_context():
        sc.pl.highly_variable_genes(adata, show=False)
        plt.savefig(Path(args.outdir, 'highly_variable_genes.png'), bbox_inches="tight")
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
        plt.savefig(Path(args.outdir, 'umap_samples.png'), bbox_inches="tight")

    # with plt.rc_context():
    #     sc.pl.umap(
    #         adata,
    #         color=["predicted_doublet", "doublet_score"],
    #         wspace=0.5, # increase horizontal space between panels
    #         size=3,
    #         show=False
    #     )
    #     plt.savefig(Path(args.outdir, 'umap_doublet.png'), bbox_inches="tight")

    with plt.rc_context():
        sc.pl.umap(
            adata,
            color=["log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
            wspace=0.5,
            ncols=2,
            show=False
        )
        plt.savefig(Path(args.outdir, 'umap_total_counts_genes_mt.png'), bbox_inches="tight")    


if __name__ == "__main__":
    sys.exit(main())
