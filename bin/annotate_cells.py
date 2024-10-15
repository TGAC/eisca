#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import celltypist as ct
import pandas as pd
from matplotlib import pyplot as plt
import argparse
import sys
from pathlib import Path
import util

logger = util.get_named_logger('ANNOTATE_CELLS')


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
        "--model_file",
        metavar="MODELFILE",
        type=Path,
        help="Path to CellTypist model file.",
    )
    parser.add_argument(
        "--model",
        default='Immune_All_Low.pkl',
        help="Specify a CellTypist model name, igored if a model file specified.",
    )    
    parser.add_argument(
        "--mode",
        default='best match',
        choices=['best match', 'prob match'],
        help="The cell-type prediction mode.",
    )
    parser.add_argument(
        "--p_thres",
        type=float,
        default=0.5,
        help="Probability threshold for the multi-label classification. Ignored if mode is 'best match'.",
    )              
    parser.add_argument(
        "--no_majority_voting",
        help="Whether to disable the majority voting classifier.",
        action='store_true',
    )
    parser.add_argument(
        "--update_models",
        help="Whether to update CellTypist models.",
        action='store_true',
    )
    parser.add_argument(
        "--meta",
        default='auto',
        choices=['auto', 'sample', 'group'],
        help="Choose a metadata column as the batch for clustering",
    )                 
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.h5ad.is_file():
        logger.error(f"The given input file {args.h5ad} was not found!")
        sys.exit(2)

    if args.model_file and not args.model_file.is_file():
        logger.error(f"The given input file {args.model_file} was not found!")
        sys.exit(2)

    util.check_and_create_folder(args.outdir)
    path_annotation = Path(args.outdir)
    util.check_and_create_folder(path_annotation)

    adata = sc.read_h5ad(args.h5ad)

    # update CellTypist cell-type models
    if args.update_models:
        ct.models.models_description(on_the_fly = True)
    

    # perform cell-type annotation
    # model = ct.models.Model.load(model = 'Immune_All_Low.pkl')
    model = args.model_file if args.model_file else args.model
    majority_voting = False if args.no_majority_voting else True
    predictions = ct.annotate(
        adata, 
        model=model, 
        majority_voting=majority_voting, 
        mode=args.mode,
        p_thres=args.p_thres,
    )
    adata = predictions.to_adata()


    # save the AnnData into a h5ad file
    adata.write_h5ad(Path(path_annotation, 'adata_annotation.h5ad'))


    # plot UMAPs to visualise the predicted cell-types
    if args.meta == 'auto':
        batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
    else:
        batch = args.meta

    label_type = 'majority_voting' if hasattr(adata.obs, 'majority_voting') else 'predicted_labels'
    for sid in sorted(adata.obs[batch].unique()):
        adata_s = adata[adata.obs[batch]==sid]   
        path_annotation_s = Path(path_annotation, f"{batch}_{sid}")
        util.check_and_create_folder(path_annotation_s)
        for res in args.resolutions:
            with plt.rc_context():
                sc.pl.umap(
                    adata_s,
                    color=label_type,
                    # legend_loc="on data",
                    show=False
                )
                plt.savefig(Path(path_annotation_s, f"umap_cell_type.png"), bbox_inches="tight")    

            with plt.rc_context():
                sc.pl.umap(
                    adata_s,
                    color='conf_score',
                    # legend_loc="on data",
                    show=False
                )
                plt.savefig(Path(path_annotation_s, f"umap_conf_score.png"), bbox_inches="tight")          


    # stacked proportion bar plot to compare between batches     
    n_cluster = len(adata.obs[label_type].unique())+1
    ncol = min((n_cluster//20 + min(n_cluster%20, 1)), 3)
    with plt.rc_context():
        prop = pd.crosstab(adata.obs[label_type], adata.obs[batch], normalize='columns').T.plot(kind='bar', stacked=True)
        prop.legend(bbox_to_anchor=(1.4+ncol*0.17, 1.02), loc='upper right', ncol=ncol)
        plt.savefig(Path(path_annotation, f"prop_{label_type}.png"), bbox_inches="tight")    



if __name__ == "__main__":
    sys.exit(main())
