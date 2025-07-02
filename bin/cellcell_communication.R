#!/usr/bin/env Rscript

library(CellChat)
library(patchwork)
library(Seurat)
library(Matrix)
library(jsonlite)
options(stringsAsFactors = FALSE)



library(argparse)

parser <- ArgumentParser()
parser$add_argument("--outdir", help="Output directory")
parser$add_argument("--count", help="Input count matrix mtx file")
parser$add_argument("--metadata", help="Input metadata csv file")
parser$add_argument("--gids", help="Input gene ID list csv file")
parser$add_argument("--cids", help="Input cell ID list csv file")
parser$add_argument("--group", default="majority_voting", help="Specify the column for grouping the cells")
parser$add_argument("--normalize", action="store_true", help="Indicates whether to normalize the counts")
parser$add_argument("--db", default="human",  choices=c('human', 'mouse'), help="Specify the species of CellChatDB.")
parser$add_argument("--dbc", help="The categories of CellChatDB, e.g. Secreted Signaling")
# parser$add_argument("--dbv", help="The version of CellChatDB, e.g. v1")
parser$add_argument("--dbenps", action="store_true", help="Use all CellChatDB excepting 'Non-protein Signaling'")
parser$add_argument("--threads", type="integer", default=4, help="Number of threads for parallel runs")
parser$add_argument("--mean_method", default="triMean",  choices=c('triMean', 'truncatedMean'), 
        help="Specify the method for calculating the average gene expression per cell group")
parser$add_argument("--mincells", type="integer", default=10, help="The minimum number of cells required in each cell group")
parser$add_argument("--meta", default="auto", choices=c('auto', 'sample', 'group'), help="Specify a metadata column to define separate subsets of cells for analysis")
parser$add_argument("--pdf", action="store_true", help="Whether to generate figure files in PDF format")
# parser$add_argument("--plotcoef", type="float", default=1.0, help="plot size equals the coef times default plot size")

args <- parser$parse_args()


dir.create(args$outdir, showWarnings = FALSE)

# create a cellchat object from counts and metadata
meta <- read.csv(args$metadata, header = TRUE, row.names = 1)
counts <- t(readMM(args$count))
genes <- readLines(args$gids)
cells <- readLines(args$cids)
rownames(counts) <- genes
colnames(counts) <- cells


# normalize the count data if input data is raw counts
if(args$normalize){
    library.size <- Matrix::colSums(counts)
    counts <- as(log1p(Matrix::t(Matrix::t(counts)/library.size) * 10000), "CsparseMatrix")
}

# set a column to define separate batches for analysis
batch <- "sample"
if (args$meta == "auto") {
  if ("group" %in% colnames(meta)) {
    batch <- "group"
  } else if ("plate" %in% colnames(meta)) {
    batch <- "plate"
  }
} else {
  batch <- args$meta
}     


# Set the ligand-receptor interaction database
if(args$db == "human"){
    CellChatDB <- CellChatDB.human
}else if(args$db == "mouse"){
    CellChatDB <- CellChatDB.mouse
}
CellChatDB.use <- CellChatDB
if(!is.null(args$dbc) && nchar(args$dbc) > 0){
    CellChatDB.use <- subsetDB(CellChatDB, search = args$dbc, key = "annotation")
}else if(args$dbenps){
    CellChatDB.use <- subsetDB(CellChatDB)
}


for(sid in unique(meta[[batch]])){
    sub_idx <- which(meta[[batch]] == sid)
    counts_s <- counts[, sub_idx]
    meta_s <- meta[sub_idx,]
    path_outdir_s <- paste0(args$outdir, '/', batch, '_', sid)
    dir.create(path_outdir_s, showWarnings = FALSE)

    cellchat <- createCellChat(object = counts_s, meta = meta_s, group.by = args$group)

    if(!is.null(CellChatDB.use)){cellchat@DB <- CellChatDB.use}

    # Preprocessing the expression data
    cellchat <- subsetData(cellchat)
    future::plan("multisession", workers = args$threads)
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # Inference of cell-cell communication network
    cellchat <- computeCommunProb(cellchat, type = args$mean_method)
    cellchat <- filterCommunication(cellchat, min.cells = args$mincells)
    # Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)

    # Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)

    # Compute the network centrality scores
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    groupSize <- as.numeric(table(cellchat@idents))
    png(file=paste0(path_outdir_s, "/aggregated_network_all.png"), width=8,height=8, units="in", res=100)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    dev.off()
    png(file=paste0(path_outdir_s, "/aggregated_network_all_weights.png"), width=8,height=8, units="in", res=100)
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()

    mat <- cellchat@net$weight
    grpnames <- rownames(mat)
    Ngrp <- nrow(mat)
    # Identify rows and columns with any non-zero interaction weights
    row_has_interaction <- rowSums(mat) > 0
    col_has_interaction <- colSums(mat) > 0
    active_groups_logical <- row_has_interaction & col_has_interaction
    actgrp_idx <- as.numeric(which(active_groups_logical))

    png(file=paste0(path_outdir_s, "/aggregated_network_groups.png"),width=8,height=3*ceiling(Ngrp/3), units = "in", res=100*ceiling(Ngrp/3))
    par(mfrow = c(ceiling(Ngrp/3),3), xpd=TRUE)
    for (i in 1:Ngrp) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()

    if(args$pdf){
        pdf(file=paste0(path_outdir_s, "/aggregated_network_all.pdf"), width=8,height=8)
        netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
        dev.off()
        pdf(file=paste0(path_outdir_s, "/aggregated_network_all_weights.pdf"), width=8,height=8)
        netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
        dev.off()

        pdf(file=paste0(path_outdir_s, "/aggregated_network_groups.pdf"),width=8,height=3*ceiling(Ngrp/3))
        par(mfrow = c(ceiling(Ngrp/3),3), xpd=TRUE)
        for (i in 1:Ngrp) {
        mat2 <- matrix(0, nrow = Ngrp, ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
        }
        dev.off()
    }

    # Visualization of cell-cell communication network for individual signalling pathways
    for(pathway in cellchat@netP$pathways){
        png(file=paste0(path_outdir_s, "/pathway_network_circle_", pathway, '.png'), width=8,height=8, units="in", res=100)
        par(mar=c(1, 3, 1, 3), cex=1.5)
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle", signaling.name = pathway)
        dev.off()
        png(file=paste0(path_outdir_s, "/pathway_network_chord_", pathway, '.png'), width=8,height=8, units="in", res=100)
        netVisual_chord_cell(cellchat, signaling = pathway, lab.cex=1)
        dev.off()
        png(file=paste0(path_outdir_s, "/pathway_network_heatmap_", pathway, '.png'), width=1.5*Ngrp,height=1*Ngrp, units="in", res=100)
        print(netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds", font.size = 12, font.size.title = 15))
        dev.off()
        if(args$pdf){
            pdf(file=paste0(path_outdir_s, "/pathway_network_circle_", pathway, '.pdf'), width=8,height=8)
            par(mar=c(1, 3, 1, 3), cex=1.5)
            netVisual_aggregate(cellchat, signaling = pathway, layout = "circle", signaling.name = pathway)
            dev.off()
            pdf(file=paste0(path_outdir_s, "/pathway_network_chord_", pathway, '.pdf'), width=8,height=8)
            netVisual_chord_cell(cellchat, signaling = pathway, lab.cex=1)
            dev.off()
            pdf(file=paste0(path_outdir_s, "/pathway_network_heatmap_", pathway, '.pdf'), width=1.5*Ngrp,height=1*Ngrp)
            print(netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds", font.size = 12, font.size.title = 15))
            dev.off()     
        }

        # Compute the contribution of each ligand-receptor pair to the overall signaling pathway
        pairLR <- extractEnrichedLR(cellchat, signaling = pathway, geneLR.return = FALSE)
        png(file=paste0(path_outdir_s, "/pathway_network_contribution_", pathway, '.png'), width=6,height=1.5*nrow(pairLR), units="in", res=100)
        print(netAnalysis_contribution(cellchat, signaling = pathway, font.size = 12, font.size.title = 15))
        dev.off()
        if(args$pdf){
            pdf(file=paste0(path_outdir_s, "/pathway_network_contribution_", pathway, '.pdf'), width=6,height=1.5*nrow(pairLR))
            print(netAnalysis_contribution(cellchat, signaling = pathway, font.size = 12, font.size.title = 15))
            dev.off()
        }

        # visualize cell-cell communication mediated by a single ligand-receptor pair
        for(i in 1:nrow(pairLR)){
            LR <- pairLR[i,]
            png(file=paste0(path_outdir_s, "/pathway_network_LR_circle_", pathway, '.png'), width=5,height=5, units="in", res=200)
            netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR, layout = "circle")
            dev.off()
            png(file=paste0(path_outdir_s, "/pathway_network_LR_chord_", pathway, '.png'), width=6,height=6, units="in", res=200)
            netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR, layout = "chord")
            dev.off()
            if(args$pdf){
                pdf(file=paste0(path_outdir_s, "/pathway_network_LR_circle_", pathway, '.pdf'), width=5,height=5)
                netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR, layout = "circle")
                dev.off()
                pdf(file=paste0(path_outdir_s, "/pathway_network_LR_chord_", pathway, '.pdf'), width=6,height=6)
                netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR, layout = "chord")
                dev.off()
            }
        }

        # Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
        for(i in actgrp_idx){
            grpname = gsub("[/ ]", "_", grpnames[i])
            try({
                gg <- netVisual_bubble(cellchat, sources.use = i, remove.isolate = FALSE)
                ggsave(filename=paste0(path_outdir_s, "/cellcell_LR_bubble_", grpname, '.png'), plot=gg, units = 'in', dpi = 200)
                if(args$pdf){
                    ggsave(filename=paste0(path_outdir_s, "/cellcell_LR_bubble_", grpname, '.pdf'), plot=gg)
                }
            },silent = TRUE)
        }
        for(i in actgrp_idx){
            grpname = gsub("[/ ]", "_", grpnames[i])
            png(file=paste0(path_outdir_s, "/cellcell_LR_chord_", grpname, '.png'), width=10,height=8, units="in", res=100)
            netVisual_chord_gene(cellchat, sources.use = i, lab.cex = 0.5,legend.pos.y = 30)
            dev.off()
            if(args$pdf){
                pdf(file=paste0(path_outdir_s, "/cellcell_LR_chord_", grpname, '.pdf'), width=10,height=8)
                netVisual_chord_gene(cellchat, sources.use = i, lab.cex = 0.5,legend.pos.y = 30)
                dev.off()
            }
        }

        # Plot the signaling gene expression distribution using violin
        png(file=paste0(path_outdir_s, "/pathway_genes_violin_", pathway, '.png'), width=1*Ngrp,height=8, units="in", res=100)
        print(plotGeneExpression(cellchat, signaling = pathway, enriched.only = TRUE))
        dev.off()
        if(args$pdf){
            pdf(file=paste0(path_outdir_s, "/pathway_genes_violin_", pathway, '.pdf'), width=1*Ngrp,height=8)
            print(plotGeneExpression(cellchat, signaling = pathway, enriched.only = TRUE))
            dev.off()
        }

        # visualize the network centrality scores
        png(file=paste0(path_outdir_s, "/pathway_network_centrality_", pathway, '.png'), width=1*Ngrp,height=4, units="in", res=100)
        netAnalysis_signalingRole_network(cellchat, signaling = pathway, width=1.5*Ngrp,height=4, font.size = 10)
        dev.off()
        if(args$pdf){
            pdf(file=paste0(path_outdir_s, "/pathway_network_centrality_", pathway, '.pdf'), width=1*Ngrp,height=4)
            netAnalysis_signalingRole_network(cellchat, signaling = pathway, width=1.5*Ngrp,height=4, font.size = 10)
            dev.off()
        }
    }

    # Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
    png(file=paste0(path_outdir_s, "/heatmap_signaling_patterns.png"), width=(1*Ngrp+3),height=1*Ngrp, units="in", res=100)
    print(ht1 + ht2)
    dev.off()
    if(args$pdf){
        pdf(file=paste0(path_outdir_s, "/heatmap_signaling_patterns.pdf"), width=(1*Ngrp+3),height=1*Ngrp)
        print(ht1 + ht2)
        dev.off()
    }


    #  save dataframe consisting of all the inferred cell-cell communications at the level of ligands/receptors
    saveRDS(subsetCommunication(cellchat), file = paste0(path_outdir_s, "/inferred_cellcell_comm.rds"))
    # Save the CellChat object
    saveRDS(cellchat, file = paste0(path_outdir_s, "/cellchat.rds"))


}


# save analysis parameters into a json file
params <- list()
params[["--count"]] <- as.character(args$count)
params[["--metadata"]] <- as.character(args$metadata)
params[["--gids"]] <- as.character(args$gids)
params[["--cids"]] <- as.character(args$cids)
params[["--group"]] <- as.character(args$group)
if(args$normalize){params[["--normalize"]] <- as.character(args$normalize)}
params[["--db"]] <- as.character(args$db)
if (!is.null(args$dbc)) {params[["--dbc"]] <- as.character(args$dbc)}
if(args$dbenps){params[["--dbenps"]] <- as.character(args$dbenps)}
params[["--mean_method"]] <- as.character(args$mean_method)
params[["--mincells"]] <- as.character(args$mincells)
params[["--threads"]] <- as.character(args$threads)
write_json(params, file.path(args$outdir, "parameters.json"), pretty = TRUE, auto_unbox = TRUE)



# yaml::write_yaml(
# list(
#     'MTX_TO_SEURAT'=list(
#         'Seurat' = paste(packageVersion('Seurat'), collapse='.')
#     )
# ),
# "versions.yml"
# )
