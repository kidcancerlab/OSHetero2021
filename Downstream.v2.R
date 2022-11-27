# Set up environment, activate library components
library("ggsci")
library("cowplot")
library("dplyr")
library("Matrix")
library("reticulate")
library("Seurat")
library("reshape2")
library("ggplot2")

# Define functions
# enricher - https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
# WRITE ANNOTATION 
# Function for msigdb and SingleR annotation

DGEA <- function(data) {
  # Set the Human Primary Cell Atlas as a reference
  hpca.se <- celldex::HumanPrimaryCellAtlasData()
  # Preparing clusterProfiler to perform pathway enrichment on msigdb signatures
  m_t2g.c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>%
    dplyr::select(gs_name, human_gene_symbol)
  m_t2n.c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>%
    dplyr::select(gs_id, gs_name)
  m_t2g.c6 <- msigdbr(species = "Homo sapiens", category = "C6") %>%
    dplyr::select(gs_name, human_gene_symbol)
  m_t2g.h <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, human_gene_symbol)
  m_t2n.h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_id, gs_name)
  m_t2g <- rbind(m_t2g.c2, m_t2g.c6)

  # getting log normalized data for specific cluster
  clust.ids <- sort(unique(data@active.ident))
  new.cluster.ids <- rep(NA, length(clust.ids))

  # store top 30 pathway enrichment analysis
  em <- NULL

  # msigdb signature to use
  msig.gene.set <- m_t2g.h
  msig.name <- m_t2n.h

  for (i in 1:length(clust.ids)) {
    clust <- GetAssayData(subset(x = data, 
      idents = clust.ids[i]),
      slot = "data")
    label <- rep(clust.ids[i], dim(clust)[2])

    # getting common genes
    common <- intersect(rownames(clust), rownames(hpca.se))

    # use only differential markers
    cluster.markers <- FindMarkers(data,
      ident.1 = clust.ids[i],
      logfc.threshold = 0.25,
      only.pos = TRUE, 
      test.type = "DESeq2")

    common <- intersect(common, rownames(cluster.markers))
    hpca.se.common <- hpca.se[common, ]

    tmp <- enricher(rownames(cluster.markers), 
      TERM2GENE = msig.gene.set,
      TERM2NAME = msig.name)

    em[[i]] <- tmp@result[, c("ID", "p.adjust")]
  }

  # Create a heatmap of enrichment
  library(pheatmap)

  # get top 10 enrichments
  em.table.top10 <- lapply(em, function(x) x[1:10, ])

  # create dataframe for heatmap
  em.hm <- NULL
  em.hm <- data.frame(
    gene_set = unique(unlist(lapply(em.table.top10, function(x) rownames(x)))))

  for (i in 1:length(em.hm$gene_set)) {
    for (j in 1:length(clust.ids)) {
      em.hm[i, j+1] <- em[[j]]$p.adjust[match(em.hm$gene_set[i], em[[j]]$ID)]
    }
  }

  rownames(em.hm) <- em.hm[, 1]
  em.hm <- em.hm[, -1]
  em.hm[is.na(em.hm)] <- 1
  colnames(em.hm) <- as.character(clust.ids)

  return(em.hm)
}

# Determine average Silhouette scores for each specified resolution.
# (calcualted using the silhouette() function (package clustter))
# References: https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html,
# https://github.com/satijalab/seurat/issues/1985
# The following variables can be defined:
# ' @param sobject A Seurat object containing all of the cells for analysis (required)
# ' @param res A character vector of resolutions to investigate (required)
# This function returns a list containing the following objects: 
# - input Seurat object [1], 
# - list of calculated silhouette scores for each specified resolution [2],
# - list of specified resolution as found in sobject metadata [3] and 
# - data frame of means of silhouette scores calculated for each specified resolution

# Example:
# sobject.nRes <- nRes(sobject, res = seq(from = 0.1, to = 0.3, by = 0.1))

nRes <- function(sobject, res) {

  sobject <- FindClusters(sobject, resolution = res)
  resolution_list <- paste("RNA_snn_res.", res, sep = "")
  resolution_list <- sort(resolution_list)

  library(cluster, quietly = TRUE)

  dist.matrix <- dist(x = Embeddings(object = sobject[["pca"]])[, 1:20])
  values <- list()
  silscore <- list()

  for (resolution in resolution_list) {
    clusters <- sobject@meta.data[[resolution]]
    sil <- silhouette(
      x = as.numeric(x = as.factor(x = clusters)),
      dist = dist.matrix)
    values[[resolution]] <- sil[, 3]
    silscore[[resolution]] <- sil
  }

  a <- names(values)
  means <- data.frame(x = a, y = 0)
  for (i in 1:length(values)) {
    means[i, 2] <- mean(values[[i]])
  }

  bestc <- means[which.max(means$y),]
  bestc <- bestc[,1]

  bplot <- boxplot(values, plot = TRUE,
    main = (paste(bestc, "is the resolution with highest Sil score")),
    xlab = "Resolution",
    ylab = "Sil Score",
    col = "gold")

  return(list(sobject, silscore, resolution_list, means))
}

# Silhouette plot to visualize silhouette score distribution of cells in each cluster.
# The following variables can be defined:
# ' @param sobject.nRes Output list from 'nRes' function (required; refer to 
# documentation above)
# ' @param res Desired resolution to use to generate silhouette plot (required)
# RStudio does not plot silhouette plot properly.
# This function does not return anything to the R interpreter instead plots the
# silhouette plot in a separate plot device
# Example:
# plot <- pSil(sobject.nRes , 0.15)

library(RColorBrewer)
n <- 60
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(
  brewer.pal,
  qual_col_pals$maxcolors,
  rownames(qual_col_pals)))

pSil <- function(sobject.nRes, res) {
  sobject <- sobject.nRes[[1]]
  silscore <- sobject.nRes[[2]]
  resolution_list <- sobject.nRes[[3]]

  sobject <- FindClusters(sobject, resolution = res)
  res <- paste("RNA_snn_res.", res, sep = "")
  k <- length(levels(sobject@meta.data[[res]]))
  col <- pal_npg("nrc")(n)
  resolution <- res
  n <- match(res, resolution_list)

  p <- plot(silscore[[n]],
    main = paste("res = ", resolution),
    do.n.k = FALSE,
    col = col[1:k])
  print(p)
}