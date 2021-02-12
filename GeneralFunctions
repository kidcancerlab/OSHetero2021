# Set up environment, activate library components
library(ggsci)
library("cowplot")
library("dplyr")
library("Matrix")
library("reticulate")
library("Seurat")
library("reshape2")
library("ggplot2")
# library("harmony")
# library("future")

# plan(strategy = "multicore", workers = 3)


## Define functions

# For tenXLoadQC, must define the path to the 10x filtered_feature_bc_matrix folder (path10x) and the type of extraction you want to do (spec).
# Valid spec values are:
#    human - use for an alignment made to a human genome
#    mouse - use for an alignment made to a murine genome
#    mixHuman - use for an alignment made to a mixed genome, return human dataset
#    mixMouse - use for an alignment made to a mixed genome, return murine dataset
tenXLoadQC <- function(path10x, spec) {
  raw.data <- Read10X(path10x)
  
  if(spec == "human") {
    seurat <- CreateSeuratObject(raw.data, min.cells = 5, min.features = 800)
    seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
  }
  else if(spec == "mouse") {
    seurat <- CreateSeuratObject(raw.data, min.cells = 5, min.features = 800)
    seurat <- PercentageFeatureSet(seurat, pattern = "^mt-", col.name = "percent.mt")
  }
  else if(spec == "mixHuman") {
    raw.data <- raw.data[grep(pattern = "^hg19", raw.data@Dimnames[[1]]),]
    raw.data@Dimnames[[1]] <- substring(raw.data@Dimnames[[1]], 6)
    seurat <- CreateSeuratObject(raw.data, min.cells = 5, min.features = 800)
    seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
  }
  else if(spec == "mixMouse") {
    raw.data <- raw.data[grep(pattern = "^mm10", raw.data@Dimnames[[1]]),]
    raw.data@Dimnames[[1]] <- substring(raw.data@Dimnames[[1]], 6)
    seurat <- CreateSeuratObject(raw.data, min.cells = 5, min.features = 800)
    seurat <- PercentageFeatureSet(seurat, pattern = "^mt-", col.name = "percent.mt")
  }
  vplot <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(vplot)
  return(seurat)
}


# The function PCAtoPlot is designed to take a QCed and trimmed Seurat object, perform a normalization/
#   scaling transformation, regress out the mitochondrial gene percentage, find the principal components,
#   jackstraw if needed for choosing significant PCs, perform dimensional reduction (UMAP/tSNE), generate
#   an initial plot, and output the Seurat object with the PCA and dimensional reduction implemented
# For process PCAtoUMAP, the following options can be set:
    # sobject - the name of the Seurat object that will be processed
    # find.pcs - number of principal components to generate in the inital PCA
    # res - resolution for input to FindClusters
    # use.pcs - number of principal components to use in the dimensional reduction (if already known), setting
    #   this paramater will skip the jackplot process, saving time
    # method - type of dimensional reduction to use, currently supports either umap or tsne (umap by default)

PCAtoPlot <- function(sobject, find.pcs = 20, res = 0.5, use.pcs = 0, method = "umap") {
  sobject <- SCTransform(sobject, vars.to.regress = "percent.mt")
  sobject <- RunPCA(sobject, npcs = find.pcs)
  if(use.pcs == 0) {
    sobject <- JackStraw(sobject, num.replicate = find.pcs)
    sobject <- ScoreJackStraw(sobject, dims = 1:find.pcs)
    plot.jack <- JackStrawPlot(sobject, dims = 1:find.pcs)
    print(plot.jack)
    use.pcs <- readline(prompt = "How many PCA components do you want to use? ")
  }
  if(method == "umap") {
    sobject <- RunUMAP(sobject, reduction = "pca", dims = 1:use.pcs)
    sobject <- FindNeighbors(sobject, reduction = "pca", dims = 1:use.pcs)
    sobject <- FindClusters(sobject, resolution = res)
    plot.umap <- DimPlot(sobject, reduction = "umap", label = TRUE, pt.size = 1)
    print(plot.umap)
  } else {
    if(method == "tsne") {
      sobject <- RunTSNE(sobject, reduction = "pca", dims = 1:use.pcs)
      sobject <- FindNeighbors(sobject, reduction = "pca", dims = 1:use.pcs)
      sobject <- FindClusters(sobject, resolution = res)
      plot.tsne <- DimPlot(sobject, reduction = "tsne", label = TRUE, pt.size = 1)
      print(plot.tsne)
    } else {
      print("Method is not a method recognized in this function.")
    }
  }
  return(sobject)
}



# The killCC function will identify cell cycle components within a dataset.  After an initial scoring 
# using the Seurat CellCycleScoring function, the user will be shown a dimensional reduction plot with cells
# labeled by cell cycle.  If indicated, the user can then trigger a process to regress out the effects of
# cell cycle within the dataset.  The function will then proceed to re-do the PCA and jackstraw if needed,
# then show a dimensional reduction plot post-regression and retun the corrected Seurat object.
# Input must be a Seurat object that already has PCA and dimensional reduction data (umap or tsne) attached.
# For process killCC, the following options can be set:
  #   sobject - the name of the Seurat object to be processed
  #   cc.regress - if set to Y, the process with run without user input and will automatically proceed to
  #     cell cycle regression
  #   find.pcs - the number of principal componets to generate in the re-do PCA post-CC regression
  #   use.pcs - the number of principal components to use in the post-regression dimensional reduction
  #     (if already known), setting this parameter will skip the jackplot process
  #   use.res - resolution to input to FindClusters
  #   method - type of dimensional reduction to use, currently supports either umap or tsne (umap by default)

killCC <- function(sobject, cc.regress = "N", find.pcs = 20, use.pcs = 0, use.res = 0.5, method = "umap") {
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  sobject <- CellCycleScoring(sobject, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  if(method == "umap") {
    plot.cc <- DimPlot(sobject, reduction = "umap", label = TRUE, pt.size = 1)
  } else {
      if(method == "tsne") {
        plot.cc <- DimPlot(sobject, reduction = "tsne", label = TRUE, pt.size = 1)
      } else {
        print ("You need to set a method supported by this function")
        }
    }
  print(plot.cc)
    
  if(cc.regress != "Y") {
    cc.regress <- readline(prompt = "Proceed with regression of cell cycle-dependent genes (Y/N)? ")
  }
  if(cc.regress == "Y") {
    sobject <- ScaleData(sobject, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(x = sobject))
    sobject <- RunPCA(sobject, npcs = find.pcs)
    if(use.pcs == 0) {
      sobject <- JackStraw(sobject, num.replicate = 50)
      sobject <- ScoreJackStraw(sobject, dims = 1:find.pcs)
      plot.jack <- JackStrawPlot(sobject, dims = 1:find.pcs)
      print(plot.jack)
      use.pcs <- readline(prompt = "How many PCA components do you want to use? ")
    }
    if(method == "umap") {
      sobject <- RunUMAP(sobject, reduction = "pca", dims = 1:use.pcs)
      sobject <- FindNeighbors(sobject, reduction = "pca", dims = 1:use.pcs)
      sobject <- FindClusters(sobject, resolution = use.res)
      plot.umap <- DimPlot(sobject, reduction = "umap", label = TRUE, pt.size = 1, group.by = "Phase")
      print(plot.umap)
      plot.umap <- DimPlot(sobject, reduction = "umap", label = TRUE, pt.size = 1)
      print(plot.umap)
    } else {
      if(method == "tsne") {
        sobject <- RunTSNE(sobject, reduction = "pca", dims = 1:use.pcs)
        sobject <- FindNeighbors(sobject, reduction = "pca", dims = 1:use.pcs)
        sobject <- FindClusters(sobject, resolution = use.res)
        plot.tsne <- DimPlot(sobject, reduction = "tsne", label = TRUE, pt.size = 1, group.by = "Phase")
        print(plot.tsne)
        plot.tsne <- DimPlot(sobject, reduction = "tsne", label = TRUE, pt.size = 1)
        print(plot.tsne)
      }
    }
  } else {
    print("No CC regression performed.")
  }
  return(sobject)
}

