# Set up environment, activate library components
#install.packages("ggsci")
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


# For processLTBC (lineage tracing barcodes), the following values should be set:
#    sobject - the name of the Seurat object that contains the cell barcodes (cids) that will be matched and integrated
#    lt.loc - the file location of the lineage tracing barcodes (extracted from the fastq files)
#    cid.loc - the file location of the cell id barcodes (extracted from the fastq files)
#    histogram - TRUE will trigger the function to generate and output a histogram plot of the top 40 most frequent lineage tracing barcodes
#    col.fill - the color that you would like to use for the barchart on the histogram
#    ymax - if set, this will define the upper limit of the y axis (ie, for creating side-by-side comparisons)
#    relative - TRUE will normalize cell counts to total number of cells containing barcodes
#    title - verbiage for the title of the histogram, if triggered
#    ret.list - TRUE will trigger the function to return the list of barcode frequencies, rather than the seurat object
# The following procedure will generate the matching lineage tracing and cell ID tables:
#   In Linux bash, navigate to the scRNA Counts directory containing the cellranger output .bam
#   Load samtools with > module load samtools-1.7
#   Extract all lines containing a read that matches the flanking sequences for the LT barcode
#     > samtools view possorted_genome_bam.bam | egrep "CGCA[ACGT]{14}TGGT[ACGT]{30}TGGT" > match.lt.sam
#     (You might want to make sure this is actually running when you submit the command, as
#       I've been burned before.  It takes forever even when it's working on a big file.  Open
#       a separate Putty window, SSH into the node that should be running samtools and egrep, 
#       then run top.  You should see one process for each of these programs running.)
#   Remove any lines that do not have a high-confidence cell ID assigned:
#     > egrep "CB:Z:[ACGT]{16}" match.lt.sam > match.both.sam
#   Extract the lineage tracing barcode and flanking sequences
#     > egrep -o "CBCA[ACGT]{14}TGGT[ACGT]{30}TGGT" match.both.sam > lt.fq
#   Extract the matching cell identity barcodes
#     > egrep -o "CB:Z:[ACGT]{16}" match.both.sam > pre-cid.fq
#   Clean up the cell identity barcodes (get rid of the prefix)
#     > egrep -o "[ACGT]{16}$" pre-cid.fq > cid.fq
#or 	> cut -d":" -f3 pre-cid.fq > cid.fq
#   Input lt.fq and cid.fq into the lt.loc and cid.loc variables in the function below

processLTBC <- function(sobject, lt.loc, cid.loc, histogram = F, col.fill = "#4DBBD5FF", ymax = NA, relative = F, title = "LT Barcode Frequency", ret.list = F) {
  
  # Read in the master Cellecta barcode tables for QC purposes
  bc14 <- read.table("R:/RESRoberts/Sanjana/Lineage Tracing/Cellecta/Cellecta-bc14s.txt")
  bc30 <- read.table("R:/RESRoberts/Sanjana/Lineage Tracing/Cellecta/Cellecta-bc30s.txt")
  names(bc14) <- c("label", "forward", "reverse")
  names(bc30) <- c("label", "forward", "reverse")
  bc14f <- setNames(bc14$label, bc14$forward)
  bc30f <- setNames(bc30$label, bc30$forward)
  
  # Import the lineage tracing and matching cell ID tables extracted from the fastqs
  lt.bc <- read.table(lt.loc)
  cid.bc <- read.table(cid.loc)
  # cid.bc <- paste(cid.bc[[1]],"-1", sep = "")
  bc <- data.frame(cid.bc, lt.bc)
  colnames(bc) <- c("cid", "lt")
  
  # Deduplicate redundant reads
  bc <- dplyr::distinct(bc)
  
  # Extract the bc14 and the bc30 reads from the lineage tracing sequences (identified by flanking regions)
  bc$bc14 <- substring(bc$lt, 5, 18)
  bc$bc30 <- substring(bc$lt, 23, 52)
  
  # Match extracted barcode reads against the Cellecta barcode tables
  bc$label14 <- bc14f[bc$bc14]
  bc$label30 <- bc30f[bc$bc30]
  bc$label14 <- substring(bc$label14, 6)
  bc$label30 <- substring(bc$label30, 6)
  
  # Eliminate barcodes that don't match the Cellecta barcode tables
  bc <- na.omit(bc, cols = c(label14, label30))
  
  # Concatenate the two barcodes into a single compound column
  bc$label <- paste(bc$label14, bc$label30, sep = "-")
  
  # Integrate the lineage tracing barcode into the Seurat object metadata
  bc <- setNames(as.character(bc$label), bc$cid)
  bc <- sort(bc)
  sobject$lt <- bc[sobject@assays$RNA@counts@Dimnames[[2]]]
  
  # Generate the frequency tables
  ylabel = "Number of Cells"
  bc.freq <- as.data.frame(table(sobject$lt))
  bc.freq <- bc.freq[order(-bc.freq$Freq),]
  if(isTRUE(relative)) {
    bc.freq$Freq = bc.freq$Freq/length(bc)*100
    ylabel = "Percentage of Cells"
  }
  bc.plot.data <- head(bc.freq, n=40)
  
  # Create histogram graphs if desired (default using blue color from npg from ggsci)
  if(isTRUE(histogram)) {
    p <- ggplot(bc.plot.data, aes(x = reorder(Var1, -Freq), y = Freq)) +
          geom_bar(fill = col.fill, stat = "identity") +
          ggtitle(title) +
          ylab(ylabel) +
          xlab("Top 40 Lineage Barcode") 
    #remove grid lines and grey background
    p <- p + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    p <- print(p + ylim(0, ymax))
    print(p)
  } 
  if (isTRUE(ret.list)) {
    return(bc.freq)
  } else {
    return_obj <- list("bc_plot_data" = bc.plot.data, "fig" = p, "sobject" = sobject)
    return(return_obj)
  }
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

# This function currently doesn't work when pulled as a function, however, the workflow
# outlined below in the function works well when used inline--you must copy it into your
# R to work.
plotLTBC <- function(sobject, lt, bc) {
  lt1 <- WhichCells(sobject, expression = lt == `bc`)
  sobject$old.ident <- Idents(sobject)
  Idents(sobject) <- "Others"
  Idents(sobject, cells = lt1) <- bc
  plot.umap <- DimPlot(sobject, reduction = "umap", pt.size = 1)
  print(plot.umap)
  Idents(sobject) <- sobject$old.ident
  return(sobject)
}

