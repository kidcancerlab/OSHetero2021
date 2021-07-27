library(rrrSingleCellUtils)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringr)
# library(future)
library(harmony)

set.seed(888)
# plan("multisession", workers = 10)
# options(future.globals.maxSize= 10 * 24000 * 1024^2)

# Make a list of sample names
s <- c("BC5", "BC6", "BC10", "BC11", "BC16",
       "BC17", "BC20", "BC21", "BC22")
qc <- c(18000, 25000, 25000, 30000, 70000, 
        40000, 70000, 50000, 50000)
path <- c("Conventional", "Conventional", "Conventional",
          "Conventional", "Conventional", "Chondroblastic",
          "Chondroblastic", "Intraosseous", "Chondroblastic")
type <- c("Primary", "Primary", "Lung Met", "Primary", "Primary",
          "Lung Met", "Primary", "Primary", "Primary")

# Download, file, and extract the files from GEO
if(!dir.exists("PrimaryTumor/GSE152048")) {
  tar_dir <- "PrimaryTumor/GSE152048"
  dir.create(tar_dir)
  geo_pre <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152048/suppl/GSE152048_"
  for(i in 1:length(s)){
    gse_path <- str_c(geo_pre, s[i], ".matrix.tar.gz")
    tar_file <- str_c(tar_dir, "/", s[i], ".tar.gz ")
    download.file(gse_path, destfile = tar_file, method = "auto")
    untar(tar_file, exdir = tar_dir)
    file.remove(tar_file)
  }
}

# Create a vector that contains normalized Seurat objects for all samples
raw <- c()
for(i in 1:9) {
  x <- tenx_load_qc(str_c("PrimaryTumor/GSE152048/",
                          s[i], "/"))
  x <- subset(x, subset = nCount_RNA < qc[i] & percent.mt <13)
  x$src <- s[i]
  x$type <- type[i]
  x$path <- path[i]
  x <- x %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.3) %>%
    RunUMAP(dims = 1:20)
  print(DimPlot(x, pt.size = 1, label = T) +
          coord_fixed() +
          theme(legend.position = "none") +
          ggtitle(str_c(s[i], " basic clustering")))
  raw <- c(raw, x)
}
rm(x)

# Merge into one Seurat object and integrate with harmony
comb <- merge(raw[[1]], y = c(raw[[2]], raw[[3]], raw[[4]], raw[[5]], 
                              raw[[6]], raw[[7]], raw[[8]], raw[[9]]),
              # merge.data = T,
              add.cell.ids = s,
              project = "GSE152048")

comb <- comb %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F)

comb <- RunHarmony(comb, group.by.vars = "src")
comb <- RunUMAP(comb, reduction = "harmony", dims = 1:30)
comb <- comb %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters()

# Save a stopping point - merged and harmony aligned
qs::qsave(comb, "PrimaryTumor/comb.qs")
# Start from this stopping point
comb <- qs::qread("PrimaryTumor/comb.qs")

DimPlot(comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("GSE152048 composite") +
  theme(legend.position = "none")

DimPlot(comb, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("GSE152048 composite") 

DimPlot(comb, reduction = "umap", group.by = "type") +
  coord_fixed() +
  ggtitle("GSE152048 composite") 


# Separate out the primary tumor lesions
primary <- subset(comb, subset = type == "Primary")

# Ensure that groups are not just dominated by cell cycle effects
# primary <- kill_cc(primary)
# 
# # Set a post-cc-regression stopping point
# saveRDS(primary, "primary-ccReg.rds")
# # Start from this stopping point
# primary <- readRDS("primary-ccReg.rds")

DimPlot(primary, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("GSE152048 composite - primaries") +
  theme(legend.position = "none")

DimPlot(primary, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("GSE152048 composite - primaries") 

# Calculate and show module scores from the paper
ms <- list()
ms$Osteoblastic <- c("RUNX2", "COL1A1", "CDH11", "IBSP")
ms$Chondroblastic <- c("SOX9", "ACAN", "PTH1R")
ms$Osteoclast <- c("ACP5", "CTSK", "MMP9")
ms$Myeloid <- c("CD74", "CD14", "FCGR3A")
ms$TCell <- c("CD3E", "IL7R", "CD8A", "CD4", "NKG7")
ms$NKCell <- c("NKG7", "GNLY")
ms$NKTCell <- c("NKG7", "GNLY", "CD3E")
ms$DCCell <- c("CD1C", "FCER1A", "CLEC9A", "CCR7", "CD14", "CD163")
ms$Fibroblast <- c("DCN", "COL1A1")
ms$Pericyte <- c("RGS5", "ACTA2")
ms$MSC <- c("MME", "THY1", "CXCL12", "SFRP2")
ms$Endothelial <- c("PECAM1", "VWF")
ms$Myoblast <- c("MYL1", "MYLPF")
ms$BCell <- c("MS4A1", "CD19", "JCHAIN")

mod_names <- names(ms)

for(i in 1:length(mod_names)) {
  primary <- AddModuleScore(primary, ms[i], name = names(ms[i]))
  p <- FeaturePlot(primary, features = str_c(names(ms[i]), "1"), pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  print(p)
}
