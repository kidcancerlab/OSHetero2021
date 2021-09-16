library(rrrSingleCellUtils)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringr)
# library(future)
library(harmony)
library(devtools)

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
# 
# tar_dir <- "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/tar_dir"
#   
#   geo_pre <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152048/suppl/GSE152048_"
#   for(i in 1:length(s)){
#     gse_path <- str_c(geo_pre, s[i], ".matrix.tar.gz")
#     tar_file <- str_c(tar_dir, "/", s[i], ".tar.gz ")
#     download.file(gse_path, destfile = tar_file, method = "auto")
#     untar(tar_file, exdir = tar_dir)
#     file.remove(tar_file)
#   }


# Create a vector that contains normalized Seurat objects for all samples
raw <- c()
for(i in 1:9) {
  x <- tenx_load_qc(str_c("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/tar_dir/",
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

# Save a stopping point - individual objects
save(raw, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/tar_dir/raw.RData")


# Merge into one Seurat object and integrate with harmony
comb <- merge(raw[[1]], y = c(raw[[2]], raw[[3]], raw[[4]], raw[[5]], 
                              raw[[6]], raw[[7]], raw[[8]], raw[[9]]),
              # merge.data = T,
              add.cell.ids = s,
              project = "GSE152048")
rm(raw)

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
save(comb, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/tar_dir/comb.RData")

# Start from this stopping point
# comb <- qs::qread("PrimaryTumor/comb.qs")

DimPlot(comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("GSE152048 composite") +
  theme(legend.position = "none")

DimPlot(comb, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("GSE152048 composite") 

DimPlot(comb, reduction = "umap", split.by = "type") +
  coord_fixed() +
  ggtitle("GSE152048 composite") 


# Separate out the primary tumor lesions
primary <- subset(comb, subset = type == "Primary")

# Ensure that groups are not just dominated by cell cycle effects
#primary <- kill_cc(primary)
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
#Added min.cutoff to isolate for cells most module like
for(i in 1:length(mod_names)) {
  comb <- AddModuleScore(comb, ms[i], name = names(ms[i]))
  p <- FeaturePlot(comb, features = str_c(names(ms[i]), "1"), min.cutoff = 2, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  pdf(paste0("FeaturePlot ",names(ms[i]), ".pdf"))
  print(p)
  dev.off()
}

#These markers don't show up in cell line and PDX derived tumors

#look for proliferation
comb<-kill_cc(comb)
DimPlot(comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("G1/S/G2M cell cycle analysis - comb tumor") +
  theme(legend.position = "none")


#Select Clusters which are tumors. Notice how one group is G2M only. 
#Maybe possible that this group is tumor and due to 0inflation there are transcripts which are not active. 
#Could be cycling back to other tumor states
Idents(comb) <- comb$seurat_clusters
DimPlot(comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("GSE152048 composite clusters - primaries")
# saveRDS(comb,'comb.rds')
# comb<-readRDS('comb.rds')

#Re-cluster after subsetting to tumor cells
tumor<- subset(comb, idents=c(1,3,28,17,15,6,5,0))
tumor<-tumor%>%
  NormalizeData()%>%
  ScaleData()%>%
  FindVariableFeatures()
tumor<-tumor%>%
  RunPCA(features=VariableFeatures(tumor))%>%

tumor <- RunHarmony(tumor, group.by.vars = "src")
tumor <- RunUMAP(tumor, reduction = "harmony", dims = 1:30)
tumor <- tumor %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters()

#Plot results
DimPlot(tumor, pt.size=1, reduction="umap", label=T, group.by="src")+
  coord_fixed()+
  ggtitle("Reclustering of Tumor Cells")
#Plot results
DimPlot(tumor, pt.size=1, reduction="umap", label=T)+
  coord_fixed()+
  ggtitle("Reclustering of Tumor Cells")
# #kill_cc adds metadata but interferes with analyses
tumor<-kill_cc(tumor)

Idents(tumor) <- tumor$seurat_clusters
DimPlot(tumor, pt.size=1, reduction="umap", label=T, group.by = "src")+
  coord_fixed()+
  ggtitle("Reclustering of Tumor Cells")

#Added min.cutoff to isolate for cells most module like
for(i in 1:length(mod_names)) {
  tumor <- AddModuleScore(tumor, ms[i], name = names(ms[i]))
  p <- FeaturePlot(tumor, features = str_c(names(ms[i]), "1"), min.cutoff = 2, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  pdf(paste0("FeaturePlot ",names(ms[i]), ".pdf"))
  print(p)
  dev.off()
}

# Save a stopping point - subset tumor cells
save(tumor, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/tar_dir/tumor.RData")

# #Total number of tumor cells identified much more than that reported in the paper
# #Re-Subset data based on expression of chondroblastic and osteoblastic markers

comb$Osteoblastic2 <- comb$Osteoblastic1>2
comb$Chondroblastic2 <- comb$Chondroblastic1>2

tumor_a<-subset(comb, subset = Osteoblastic2 == "TRUE")
tumor_b<-subset(comb, subset = Chondroblastic2 == "TRUE")

names_a <- names(Idents(tumor_a))
names_b <- names(Idents(tumor_b))
names <- c(names_a, names_b)
names <- unique(names) ##Looks like there is no overlap

tumor_2 <- subset(comb, cells = names)

#Cutoff of 2 is too stringent, reducing to 1
#Added min.cutoff to isolate for cells most module like
for(i in 1:2) {
  # comb <- AddModuleScore(comb, ms[i], name = names(ms[i]))
  p <- FeaturePlot(comb, features = str_c(names(ms[i]), "1"), min.cutoff = 1.5, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))+coord_fixed()
  pdf(paste0("FeaturePlot ",names(ms[i]), ".pdf"))
    print(p)
    dev.off()
}

#Recluster new tumor
tumor_2<-tumor_2%>%
  NormalizeData()%>%
  ScaleData()%>%
  FindVariableFeatures()
tumor_2<-tumor_2%>%
  RunPCA(features=VariableFeatures(tumor_2))%>%
  RunUMAP(dims=1:20)%>%
  FindNeighbors()%>%
  FindClusters(resolution=0.3)

#Plot results
DimPlot(tumor_2, pt.size=.1, reduction="umap", label=T, group.by="src")+
  coord_fixed()+
  ggtitle("Reclustering of tumor_2 Cells")
#Plot results
DimPlot(tumor_2, pt.size=.1, reduction="umap", label=T)+
  coord_fixed()+
  ggtitle("Reclustering of tumor_2 Cells")
DimPlot(tumor_2, pt.size=.1, reduction="umap", label=T,group.by="seurat_clusters")+
  coord_fixed()+
  ggtitle("Reclustering of tumor_2 Cells")
#kill_cc adds metadata but interferes with analyses
tumor_2<-kill_cc(tumor_2)
Idents(tumor_2)<-tumor$seurat_clusters
for(i in 1:length(mod_names)) {
  tumor_2 <- AddModuleScore(tumor_2, ms[i], name = names(ms[i]))
  p <- FeaturePlot(tumor_2, features = str_c(names(ms[i]), "1"), pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  print(p)
}
saveRDS(tumor_2,'tumor_2.rds')
tumor_2<-readRDS('tumor_2.rds')

#Now that sure of tumor cells. Take the tumor subsets from "tumor" and combine across data with harmony
tumor_comb<-subset(tumor, idents=c(1,2,3,4,5,6,9,10,11,12,13,14))
tumor_comb <- tumor_comb %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F)

tumor_comb <- RunHarmony(tumor_comb, group.by.vars = "src")
tumor_comb <- RunUMAP(tumor_comb, reduction = "harmony", dims = 1:30)
tumor_comb <- tumor_comb %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution=0.3)

# Save a stopping point - merged and harmony aligned
library(qs)
# qs::qsave(tumor_comb, "PrimaryTumor/tumor_comb.qs")


# Start from this stopping point
tumor_comb <- qread("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/tar_dir/tumor_comb.qs")

DimPlot(tumor_comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("Tumor Clusters post-integration") +
  theme(legend.position = "none")

DimPlot(tumor_comb, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("Tumor by sample post-integration") 

tumor_comb<-kill_cc(tumor_comb)

#Run cell type module as before
Idents(tumor_comb)<-tumor_comb$seurat_clusters
for(i in 1:length(mod_names)) {
  tumor_comb <- AddModuleScore(tumor_comb, ms[i], name = names(ms[i]))
  p <- FeaturePlot(tumor_comb, features = str_c(names(ms[i]), "1"), pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  print(p)
}

#
DimPlot(tumor_comb, reduction = "umap", group.by = "path") +
  coord_fixed() +
  ggtitle("Tumor by pathology type") 

DimPlot(tumor_comb, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("Primary or Secondary Tumor") 

###Get tumor_comb from Ryan, run analysis on individual samples
(table(Idents(tumor_comb)))

#Subset to each sample

s <- c("BC5", "BC6", "BC11", "BC16",
       "BC20", "BC21", "BC22")
raw <- c()
for(i in 1:9) {
  x <- subset(tumor_comb, subset = src == s[i])
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
names(raw) <- s
################################Set optimum number of clusters using Silhoutte scoring

data <-raw[[7]] 
data.nC <- nRes(data, 
                 res = seq(from = 0.1, to = 0.9, by = 0.1))
# Optimum resolution for all samples is 0.1, except for BC21 (res = 0.7)
for (i in c(1:5,7)){
  # raw[[i]] <- raw[[i]] %>% FindClusters(resolution = 0.1)
  print(DimPlot(raw[[i]], pt.size = 1, label = T) +
          coord_fixed() +
          theme(legend.position = "none") +
          ggtitle(str_c(s[i], " basic clustering")))
}

for(i in 6){
  raw[[i]] <- raw[[i]] %>% FindClusters(resolution = 0.7)
  print(DimPlot(raw[[i]], pt.size = 1, label = T) +
          coord_fixed() +
          theme(legend.position = "none") +
          ggtitle(str_c(s[i], " basic clustering")))
}

#Apply Glycolysis module score to each of the primary tumor cells
gl <- list()
Glycolysis <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/Glycolysis.txt", header = FALSE, sep = ",")
Glycolysis <- Glycolysis[,1]

# OxPhos <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/OxPhos.txt", header = FALSE, sep = ",")
# OxPhos <- OxPhos[,1]

Hypoxia <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/Hypoxia.txt", header = FALSE, sep = ",")
Hypoxia <- Hypoxia[,1]

# SASP <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/SASP.txt", header = FALSE, sep = ",")
# SASP <- SASP[,1]

EMT <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/EMT.txt", header = FALSE, sep = ",")
EMT <- EMT[,1]

TNFA <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/TNFA.txt", header = FALSE, sep = ",")
TNFA <- TNFA[,1]

# growth <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/growth.txt", header = FALSE, sep = ",")
# growth <- growth[,1]

gl$Glycolysis <- Glycolysis
gl$Hypoxia <- Hypoxia
gl$EMT <- EMT
gl$TNFA <- TNFA
# gl$growth <- growth

mod_names <- names(gl)

#Glycolysis - Added min.cutoff to isolate for cells most module like
for(i in 1:length(raw)) {
  raw[[i]] <- AddModuleScore(raw[[i]], gl[1], name = names(gl[1]))
  p <- FeaturePlot(raw[[i]], features = str_c(names(gl[1]), "1"), min.cutoff = 0.1, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  pdf(paste0(names(gl[1]), " FeaturePlot ", names(raw[i]), ".pdf"))
  print(p)
  dev.off()
}

#Hypoxia - Added min.cutoff to isolate for cells most module like
for(i in 1:length(raw)) {
  raw[[i]] <- AddModuleScore(raw[[i]], gl[2], name = names(gl[2]))
  p <- FeaturePlot(raw[[i]], features = str_c(names(gl[2]), "1"), min.cutoff = 0.1, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  pdf(paste0(names(gl[2]), " FeaturePlot ", names(raw[i]), ".pdf"))
  print(p)
  dev.off()
}

#EMT - Added min.cutoff to isolate for cells most module like
for(i in 1:length(raw)) {
  raw[[i]] <- AddModuleScore(raw[[i]], gl[3], name = names(gl[3]))
  p <- FeaturePlot(raw[[i]], features = str_c(names(gl[3]), "1"), min.cutoff = 0.1, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  pdf(paste0(names(gl[3]), " FeaturePlot ", names(raw[i]), ".pdf"))
  print(p)
  dev.off()
}

#TNFA - Added min.cutoff to isolate for cells most module like
for(i in 1:length(raw)) {
  raw[[i]] <- AddModuleScore(raw[[i]], gl[4], name = names(gl[4]))
  p <- FeaturePlot(raw[[i]], features = str_c(names(gl[4]), "1"), min.cutoff = 0.1, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  pdf(paste0(names(gl[4]), " FeaturePlot ", names(raw[i]), ".pdf"))
  print(p)
  dev.off()
}
# 
# #Growth - Added min.cutoff to isolate for cells most module like
# for(i in 1:length(raw)) {
#   raw[[i]] <- AddModuleScore(raw[[i]], gl[5], name = names(gl[5]))
#   p <- FeaturePlot(raw[[i]], features = str_c(names(gl[5]), "1"), min.cutoff = 0.1, pt.size = 1,
#                    order = T, cols = c("lightgoldenrod", "darkred"))
#   pdf(paste0(names(gl[5]), " FeaturePlot ", names(raw[i]), ".pdf"))
#   print(p)
#   dev.off()
# }

# ##Perform cell cycle regression on each of the samples
# for(i in 1:length(raw)) {
#   raw[[i]] <- killCC(raw[[i]], cc.regress = "Y", use.pcs = 20)
# }

##Glycolysis and Hypoxia seem to have very high correlation. 
#Output numbe rof cells in each primary tumor sample

# sum(table(Idents(raw[[1]])))
# 8872
# sum(table(Idents(raw[[2]])))
# 825
# sum(table(Idents(raw[[3]])))
# 6357
# sum(table(Idents(raw[[4]])))
# 3115
# sum(table(Idents(raw[[5]])))
# 6886
# sum(table(Idents(raw[[6]])))
# 986
# sum(table(Idents(raw[[7]])))
# 2812
