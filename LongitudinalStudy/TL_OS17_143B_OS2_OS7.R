source("R:/RESRoberts/Bioinformatics/Analysis/scSeurat.R")
source("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/Downstream.v2.R")
#loading libraries
library(Seurat)
library(future)
library(ggplot2)
# library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(SingleR)
library(dplyr) # %<%
library(pheatmap)
library(RColorBrewer)
library(viridis) # inferno color palette
library(grid)

#Merge by condition
#OS17
os17.tib.raw <- tenXLoadQC("R:/RESRoberts/Bioinformatics/scRNAOuts/S0018xS0028/filtered_feature_bc_matrix/", spec = "mixHuman")
os17.tib.raw <- subset(os17.tib.raw, subset = nFeature_RNA >3000 & nCount_RNA <60000 & percent.mt <18)
os17.tib.raw$src <- "OS17_Tibia"
os17.tib.raw$cond <- "OS17"

os17.lung.raw <- tenXLoadQC("R:/RESRoberts/Bioinformatics/scRNAOuts/S0024xS0029/filtered_feature_bc_matrix/", spec = "mixHuman")
os17.lung.raw <- subset(os17.lung.raw, subset = nFeature_RNA >1250 & nCount_RNA <60000 & percent.mt <25)
os17.lung.raw$src <- "OS17_Lung"
os17.lung.raw$cond <- "OS17"

  #Subset number of cells in each Seurat object to the least in the group
  a <- table(os17.tib.raw$orig.ident)
  b <- table(os17.lung.raw$orig.ident)
  
  list <- c(a, b)
  min <- min(list)
  
  os17.tib.raw <- subset(os17.tib.raw, cells = sample(Cells(os17.tib.raw), min))
  os17.lung.raw <- subset(os17.lung.raw, cells = sample(Cells(os17.lung.raw), min))

#t143b
t143b.cx.raw <- tenXLoadQC(path10x = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0017-143b/filtered_feature_bc_matrix/", spec = "mixHuman")
t143b.cx.raw <- subset(t143b.cx.raw, subset = nFeature_RNA >4000 & nCount_RNA <74000 & percent.mt <17)
VlnPlot(t143b.cx.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t143b.cx.raw$src <- "t143b_Cx"
t143b.cx.raw$cond <- "t143b"

t143b.tib.raw <- tenXLoadQC(path10x = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0052-143B-Tibia/filtered_feature_bc_matrix/", spec = "mixHuman")
t143b.tib.raw <- subset(t143b.tib.raw, subset = nFeature_RNA >1000 & nCount_RNA <30000 & percent.mt <22)
VlnPlot(t143b.tib.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t143b.tib.raw$src <- "t143b_Tibia"
t143b.tib.raw$cond <- "t143b"

t143b.lung.raw <- tenXLoadQC(path10x = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0019-143B-lung/filtered_feature_bc_matrix/", spec = "mixHuman")
t143b.lung.raw <- subset(t143b.lung.raw, subset = nFeature_RNA >1000 & nCount_RNA <30000 & percent.mt <22)
VlnPlot(t143b.lung.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t143b.lung.raw$src <- "t143b_Lung"
t143b.lung.raw$cond <- "t143b"

  #Subset number of cells in each Seurat object to the least in the group
  a <- table(t143b.tib.raw$orig.ident)
  b <- table(t143b.lung.raw$orig.ident)
  
  list <- c(a, b)
  min <- min(list)
  
  t143b.tib.raw <- subset(t143b.tib.raw, cells = sample(Cells(t143b.tib.raw), min))
  t143b.lung.raw <- subset(t143b.lung.raw, cells = sample(Cells(t143b.lung.raw), min))

#NCHOS2
OS2.tib.raw<- tenXLoadQC(path10x = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0042-NCHOS2-tibia/filtered_feature_bc_matrix/", spec = "mixHuman")
OS2.tib.raw <- subset(OS2.tib.raw, subset = nFeature_RNA >2500 & nCount_RNA <40000 & percent.mt <20)
OS2.tib.raw$src <- "NCHOS2_Tibia"
VlnPlot(OS2.tib.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS2.tib.raw$cond <- "NCHOS2"

OS2.lung.raw<- tenXLoadQC(path10x = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0041-NCHOS2-lung/filtered_feature_bc_matrix/", spec = "mixHuman")
OS2.lung.raw <- subset(OS2.lung.raw, subset = nFeature_RNA >2500 & nCount_RNA <60000 & percent.mt <30)
OS2.lung.raw$src <- "NCHOS2_Lung"
VlnPlot(OS2.lung.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS2.lung.raw$cond <- "NCHOS2"

  #Subset number of cells in each Seurat object to the least in the group
  a <- table(OS2.tib.raw$orig.ident)
  b <- table(OS2.lung.raw$orig.ident)
  
  list <- c(a, b)
  min <- min(list)
  
  OS2.tib.raw <- subset(OS2.tib.raw, cells = sample(Cells(OS2.tib.raw), min))
  OS2.lung.raw <- subset(OS2.lung.raw, cells = sample(Cells(OS2.lung.raw), min))

#NCHOS7
OS7.tib.raw<- tenXLoadQC(path10x = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0034-NCHOS7-tib/filtered_feature_bc_matrix/", spec = "mixHuman")
OS7.tib.raw <- subset(OS7.tib.raw, subset = nFeature_RNA >2000 & nCount_RNA <35000 & percent.mt <14)
OS7.tib.raw$src <- "NCHOS7_Tibia"
VlnPlot(OS7.tib.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS7.tib.raw$cond <- "NCHOS7"

OS7.lung.raw<- tenXLoadQC(path10x = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0055-NCHOS7-lung/filtered_feature_bc_matrix/", spec = "mixHuman")
OS7.lung.raw <- subset(OS7.lung.raw, subset = nCount_RNA <42000 & percent.mt <20)
OS7.lung.raw$src <- "NCHOS7_Lung"
VlnPlot(OS7.lung.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS7.lung.raw$cond <- "NCHOS7"

  #Subset number of cells in each Seurat object to the least in the group
  a <- table(OS7.tib.raw$orig.ident)
  b <- table(OS7.lung.raw$orig.ident)
  
  list <- c(a, b)
  min <- min(list)
  
  OS7.tib.raw <- subset(OS7.tib.raw, cells = sample(Cells(OS7.tib.raw), min))
  OS7.lung.raw <- subset(OS7.lung.raw, cells = sample(Cells(OS7.lung.raw), min))


# Merge into a single Seurat object
OS17.TL <- merge(os17.tib.raw, y = c(os17.lung.raw),
                   add.cell.ids = c("OS17_Tibia", "OS17_Lung"),
                   project = "Heterogeneity")

t143b <- merge(t143b.tib.raw, y = c(t143b.lung.raw),
               add.cell.ids = c("t143b_Tibia", "t143b_Lung"),
               project = "Heterogeneity")
  
OS2 <- merge(OS2.tib.raw, y = c(OS2.lung.raw),
            add.cell.ids = c("OS2_Tibia", "OS2_Lung"),
            project = "Heterogeneity")

OS7 <- merge(OS7.tib.raw, y = c(OS7.lung.raw),
             add.cell.ids = c("OS7_Tibia", "OS7_Lung"),
             project = "Heterogeneity")

# Process and cluster
OS17.TL <- NormalizeData(OS17.TL) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(pc.genes = OS17.TL@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

t143b <- NormalizeData(t143b) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(pc.genes = t143b@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

OS2 <- NormalizeData(OS2) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(pc.genes = OS2@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

OS7 <- NormalizeData(OS7) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(pc.genes = OS7@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

# CCR
# Attempt to regress out the effects of cell cycle on these tumor cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#OS17
OS17.TL <- CellCycleScoring(object = OS17.TL, s.features = s.genes,
                            g2m.features = g2m.genes, set.ident = TRUE)

OS17.TL <- ScaleData(object = OS17.TL, vars.to.regress = c("S.Score", "G2M.Score"),
                     features = rownames(x = OS17.TL)) 
OS17.TL <- RunPCA(OS17.TL, pc.genes = OS17.TL@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

#t143b
t143b <- CellCycleScoring(object = t143b, s.features = s.genes,
                            g2m.features = g2m.genes, set.ident = TRUE)

t143b <- ScaleData(object = t143b, vars.to.regress = c("S.Score", "G2M.Score"),
                     features = rownames(x = t143b)) 
t143b <- RunPCA(t143b, pc.genes = t143b@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

#OS2
OS2 <- CellCycleScoring(object = OS2, s.features = s.genes,
                       g2m.features = g2m.genes, set.ident = TRUE)

OS2 <- ScaleData(object = OS2, vars.to.regress = c("S.Score", "G2M.Score"),
                features = rownames(x = OS2)) 
OS2 <- RunPCA(OS2, pc.genes = OS2@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

#OS7
OS7 <- CellCycleScoring(object = OS7, s.features = s.genes,
                       g2m.features = g2m.genes, set.ident = TRUE)

OS7 <- ScaleData(object = OS7, vars.to.regress = c("S.Score", "G2M.Score"),
                features = rownames(x = OS7)) 
OS7 <- RunPCA(OS7, pc.genes = OS7@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

# Plot the data 
DimPlot(OS17.TL, reduction = "umap", group.by = "src", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS17.TL by Source") + 
  scale_color_npg(alpha = 0.7)
DimPlot(OS17.TL, reduction = "umap", pt.size = 1, label = T, split.by = "src") + 
  coord_fixed() + 
  ggtitle("OS17 Clusters") + 
  scale_color_npg(alpha = 0.7)

DimPlot(t143b, reduction = "umap", group.by = "src", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("t143b by Source") + 
  scale_color_npg(alpha = 0.7)
DimPlot(t143b, reduction = "umap", pt.size = 1, label = T, split.by = "src") + 
  coord_fixed() + 
  ggtitle("t143b Clusters") + 
  scale_color_npg(alpha = 1)

DimPlot(OS2, reduction = "umap", group.by = "src", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS by Source") + 
  scale_color_npg(alpha = 0.7)
DimPlot(OS2, reduction = "umap", pt.size = 1, label = T, split.by = "src") + 
  coord_fixed() + 
  ggtitle("OS Clusters") + 
  scale_color_npg(alpha = 0.7)

DimPlot(OS7, reduction = "umap", group.by = "src", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS by Source") + 
  scale_color_npg(alpha = 0.7)
DimPlot(OS7, reduction = "umap", pt.size = 1, label = T, split.by = "src") + 
  coord_fixed() + 
  ggtitle("OS Clusters") + 
  scale_color_npg(alpha = 0.7)

save.image(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/OS_All.RData")

#Plot the cluster distribution in each sample
temp <- colSums(table(Idents(OS2), OS2$src))
min <- min(temp)
cell.ids <- sample(colnames(OS2["OS2"]))

table(Idents(OS17.TL), OS17.TL$src)
prop <- prop.table(table(Idents(OS17.TL), OS17.TL$src), margin = 2)*100
colSums(prop)

#     Tibia       Lung
# 0   50.263250   46.577747
# 1   43.453843   10.179010
# 2   1.053001   39.698140
# 3   5.229905   3.545104

table(Idents(t143b), t143b$src)
prop <- prop.table(table(Idents(t143b), t143b$src), margin = 2)*100
colSums(prop)

#     Tibia        Lung
# 0   97.8963771    2.2204908
# 1   1.1492014   92.3061940
# 2   0.4479938   3.3502143
# 3   0.5064277   2.1231009

table(Idents(OS2), OS2$src)
prop <- prop.table(table(Idents(OS2), OS2$src), margin = 2)*100
colSums(prop)

#     Tibia       Lung
# 0   82.281553   20.327670
# 1   5.339806    71.480583
# 2   12.378641   8.191748

table(Idents(OS7), OS7$src)
prop <- prop.table(table(Idents(OS7), OS7$src), margin = 2)*100
colSums(prop)

#     Tibia       Lung
# 0    0.000000   71.013490
# 1   59.425804    1.383604
# 2   40.574196    1.176064
# 3    0.000000   26.426842

save(list = c(cx.sub, lung.sub, tib.sub, os17), 
           file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/OS17CTL.sub.RData")
save(list = c(os17.tib.raw, os17.lung.raw, OS17.TL), 
     file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Fig3_OS17TL.RData")
save(list = c(t143b.cx.raw, t143b.lung.raw, t143b.tib.raw, t143b), 
     file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Fig3_t143bTL.RData")
save(list = c(OS2.tib.raw, OS2.lung.raw, OS2), 
     file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Fig3_OS2TL.RData")
save(list = c(OS7.tib.raw, OS7.lung.raw, OS7), 
     file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Fig3_OS7TL.RData")

#remove these objects and then run CCR on OS and save that file separately; create OS.listonce cleared and re-started
save(list = c(OS, OS.CCR), 
     file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Fig3_OS7TL.RData")

#Modify data to account for optimal number of clusters using nRes and pSil functions
OS17.TL.nRes <- nRes(OS17.TL, 
                     res = seq(from = 0.1, to = 0.2, by = 0.05))

plot <- pSil(OS17.TL.nRes , 0.15)

OS17.TL <- FindClusters(OS17.TL, resolution = 0.20)

DimPlot(OS17.TL, reduction = "umap", pt.size = 1, label = T, split.by = "src") +
  coord_fixed() +
  NoLegend() + NoAxes()

#Modify data to account for optimal number of clusters using nRes and pSil functions
t143b.nRes <- nRes(t143b, 
                     res = seq(from = 0.1, to = 0.4, by = 0.05))

plot <- pSil(t143b.nRes , 0.15)

t143b <- FindClusters(t143b, resolution = 0.15)

DimPlot(t143b, reduction = "umap", pt.size = 1, label = T, split.by = "src") +
  coord_fixed()

#Modify data to account for optimal number of clusters using nRes and pSil functions
OS2.nRes <- nRes(OS2, 
                     res = seq(from = 0.1, to = 0.4, by = 0.05))

plot <- pSil(OS2.nRes , 0.15)

OS2 <- FindClusters(OS2, resolution = 0.15)

DimPlot(OS2, reduction = "umap", pt.size = 1, label = T, split.by = "src") +
  coord_fixed() +
  NoLegend() + NoAxes()

#Modify data to account for optimal number of clusters using nRes and pSil functions
OS7.nRes <- nRes(OS7, 
                     res = seq(from = 0.1, to = 0.2, by = 0.05))

plot <- pSil(OS7.nRes , 0.1)

OS7 <- FindClusters(OS7, resolution = 0.1)

DimPlot(OS7, reduction = "umap", pt.size = 1, label = T, split.by = "src") +
  coord_fixed() +
  NoLegend() + NoAxes()


#For lab meeting
cx.raw.nC <- nRes(cx.raw, 
                  res = seq(from = 0.1, to = 0.5, by = 0.1))
plot <- pSil(cx.raw.nC, 0.5)

cx.raw <- FindClusters(cx.raw, resolution = 0.2)

plot <- DimPlot(cx.raw, reduction = "umap", pt.size = 1, label = T, split.by = "src") +
  coord_fixed()
# 
# #Discussion from lab meeting - manually assign cluster IDs for OS17-culture
# cx.raw.new <- cx.raw
# library(grDevices)
# library(grDevices, lib.loc = "C:/Program Files/R/R-3.6.2/library")
# windows()
# cx.raw.new <- CellSelector(plot = plot, object = cx.raw.new, ident = "4")
# 
# Idents(cx.raw) <- factor(x = Idents(cx.raw), levels = sort(levels(cx.raw)))
# Clusters <- Idents(object = cx.raw)
# names(Clusters) <- colnames(x = cx.raw)
# cx.raw.new <- AddMetaData(
#   object = cx.raw.new,
#   metadata = Clusters,
#   col.name = 'cluster.idents'
# )
# head(x = cx.raw.new[[]])
# 
# clusters <- cx.raw.new@meta.data[["cluster.idents"]]
# dist.matrix <- dist(x = Embeddings(object = cx.raw.new[["pca"]])[, 1:20])
# sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
# # cx.raw.new$sil <- sil[, 3]
# 
# windows()
# # col= pal_npg("nrc")(n)
# plot(sil,
#      col = pal_npg("nrc")(4)) # with cluster-wise coloring
# 
# pdf("temp.pdf")
# DimPlot(cx.raw.new, reduction = "umap", pt.size = 1, label = T) +
#   coord_fixed()
# dev.off()

#subset after finding clusters
# re-arrange and re-name levels 

# temp <- OS17.TL
# temp$src <- factor(as.factor(temp$src), levels = c("OS17_Tibia", "OS17_Lung"))
library(plyr)
# temp$src <- revalue(temp$src, c("OS17_Tibia"="Tibia", "OS17_Lung"="Lung"))
OS17.TL$src <- as.factor(OS17.TL$src)
OS17.TL$src <- factor(as.factor(OS17.TL$src), levels = c("OS17_Tibia", "OS17_Lung"))
OS17.TL$src <- revalue(OS17.TL$src, c("OS17_Tibia"="Tibia", "OS17_Lung"="Lung"))

OS2$src <- as.factor(OS2$src)
OS2$src <- factor(as.factor(OS2$src), levels = c("NCHOS2_Tibia", "NCHOS2_Lung"))
OS2$src <- revalue(OS2$src, c("NCHOS2_Tibia"="Tibia", "NCHOS2_Lung"="Lung"))

OS7$src <- as.factor(OS7$src)
OS7$src <- factor(as.factor(OS7$src), levels = c("NCHOS7_Tibia", "NCHOS7_Lung"))
OS7$src <- revalue(OS7$src, c("NCHOS7_Tibia"="Tibia", "NCHOS7_Lung"="Lung"))

OS17.TL <- OS17.TL %>% FindClusters(resolution = 0.2)
t143b <- t143b %>% FindClusters(resolution = 0.15)
OS2 <- OS2 %>% FindClusters(resolution = 0.15)
OS7 <- OS7 %>% FindClusters(resolution = 0.2)

pdf("OS17.TL_fig3.pdf", width = 5, height = 5)
OS17.TL %>% FindClusters(resolution = 0.2) %>% 
  subset(cells = sample(Cells(OS17.TL), 3000)) %>% 
  DimPlot(reduction = "umap", pt.size = 1, label = T, split.by = "src") + 
  coord_fixed() + NoLegend() + NoAxes()+
  scale_color_npg(alpha = 1)
dev.off()

pdf("t143b_fig3.pdf", width = 5, height = 5)
t143b %>% FindClusters(resolution = 0.15) %>% 
  subset(cells = sample(Cells(t143b), 3000)) %>% 
  DimPlot(reduction = "umap", pt.size = 1, label = T, split.by = "src") + 
  coord_fixed() + NoLegend() + NoAxes()+
  scale_color_npg(alpha = 1)
dev.off()

pdf("OS2_fig3.pdf", width = 5, height = 5)
OS2 %>% FindClusters(resolution = 0.15) %>% 
  subset(cells = sample(Cells(OS2), 3000)) %>% 
  DimPlot(reduction = "umap", pt.size = 1, label = T, split.by = "src") + 
  coord_fixed() + NoLegend() + NoAxes()+
  scale_color_npg(alpha = 1)
dev.off()

pdf("OS7_fig3.pdf", width = 10, height = 10)
OS7 %>% FindClusters(resolution = 0.2) %>% 
  subset(cells = sample(Cells(OS7), 3000)) %>% 
  DimPlot(reduction = "umap", pt.size = 1, label = T, split.by = "src") +
  ggtitle("OS7")+
  coord_fixed() + NoLegend() + NoAxes()+
  scale_color_npg(alpha = 1)
dev.off()

#create heatmap 
B.list <- list(OS17 = OS17.TL,
               t143b = t143b,
               OS2 = OS2,
               OS7 = OS7)

#Heatmap and em.hm files
em.hm.list <- list()
for (i in 1:length(B.list)) {
  em.hm.list[[i]] <- DGEA(B.list[[i]]) 
}

library(data.table)
for(i in 1:(length(em.hm.list))) {
  em.hm.list[[i]] <- setDT(em.hm.list[[i]], keep.rownames = TRUE)[]
}

temp1 <- em.hm.list[[1]]
temp2 <- em.hm.list[[2]]
temp3 <- em.hm.list[[3]]
temp4 <- em.hm.list[[4]]

cx <- inner_join(temp1, temp2, by = "rn")
pd <- inner_join(temp3, temp4, by = "rn")
cx.pd <- inner_join(cx, pd, by = "rn")

##Write table to edit rownames to easily convert to dataframe
write.table(cx.pd, file=('R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2021/cxpd.tsv'), quote=FALSE, sep='\t')
cx.pd <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2021/cxpd.tsv", header=T)

#remove "HALLMARK_"
c <- rownames(cx.pd)
c <- gsub("HALLMARK_", "", c)
#remove "_" by removing special characters
c <- gsub("_", " ", c)
rownames(cx.pd) <- c

cx.pd.log <- -log10(cx.pd) #log transform
library(pheatmap)   
cx.pd_transpose <- as.data.frame(t(cx.pd.log))
df <- as.matrix(cx.pd_transpose)
#Write table to print supplemental information
write.table(df, file=('R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2021/cxpd_plot.tsv'), quote=FALSE, sep='\t')

subj<-c(paste0("A", 0:3),
        paste0("B", 0:3),
        paste0("C", 0:2),
        paste0("D", 0:3))
rownames(df)<-subj
aka2 = data.frame(ID = factor(c(rep("OS17", 4),
                                rep("t143B", 4),
                                rep("OS2", 3),
                                rep("OS7", 4)
                                )))
rownames(aka2)<-subj
aka3 = list(ID = c(OS17 = "#E64B35FF", t143B = "#4DBBD5FF", OS2 = "#00A087FF", OS7 = "#3C5488FF"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)
#calculate percentage of cells in lung for each cluster
aka4 <- data.frame (Lpercent = c(46.577747,
                                 10.17901,
                                 39.69814,
                                 3.545104,
                                 2.2204908,
                                 92.306194,
                                 3.3502143,
                                 2.1231009,
                                 20.32767,
                                 71.480583,
                                 8.191748,
                                 71.01349,
                                 1.383604,
                                 1.176064,
                                 26.426842))
rownames(aka4)<-subj
aka5 = (c(rep("OS17", 4), rep("t143B", 4), rep("OS2", 3), rep("OS7", 4)))
column_ha = HeatmapAnnotation(model = aka5, bar1 = anno_barplot(aka4))
pdf("temp.pdf", height = 7, width = 8)
Heatmap(t(df), 
        col = col, 
        name = "-log10(p)", 
        cluster_columns = FALSE,
        top_annotation = column_ha)+ geom_text(aes(label = pvalue))
dev.off()

col.breaks=seq(-log10(1),min(max(-log10(cx.pd))+1,18),by=0.5)
col=inferno(length(col.breaks)) # library(viridis)
col=c("white",colorRampPalette(brewer.pal(n = 7, name ="Reds"))(50))
pheatmap(-log10(cx.pd,cluster_rows = TRUE,cluster_cols = TRUE,
         cellwidth = 5,cellheight = 7,treeheight_row = 0,treeheight_col=0,
         color = col,scale='none',breaks=col.breaks,fontsize = 8))


#IPA analysis
# Decided to proceed with A2, B1, (B2, B3), C1, D0 and D3
# Find differentially expressed features selected cluster and remaining cells
#OS17
A0.markers <- FindMarkers(OS17.TL, ident.1 = "0", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
A1.markers <- FindMarkers(OS17.TL, ident.1 = "1", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
A2.markers <- FindMarkers(OS17.TL, ident.1 = "2", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
A3.markers <- FindMarkers(OS17.TL, ident.1 = "3", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)

#remove unnecessary columns
myList <- list(A0.markers, A1.markers, A2.markers, A3.markers)
myList <- lapply(myList, function(x) { x["p_val"] <- NULL; x[, 2:3] <- NULL; x })
A0.markers <- myList[[1]]
A1.markers <- myList[[2]]
A2.markers <- myList[[3]]
A3.markers <- myList[[4]]
  
# view results
# head(A2.markers)
# FeaturePlot(object = OS17.TL, features = 'NFKBIA')
B0.markers <- FindMarkers(t143b, ident.1 = "0", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
B1.markers <- FindMarkers(t143b, ident.1 = "1", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
B2.markers <- FindMarkers(t143b, ident.1 = "2", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
B3.markers <- FindMarkers(t143b, ident.1 = "3", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)

myList <- list(B0.markers, B1.markers, B2.markers, B3.markers)
myList <- lapply(myList, function(x) { x["p_val"] <- NULL; x[, 2:3] <- NULL; x })
B0.markers <- myList[[1]]
B1.markers <- myList[[2]]
B2.markers <- myList[[3]]
B3.markers <- myList[[4]]

C0.markers <- FindMarkers(OS2, ident.1 = "0", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
C1.markers <- FindMarkers(OS2, ident.1 = "1", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
C2.markers <- FindMarkers(OS2, ident.1 = "2", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)

myList <- list(C0.markers, C1.markers, C2.markers)
myList <- lapply(myList, function(x) { x["p_val"] <- NULL; x[, 2:3] <- NULL; x })
C0.markers <- myList[[1]]
C1.markers <- myList[[2]]
C2.markers <- myList[[3]]

D0.markers <- FindMarkers(OS7, ident.1 = "0", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
D1.markers <- FindMarkers(OS7, ident.1 = "1", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
D2.markers <- FindMarkers(OS7, ident.1 = "2", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
D3.markers <- FindMarkers(OS7, ident.1 = "3", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)

myList <- list(D0.markers, D1.markers, D2.markers, D3.markers)
myList <- lapply(myList, function(x) { x["p_val"] <- NULL; x[, 2:3] <- NULL; x })
D0.markers <- myList[[1]]
D1.markers <- myList[[2]]
D2.markers <- myList[[3]]
D3.markers <- myList[[4]]

library(xlsx)
write.xlsx(A0.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/A0.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(A1.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/A1.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(A2.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/A2.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(A3.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/A3.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)


write.xlsx(B0.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/B0.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(B1.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/B1.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(B2.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/B2.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(B3.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/B3.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

write.xlsx(C0.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/C0.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(C1.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/C1.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(C2.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/C2.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

write.xlsx(D0.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/D0.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(D1.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/D1.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(D2.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/D2.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(D3.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/D3.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

#Uploaded individual datasets on IPA; ran Expression analysis with adj p value cutoff at 0.01
df <- read.delim("IPA_updown.txt", header = TRUE)
rownames(df) <- df[,1]
df[,1] <- NULL


# cx.pd <- cx.pd[-(23:30),] #remove rows with names NA - figure out why we have NAs!!
cx.pd.log <- -log10(cx.pd) #log transform
library(pheatmap)   

cx.pd_transpose <- as.data.frame(t(cx.pd.log))
df <- as.matrix(cx.pd_transpose)

subj<-c(paste0("A", 0:3),
        paste0("B", 0:3),
        paste0("C", 0:2),
        paste0("D", 0:3))
rownames(df)<-subj
aka2 = data.frame(ID = factor(c(rep("OS17", 4),
                                rep("t143B", 4),
                                rep("OS2", 3),
                                rep("OS7", 4)
)))
rownames(aka2)<-subj
aka3 = list(ID = c(OS17 = "#E64B35FF", t143B = "#4DBBD5FF", OS2 = "#00A087FF", OS7 = "#3C5488FF"))
# 
# df[] <- lapply(df, gsub, pattern='HALLMARK_', replacement='')
# df

cols <- makeColorRampPalette(c("snow1", "red"), # distances 3 to max(distmat) colored from green to black
                             100)

pheatmap(t(scale(df)),
         color = cols, 
         annotation_col = aka2, 
         annotation_colors = aka3,
         annotation_legend = TRUE,
         gaps_col =  4,
         show_colnames = T, show_rownames = T, cluster_rows = T, 
         cluster_cols = T, legend = TRUE, 
         clustering_distance_rows = "euclidean", border_color = FALSE)

Heatmap(t(scale(df)), 
        col = col, 
        name = "-log10(p)", 
        top_annotation = column_ha)

color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)

##Apply glycolysis, Hypoxia, tnfa, emt modules and compare between clusters
data <- OS17.TL

for(i in 1:length(mod_names)) {
  data <- AddModuleScore(data, gl[i], name = names(gl[i]))
}
RidgePlot(data, features = str_c(names(gl[6]), "1"))
# VlnPlot(data, features = str_c(names(gl[3]), "1"))
DotPlot(data, features = TNFA, cluster.idents = TRUE)

#noted a clear differnece in the TNFA_VIA_NFkB module
#Identify genes that are most differntially expressed in this module

data <- OS17.TL
pdf("Dotplot_OS17.pdf", height = 30, width = 7)
DotPlot(data, features = TNFA, cluster.idents = TRUE)+coord_flip()
dev.off()

data <- t143b
pdf("Dotplot_143B.pdf", height = 30, width = 7)
DotPlot(data, features = TNFA, cluster.idents = TRUE)+coord_flip()
dev.off()

data <- OS2
pdf("Dotplot_OS2.pdf", height = 30, width = 7)
DotPlot(data, features = TNFA, cluster.idents = TRUE)+coord_flip()
dev.off()

data <- OS7
pdf("Dotplot_OS7.pdf", height = 30, width = 7)
DotPlot(data, features = TNFA, cluster.idents = TRUE)+coord_flip()
dev.off()

################################Pathway enrichment analysis 

P1 <- df %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  full_join(aka4 %>% rownames_to_column(var = "Sample")) %>%
  pivot_longer(c(-Sample, -Lpercent),
               names_to = "Pathway",
               values_to = "pval",
               names_repair = "minimal") %>%
  mutate(Signif = pval >= (-1 * log10(0.05)),
         sample_order = str_remove(Sample, "[0-9]") %>%
           rank() * 1000 - Lpercent,
         Sample = reorder(Sample, sample_order)) %>%
  ggplot(., aes(x = Sample, y = Pathway, fill = Signif)) + 
  geom_tile() +
  geom_text(aes(label = sprintf("%0.2f", pval))) +
  scale_fill_manual(values = c("gray", "#CC3333")) +
  scale_y_discrete(limits = rev) +
  theme(legend.position = "none") +
  ylab("") 
P2 <- aka4 %>%
  rownames_to_column(var = "Sample") %>%
  full_join(aka2 %>%
              rownames_to_column(var = "Sample")) %>%
  mutate(sample_order = str_remove(Sample, "[0-9]") %>%
           rank() * 1000 - Lpercent,
         Sample = reorder(Sample, sample_order),
         ID = reorder(ID, sample_order)) %>%
  ggplot(., aes(x = Sample, y = Lpercent, fill = ID)) +
  geom_bar(stat = "identity") +
  ylab("Percent\nin lung") +
  xlab("") + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "top", ) +
  labs(fill = "") 

pdf("Pvalueplot.pdf", height = 12, width = 12)
cowplot::plot_grid(P2, P1, ncol = 1, align = "v",rel_heights = c(2, 10))
dev.off()
