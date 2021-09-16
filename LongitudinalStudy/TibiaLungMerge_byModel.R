source("R:/RESRoberts/Bioinformatics/Analysis/scSeurat.R")

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


# Merge (by model) tibi and lung datasets into a single Seurat object
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

save(list = c(os17.tib.raw, os17.lung.raw, OS17.TL), 
     file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Fig3_OS17TL.RData")
save(list = c(t143b.cx.raw, t143b.lung.raw, t143b.tib.raw, t143b), 
     file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Fig3_t143bTL.RData")
save(list = c(OS2.tib.raw, OS2.lung.raw, OS2), 
     file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Fig3_OS2TL.RData")
save(list = c(OS7.tib.raw, OS7.lung.raw, OS7), 
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

#optimal clsuter resolution set usign Silhouette scoring
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
