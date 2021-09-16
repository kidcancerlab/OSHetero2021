library("Seurat")
library(ggplot2)
library(future)
library(ggplot2)
source("/gpfs0/home/gdrobertslab/rssxr002/Analysis/scSeurat.R")

# Create Seurat objects and perform initial QC.  Label original source.
#Osteoblast
OB <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0031/outs/filtered_feature_bc_matrix/", spec = "human")
OB <- subset(OB, subset = nCount_RNA <20000 & percent.mt <10)
OB <- OB %>% subset(cells = sample(Cells(OB), 1500))
OB$src <- "OB"
OB$cond <- "OB"


#OS17
os17.cx.raw <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0016/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
os17.cx.raw <- subset(os17.cx.raw, subset = nFeature_RNA >3500 & nCount_RNA <50000 & percent.mt <15)
os17.cx.raw <- os17.cx.raw %>% subset(cells = sample(Cells(os17.cx.raw), 1500))
os17.cx.raw$src <- "OS17_Culture"
os17.cx.raw$cond <- "OS17"

os17.tib.raw <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0018xS0028/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
os17.tib.raw <- subset(os17.tib.raw, subset = nFeature_RNA >3000 & nCount_RNA <60000 & percent.mt <18)
os17.tib.raw <- os17.tib.raw %>% subset(cells = sample(Cells(os17.tib.raw), 1500))
os17.tib.raw$src <- "OS17_Tibia"
os17.tib.raw$cond <- "OS17"

os17.lung.raw <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0018xS0028/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
os17.lung.raw <- subset(os17.lung.raw, subset = nFeature_RNA >1250 & nCount_RNA <60000 & percent.mt <25)
os17.lung.raw <- os17.lung.raw %>% subset(cells = sample(Cells(os17.lung.raw), 1500))
os17.lung.raw$src <- "OS17_Lung"
os17.lung.raw$cond <- "OS17"

#t143b
t143b.cx.raw <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0017/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
t143b.cx.raw <- subset(t143b.cx.raw, subset = nFeature_RNA >4000 & nCount_RNA <74000 & percent.mt <17)
VlnPlot(t143b.cx.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t143b.cx.raw <- t143b.cx.raw %>% subset(cells = sample(Cells(t143b.cx.raw), 1500))
t143b.cx.raw$src <- "t143b_Culture"
t143b.cx.raw$cond <- "t143b"

t143b.tib.raw <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0052/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
t143b.tib.raw <- subset(t143b.tib.raw, subset = nFeature_RNA >1000 & nCount_RNA <30000 & percent.mt <22)
VlnPlot(t143b.tib.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t143b.tib.raw <- t143b.tib.raw %>% subset(cells = sample(Cells(t143b.tib.raw), 1500))
t143b.tib.raw$src <- "t143b_Tibia"
t143b.tib.raw$cond <- "t143b"

t143b.lung.raw <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0019-143B-lung/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
t143b.lung.raw <- subset(t143b.lung.raw, subset = nFeature_RNA >1000 & nCount_RNA <30000 & percent.mt <22)
VlnPlot(t143b.lung.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t143b.lung.raw <- t143b.lung.raw %>% subset(cells = sample(Cells(t143b.lung.raw), 1500))
t143b.lung.raw$src <- "t143b_Lung"
t143b.lung.raw$cond <- "t143b"

#NCHOS2
OS2.cx.raw <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0076/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
OS2.cx.raw <- subset(OS2.cx.raw, subset = nCount_RNA <36000 & percent.mt <50)
VlnPlot(OS2.cx.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS2.cx.raw <- OS2.cx.raw %>% subset(cells = sample(Cells(OS2.cx.raw), 1500))
OS2.cx.raw$src <- "NCHOS2_Flank"
OS2.cx.raw$cond <- "NCHOS2"

OS2.tib.raw<- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0042/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
OS2.tib.raw <- subset(OS2.tib.raw, subset = nFeature_RNA >2500 & nCount_RNA <40000 & percent.mt <20)
VlnPlot(OS2.tib.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS2.tib.raw <- OS2.tib.raw %>% subset(cells = sample(Cells(OS2.tib.raw), 1500))
OS2.tib.raw$src <- "NCHOS2_Tibia"
OS2.tib.raw$cond <- "NCHOS2"

OS2.lung.raw<- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0041/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
OS2.lung.raw <- subset(OS2.lung.raw, subset = nFeature_RNA >2500 & nCount_RNA <60000 & percent.mt <30)
VlnPlot(OS2.lung.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS2.lung.raw <- OS2.lung.raw %>% subset(cells = sample(Cells(OS2.lung.raw), 1500))
OS2.lung.raw$src <- "NCHOS2_Lung"
OS2.lung.raw$cond <- "NCHOS2"

#NCHOS7
OS7.cx.raw <- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts//S0043/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
OS7.cx.raw <- subset(OS7.cx.raw, subset = nCount_RNA <50000 & percent.mt <25)
VlnPlot(OS7.cx.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS7.cx.raw <- OS7.cx.raw %>% subset(cells = sample(Cells(OS7.cx.raw), 1500))
OS7.cx.raw$src <- "NCHOS7_Flank"
OS7.cx.raw$cond <- "NCHOS7"

OS7.tib.raw<- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0034/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
OS7.tib.raw <- subset(OS7.tib.raw, subset = nFeature_RNA >2000 & nCount_RNA <35000 & percent.mt <14)
VlnPlot(OS7.tib.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS7.tib.raw <- OS7.tib.raw %>% subset(cells = sample(Cells(OS7.tib.raw), 1500))
OS7.tib.raw$src <- "NCHOS7_Tibia"
OS7.tib.raw$cond <- "NCHOS7"

OS7.lung.raw<- tenXLoadQC(path10x = "/gpfs0/home/gdrobertslab/lab/Counts/S0055/outs/filtered_feature_bc_matrix/", spec = "mixHuman")
OS7.lung.raw <- subset(OS7.lung.raw, subset = nCount_RNA <42000 & percent.mt <20)
VlnPlot(OS7.lung.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
OS7.lung.raw <- OS7.lung.raw %>% subset(cells = sample(Cells(OS7.lung.raw), 1500))
OS7.lung.raw$src <- "NCHOS7_Lung"
OS7.lung.raw$cond <- "NCHOS7"

# Merge into a single Seurat object
OS <- merge(OB, y = c(os17.cx.raw, os17.tib.raw, os17.lung.raw, 
                               t143b.cx.raw, t143b.tib.raw, t143b.lung.raw, 
                               OS2.cx.raw, OS2.tib.raw, OS2.lung.raw, 
                               OS7.cx.raw, OS7.tib.raw, OS7.lung.raw),
            add.cell.ids = c("OB", "OS17_Culture", "OS17_Tibia", "OS17_Lung", 
                             "t143b_Culture", "t143b_Tibia", "t143b_Lung",
                             "NCHOS2_Flank", "NCHOS2_Tibia", "NCHOS2_Lung", 
                             "NCHOS7_Flank", "NCHOS7_Tibia", "NCHOS7_Lung"),
            project = "Heterogeneity")

# Process and cluster
OS <- NormalizeData(OS) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(pc.genes = os.17@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

save(OS, file = "OS.combined.OBpreCCR.RData")
# CCR
# Attempt to regress out the effects of cell cycle on these tumor cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
OS <- CellCycleScoring(object = OS, s.features = s.genes,
                         g2m.features = g2m.genes, set.ident = TRUE)

OS <- ScaleData(object = OS, vars.to.regress = c("S.Score", "G2M.Score"),
                  features = rownames(x = OS)) 
OS <- RunPCA(OS, pc.genes = OS@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

save(OS, file = "OS.combined.OBpostCCR.RData")

OS$cond <- as.factor(OS$cond)
OS$cond <- factor(as.factor(OS$cond), levels = c("OS17", "t143b", "NCHOS2", "NCHOS7", "OB"))

# Plot the data (what do you think of these tweaks in terms of color, transparency, and labels?)
DimPlot(OS, reduction = "umap", group.by = "cond", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS by Source") + 
  scale_color_npg(alpha = 1)
DimPlot(OS, reduction = "umap", group.by = "src", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS by Source")
DimPlot(OS, reduction = "umap", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS Clusters") + 
  scale_color_npg(alpha = 0.7)


