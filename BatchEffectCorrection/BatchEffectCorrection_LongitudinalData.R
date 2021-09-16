#Harmony on different samples
source("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/scSeurat.R")
set.seed(108)
library("harmony")

load("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure2/OS.combined.postCCR.RData")

DimPlot(OS, reduction = "umap", group.by = "Phase", pt.size = 1, label = F) + 
  coord_fixed() + 
  ggtitle("OS by Source") + 
  scale_color_npg(alpha = 1)

DimPlot(OS, reduction = "umap", pt.size = 1, label = F, group.by = "cond") +
  coord_fixed() + 
  ggtitle("OS Clusters") +
scale_color_npg(alpha = 0.9) 

#Run Harmony
OS <- OS %>%
  RunPCA(pc.genes = OS@var.genes, npcs = 20) %>%
  RunHarmony(group.by.vars = "src", plot_convergence = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.2) ##Note change in resolution

#re-order levels
#Feature Plot for INPP4B, FOSL1, EGR1, FOS
OS$cond <- as.factor(OS$cond)
OS$cond <- factor(as.factor(OS$cond), levels = c("OB", "OS17", "t143b", "NCHOS2", "NCHOS7"))
OS$cond <- revalue(OS$cond, c("OB"="Osteoblast", 
                              "OS17"="OS17", 
                              "t143b" = "143B", 
                              "NCHOS2" = "NCH-OS-2", 
                              "NCHOS7" = "NCH-OS-7"))

#Use FastMNN for batch effect correction; cells still cluster based on cell cycle post-correction
#Note: Running FastMNN after running Harmony gives the same result as without

OS <- RunFastMNN(object.list = SplitObject(OS, split.by = "cond"))
OS <- RunUMAP(OS, reduction = "mnn", dims = 1:20)
OS <- FindNeighbors(OS, reduction = "mnn", dims = 1:20)
OS <- FindClusters(OS, resolution = 0.2)
DimPlot(OS, group.by = c("cond", "src"), ncol = 3)

DimPlot(OS, reduction = "umap", group.by = "Phase", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS by Source") + 
  scale_color_npg(alpha = 1)

