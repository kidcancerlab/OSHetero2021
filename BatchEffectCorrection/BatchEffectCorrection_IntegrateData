#BatchEffectCorrectionSCTransform - all
source("/gpfs0/home/gdrobertslab/rssxr002/Analysis/scSeurat.R")
library(dplyr)
library(ggplot2)
library(harmony)
library(cowplot)
library(Seurat)
set.seed(100)

load("/gpfs0/home/gdrobertslab/rssxr002/Analysis/OS.combined.postCCR.RData")
#Integrate using Seurats SCTransform Integration workflow
OSrep.list <- SplitObject(OS, split.by = "src")

for (i in 1:length(OSrep.list)) {
  OSrep.list[[i]] <- SCTransform(OSrep.list[[i]], verbose = FALSE)
}

OSrep.features <- SelectIntegrationFeatures(object.list = OSrep.list, nfeatures = 3000)
OSrep.list <- PrepSCTIntegration(object.list = OSrep.list, anchor.features = OSrep.features, 
                                 verbose = FALSE)
OSrep.anchors <- FindIntegrationAnchors(object.list = OSrep.list, normalization.method = "SCT", 
                                        anchor.features = OSrep.features, verbose = FALSE)
OSrep.integrated <- IntegrateData(anchorset = OSrep.anchors, normalization.method = "SCT", 
                                  verbose = FALSE)

OSrep.integrated <- RunPCA(OSrep.integrated, verbose = FALSE)
OSrep.integrated <- RunUMAP(OSrep.integrated, dims = 1:30)

saveRDS(OSrep.integrated, file = "OSrep.integrated.rds")

pdf("OSrep.integrated.Phase.pdf")
DimPlot(OS, reduction = "umap", group.by = "Phase", pt.size = 1, label = F) + 
  coord_fixed() + 
  ggtitle("OS by Source") + 
  scale_color_npg(alpha = 1)
dev.off()

pdf("OSrep.integrated.pdf")
DimPlot(OS, reduction = "umap", pt.size = 1, label = T) +
  coord_fixed() + 
  ggtitle("OS Clusters") 
scale_color_npg(alpha = 1) 
dev.off()
