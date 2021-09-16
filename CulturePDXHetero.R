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

#OS17
# Create Seurat objects and perform initial QC. Label original source.
cx.raw <- tenXLoadQC("R:/RESRoberts/Bioinformatics/scRNAOuts/S0016-zymo/filtered_feature_bc_matrix/", spec = "mixHuman")
cx.raw <- subset(cx.raw, subset = nFeature_RNA >3500 & nCount_RNA <50000 & percent.mt <15)
cx.raw$src <- "Culture"

cx.raw <- NormalizeData(cx.raw) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(pc.genes = cx.raw@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

DimPlot(cx.raw, reduction = "umap", pt.size = 1, label = T) + 
  coord_fixed() + 
  scale_color_npg(alpha = 0.7)

# Attempt to regress out the effects of cell cycle on these tumor cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cx.raw <- CellCycleScoring(object = cx.raw, s.features = s.genes,
                           g2m.features = g2m.genes, set.ident = TRUE)

cx.raw <- ScaleData(object = cx.raw, vars.to.regress = c("S.Score", "G2M.Score"),
                    features = rownames(x = cx.raw)) 
cx.raw <- RunPCA(cx.raw, pc.genes = cx.raw@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

#Find optimal clustering resolution
cx.raw.nC <- nRes(cx.raw, 
                  res = seq(from = 0.1, to = 0.2, by = 0.05))
plot <- pSil(cx.raw.nC, 0.15)

#With res = 0.15, cells in cluster 3 do not have any significantly different genes that are deferentially regulated
#Therefore, we decided to proceed with res = 0.1
cx.raw <- FindClusters(cx.raw, resolution = 0.1)
#Culture plots
pdf("OS17.Culture.pdf", width = 5, height = 5)
pcx <- DimPlot(cx.raw, reduction = "umap", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS17 in Culture") + NoAxes()
pcx +  scale_color_npg(alpha = 1)
dev.off()

#cell cycle distribution of clusters

cid <- sort(unique(cx.raw@active.ident))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("OS17_pie.pdf", width = 15, height = 15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
for(i in 1:4){
  temp <- WhichCells(cx.raw, idents = cid[[i]]) 
  LT.cells <- subset(cx.raw, cells = temp)
  df <- table(LT.cells$Phase)
  bp <- ggplot(as.data.frame(df), aes(x="", y=Freq, fill=Var1))+
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0) + scale_fill_brewer(palette="Blues")+
    theme_minimal() + NoAxes() +
    ggtitle(paste("Cluster", cid[[i]]))
  print(pie, vp=vplayout(ceiling(i/4), i))
}
dev.off()

#Total number of cells
# table(Idents(cx.raw)) %>%sum()
# [1] 3178

# NCHOS7
# Create Seurat objects and perform initial QC. Label original source.

S0043 <- tenXLoadQC(path10x = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0043-NCHOS7-flank/filtered_feature_bc_matrix/", spec = "mixHuman")
S0043 <- subset(S0043, subset = nCount_RNA <50000 & percent.mt <25)
VlnPlot(S0043, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

S0043 <- NormalizeData(S0043) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(pc.genes = S0043@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

DimPlot(S0043, reduction = "umap", pt.size = 1, label = T) + 
  coord_fixed() + 
  scale_color_npg(alpha = 0.7)

# Attempt to regress out the effects of cell cycle on these tumor cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
S0043 <- CellCycleScoring(object = S0043, s.features = s.genes,
                          g2m.features = g2m.genes, set.ident = TRUE)

S0043 <- ScaleData(object = S0043, vars.to.regress = c("S.Score", "G2M.Score"),
                   features = rownames(x = S0043)) 
S0043 <- RunPCA(S0043, pc.genes = S0043@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

#Find optimal clustering resolution
S0043.CCR.nC <- nRes(S0043.CCR, 
                     res = seq(from = 0.1, to = 0.15, by = 0.01))
plot <- pSil(S0043.CCR.nC, 0.15)

S0043.CCR <- FindClusters(S0043.CCR, resolution = 0.15)

save(S0043.CCR, file ="R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure1/NCHOS7-S0043_CCR_v1.RData")

#Culture plots
pdf("NCHOS7.Culture.pdf", width = 5, height = 5)
pcx <- DimPlot(S0043.CCR, reduction = "umap", pt.size = 1, label = F) + 
  coord_fixed() + 
  ggtitle("OS17 in Culture") + NoAxes()
pcx +  scale_color_npg(alpha = 1)
dev.off()

#cell cycle distribution of clusters

cid <- sort(unique(S0043.CCR@active.ident))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("NCHOS7_pie.pdf", width = 15, height = 15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
for(i in 1:4){
  temp <- WhichCells(S0043.CCR, idents = cid[[i]]) 
  LT.cells <- subset(S0043.CCR, cells = temp)
  df <- table(LT.cells$Phase)
  bp <- ggplot(as.data.frame(df), aes(x="", y=Freq, fill=Var1))+
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0) + scale_fill_brewer(palette="Blues")+
    theme_minimal() + NoAxes() + NoLegend()
  ggtitle(paste("Cluster", cid[[i]]))
  print(pie, vp=vplayout(ceiling(i/4), i))
}
dev.off()

#Total number of cells
# table(Idents(S0043.CCR)) %>%sum()
# [1] 1998

# library("scales")
# show_col(pal_nejm("default")(8))
# show_col(pal_nejm("default", alpha = 0.6)(8))

##########################Pathway enrichment analysis; display adjust p-values
B.list <- list(OS17 = cx.raw,
               OS7 = S0043.CCR)

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

cx.pd <- full_join(temp1, temp2, by = "rn")

##Write table to save for supplemental data and for generating Heatmap
write.table(cx.pd, file=('R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2021/Fig1.tsv'), quote=FALSE, sep='\t')
save.image(file ="R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2021/ProcessedData/CulturePDXHetero.RData")

cx.pd <- column_to_rownames(cx.pd, var = "rn")
cx.pd[is.na(cx.pd)] <- 1

#remove "HALLMARK_"
c <- rownames(cx.pd)
c <- gsub("HALLMARK_", "", c)
#remove "_" by removing special characters
c <- gsub("_", " ", c)
rownames(cx.pd) <- c

cx.pd.log <- -log10(cx.pd)
colnames(cx.pd.log) <- c("A0", "A1", "A2", "B0", "B1", "B2", "B3")
##To prevent values that show up as -0.00, we take the absolute values
cx.pd.log <- abs(cx.pd.log)

cx.pd.log %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  left_join(aka4 %>% rownames_to_column(var = "Sample")) %>%
  pivot_longer(c(-Sample, -Lpercent),
               names_to = "Pathway",
               values_to = "pval",
               names_repair = "minimal") %>%
  mutate(Signif = pval >= (-1 * log10(0.05))) %>%
  ggplot(., aes(x = Sample, y = Pathway, fill = Signif)) + 
  geom_tile() +
  geom_text(aes(label = sprintf("%0.2f", pval))) +
  scale_fill_manual(values = c("gray", "red")) +
  scale_y_discrete(limits = rev) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("") 
