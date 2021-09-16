####################subset to cells that have a lineage tag
memory.limit(size = 40000)
#Culture
Culture.LTcells <- vector(mode = "character")
for(i in 1:nrow(cx.lt.list)){
  temp <- WhichCells(cx.sub, expression = lt == cx.lt.list[[i,1]]) 
  Culture.LTcells <- append(Culture.LTcells, temp)
}
cx.tagged <- subset(cx.sub, cells = Culture.LTcells)

#Tibia
tib.LTcells <- vector(mode = "character")
for(i in 1:nrow(tib.lt.list)){
  temp <- WhichCells(tib.sub, expression = lt == tib.lt.list[[i,1]]) 
  tib.LTcells <- append(tib.LTcells, temp)
}
tib.tagged <- subset(tib.sub, cells = tib.LTcells)

#Lung
lung.LTcells <- vector(mode = "character")
for(i in 1:515){
  temp <- WhichCells(lung.sub, expression = lt == lung.lt.list[[i,1]]) 
  lung.LTcells <- append(lung.LTcells, temp)
}
lung.tagged <- subset(lung.sub, cells = lung.LTcells)

lung.sub <- subset(lung.sub, idents = c("0", "1", "2", "6")) ##Remove outlier cells in cluster 3 and 5

pdf("LungOutline.pdf", height = 5, width = 5)
DimPlot(lung.sub, reduction = "umap", pt.size = 1, label = TRUE) + 
  coord_fixed() + 
  ggtitle("Lung-derived") +
  xlim(0,9) + NoAxes()+ 
  scale_color_npg(alpha = 1)
dev.off()

#Lung
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("grid_lung_v4.pdf", width = 10, height = 10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
for (i in 1:4) {
  p <- DimPlot(lung.tagged, 
               reduction = "umap",
               pt.size = 1,
               cells.highlight = WhichCells(lung.sub, expression = lt == lung.lt.list[[i,1]]), 
               cols.highlight = "#4DBBD5FF",
               sizes.highlight = 3) + 
    coord_fixed() +
    theme(legend.position = "none") +
    ggtitle(paste("Clone", lung.lt.list[[i,1]])) +
    xlim(0,9) + NoLegend() + NoAxes()
  print(p, vp=vplayout(ceiling(i/4), i))
}
dev.off()

#cell cycle distribution of top 4 lineages - change to Pie chart and percentage

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("pie.pdf", width = 15, height = 15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
for(i in 1:4){
  temp <- WhichCells(lung.sub, expression = lt == lung.lt.list[[i,1]]) 
  LT.cells <- subset(lung.sub, cells = temp)
  df <- table(LT.cells$Phase)
  bp <- ggplot(as.data.frame(df), aes(x="", y=Freq, fill=Var1))+
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0) + scale_fill_brewer(palette="Blues")+
    theme_minimal() + NoAxes() +
    ggtitle(paste("Clone", lung.lt.list[[i,1]]))
  print(pie, vp=vplayout(ceiling(i/4), i))
}
dev.off()

#barplot
bp <-  ggplot(as.data.frame(df), aes(x = Var1, y = Freq)) +
  geom_bar(fill = "#4DBBD5FF", stat = "identity") +
  ggtitle("Cell cycle distribution") +
  ylab("Count") +
  xlab("Phase")

#######################Pathways activated in top 10 lineages

lung.lt.list[1:10,1]
# 052-082154 
# 039-074146 
# 037-078792 
# 039-052640 
# 095-044696 
# 039-098388 
# 064-065939 
# 013-044371 
# 063-007191 
# 018-082025

data <- lung.tagged

#Subset cells of each lineage into a separate cluster in lung tagged
for(i in 1:10) {
  cells <- WhichCells(data, expression = lt == lung.lt.list[[i,1]])
  Idents(object = data, cells = cells) <- paste("Enriched Lineage", i, sep = " ")
}

library(forcats)
Idents(data) <- fct_rev(Idents(data))

# Idents(data) <- data$seurat_clusters
# Stash cell identity classes
data[["Lineage"]] <- Idents(data)

DimPlot(data)
#reorganize Lineage levels


#Top markers for each lineage

# find markers for every cluster compared to all remaining cells, report positive and negative ones
data.markers <- FindAllMarkers(data, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers <- FindAllMarkers(data, only.pos = TRUE)

#Run through function
em.hm <- DGEA(data)

#Remove non-LT clusters
em.hm.lt <- em.hm[-c(1:4)]
clust.ids = sort(unique(data@active.ident))
clust.ids <- clust.ids[-c(1:4)]
c <- rownames(em.hm.lt)
#remove "HALLMARK_"
c <- gsub("HALLMARK_", "", c)
#remove "_" by removing special characters
c <- gsub("_", " ", c)
rownames(em.hm.lt) <- c
#remove rows with low pvalues
row_names_df_to_remove<-c("FATTY ACID METABOLISM",
                          "ADIPOGENESIS",
                          "PEROXISOME",
                          "SPERMATOGENESIS",
                          "XENOBIOTIC METABOLISM",
                          "PI3K AKT MTOR SIGNALING",
                          "HEDGEHOG SIGNALING")
em.hm.lt <- em.hm.lt[!(row.names(em.hm.lt) %in% row_names_df_to_remove),]

col.breaks=seq(-log10(1),min(max(-log10(em.hm.lt))+1,20),by=0.5)
col=inferno(length(col.breaks)) # library(viridis)
col=c("white",colorRampPalette(brewer.pal(n = 7, name ="Reds"))(50))
p <- pheatmap(-log10(em.hm.lt[,1:length(clust.ids)]),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 10,
         treeheight_row = 0,
         treeheight_col=0,
         color = col,
         scale='none',
         breaks=col.breaks,
         fontsize = 11)

pdf("EnrichedLTHeatmap.pdf", width = 9, height = 10)
print(p)
dev.off()

########################Run IPA on the DEGs from top ten lineages
#NOTE: After Enriched lineage 5, not more than 25 DEGs with P_adj_val < 0.01

# L1.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 1", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L2.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 2", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L3.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 3", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L4.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 4", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L5.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 5", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L6.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 6", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L7.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 7", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L8.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 8", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L9.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 9", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# L10.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 10", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
# 
# #remove unnecessary columns
# myList <- list(L1.markers, L2.markers, L3.markers, L4.markers, L5.markers, L6.markers, L7.markers, L8.markers, L9.markers, L10.markers)
# myList <- lapply(myList, function(x) { x["p_val"] <- NULL; x[, 2:3] <- NULL; x })
# 
# #write markers to excel sheet
# names <- c("L1.markers", "L2.markers", "L3.markers", "L4.markers", "L5.markers", "L6.markers", "L7.markers", "L8.markers", "L9.markers", "L10.markers")
# library(xlsx)
# for (i in 1:10){
#   location <- paste0("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/L",i,".markers.xlsx", sep = "")
#   write.xlsx(myList[[i]], file = location, 
#              col.names = TRUE, row.names = TRUE, append = FALSE)
# } 

#######################Gllycolysis and Hypoxia in LT enriched cells
load("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/Preprocessed_Seurat_objects/OS17_Zymo_CTL.RData")
load("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/OS17CTL.sub.RData")

DimPlot(lung.tagged, reduction = "umap", pt.size = 1, label = TRUE) + 
  coord_fixed() + 
  ggtitle("Lung-derived") +
  xlim(0,9) + NoAxes()+ 
  scale_color_npg(alpha = 1)

data <- lung.tagged

#Apply Glycolysis module score to each of the primary tumor cells
gl <- list()
Glycolysis <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/Glycolysis.txt", header = FALSE, sep = ",")
Glycolysis <- Glycolysis[,1]
OxPhos <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/OxPhos.txt", header = FALSE, sep = ",")
OxPhos <- OxPhos[,1]
Hypoxia <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/Hypoxia.txt", header = FALSE, sep = ",")
Hypoxia <- Hypoxia[,1]
SASP <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/SASP.txt", header = FALSE, sep = ",")
SASP <- SASP[,1]
EMT <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/EMT.txt", header = FALSE, sep = ",")
EMT <- EMT[,1]
TNFA <- read.table(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/TNFA.txt", header = FALSE, sep = ",")
TNFA <- TNFA[,1]
  
gl$Glycolysis <- Glycolysis
gl$OxPhos <- OxPhos
gl$Hypoxia <- Hypoxia
gl$SASP <- SASP
gl$EMT <- EMT
gl$TNFA <- TNFA
mod_names <- names(gl)

#Added min.cutoff to isolate for cells most module like
for(i in 1:length(mod_names)) {
  data <- AddModuleScore(data, gl[i], name = names(gl[i]))
  p <- FeaturePlot(data, features = str_c(names(gl[i]), "1"), min.cutoff = 0.1, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  pdf(paste0("FeaturePlot ",names(gl[i]), ".pdf"))
  print(p)
  dev.off()
}

#Compare Enriched vs non-enriched
#Glycolysis and Hypoxia seem to be trending together

#Group first 12 against the tagged Lung
data <- lung.tagged
LT <- vector(mode = "character")
for(i in 1:10){
  temp <- WhichCells(lung.sub, expression = lt == lung.lt.list[[i,1]]) 
  LT <- append(LT, temp)
}
Idents(data, cells = LT) <- 'Enriched'

#Non-enriched cells are clones with a frequency of 1
# lung.lt.list[134,] 
# 2
# lung.lt.list[134,] 
# 1
LT <- vector(mode = "character")
for(i in 135:377){
  temp <- WhichCells(lung.sub, expression = lt == lung.lt.list[[i,1]])
  LT <- append(LT, temp)
}
Idents(data, cells = LT) <- 'Non-enriched'

cell.ids <- sample(colnames(data))

#Culture plots
# pdf("OS17_Enr.pdf", width = 5, height = 5)
DimPlot(data, 
             reduction = "umap",
             pt.size = 1,
             cells.highlight = WhichCells(data, idents = "Enriched"), 
             cols.highlight = "#E64B35FF",
             sizes.highlight = 2, order = cell.ids) + 
  coord_fixed() 
# dev.off()
####################Pathways associated with Enriched clonal families
library(tibble)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(stringr)

B.list <- list(OS17 = data)

em.hm.list <- list()
for (i in 1:length(B.list)) {
  em.hm.list[[i]] <- DGEA(B.list[[i]]) 
}

temp1 <- em.hm.list[[1]]
temp1 <- -log10(temp1)
#remove "_" by removing special characters
c <- rownames(temp1)
c <- gsub("HALLMARK ", "", c)
c <- gsub("_", " ", c)
rownames(temp1) <- c
#remove negative 0 values due to number of decimal points displayed
temp1 <- abs(temp1)
temp1 <- temp1["Enriched"]

#Subset to only pathways with siginificant pvalue (-1 * log10(0.05)) = 1.30103
temp1["Enriched"] >= 1.30103 -> temp1$logical
colnames(temp1) <- c("pvalue", "significant")
logical <- temp1$significant
temp1 <- temp1[logical,]
temp1 <- rownames_to_column(temp1, var = "Pathway")
position <- order(temp1$pvalue, decreasing = TRUE)
temp1 <- temp1[position,]
plot <- ggplot(temp1, aes(x = pvalue,
                         y = Pathway)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_bw() +
  ylab("") 

pdf("Enr_Pathway.pdf", height = 5, width = 10)
plot(plot)
dev.off()

# heatmap <- temp1 %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "Sample") %>%
#   left_join(aka4 %>% rownames_to_column(var = "Sample")) %>%
#   pivot_longer(c(-Sample, -Lpercent),
#                names_to = "Pathway",
#                values_to = "pval",
#                names_repair = "minimal") %>%
#   mutate(Signif = pval >= (-1 * log10(0.05))) %>%
#   ggplot(., aes(x = Sample, y = Pathway, fill = Signif)) + 
#   geom_tile() +
#   geom_text(aes(label = sprintf("%0.2f", pval))) +
#   scale_fill_manual(values = c("gray", "red")) +
#   scale_y_discrete(limits = rev) +
#   theme(legend.position = "none",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("") 
# 
# pdf("Enr_Pathway.pdf", height = 8, width = 8)
# plot(heatmap)
# dev.off()

# Find differentially expressed features between Enriched and remaining Lung
enr.markers <- FindMarkers(data, 
                           ident.1 = "Enriched",
                           logfc.threshold = 0.25, 
                           only.pos = TRUE, 
                           test.type= "DESeq2")

lung.markers <- FindMarkers(data, 
                            ident.1 = "Non-enriched", 
                            logfc.threshold = 0.25, 
                            only.pos = TRUE, 
                            test.type= "DESeq2")


# view results
write.xlsx(enr.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/enr.10_markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(lung.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/lung_markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

VlnPlot(object = data, idents = c('Enriched', 'Non-enriched'), features = c('ENO2','PGK1','ENO1'), combine = TRUE, ncol = 3)
VlnPlot(object = data, idents = c('Enriched', 'Non-enriched'), features = c('IFITM3','ISG15','IFI6'), combine = TRUE, ncol = 3)

