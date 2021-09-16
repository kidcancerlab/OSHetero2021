############################Overlay osteoblast markers on Culture/Flank cells 
load("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure3/OB-S0031_CCR_v2.RData")
library(tibble)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(stringr)

data <- S0031.CCR
DimPlot(data, reduction = "umap", pt.size = 1, label = T) + 
  coord_fixed() + 
  scale_color_npg(alpha = 0.7)

#Set optimum resolution for clustering
data.nC <- nRes(data, 
                res = seq(from = 0.1, to = 0.9, by = 0.1))
#Optimal resolution identified as 0.1
plot <- pSil(data.nC, 0.9)
data <- data %>% FindClusters(resolution = 0.05)

#Culture plots
pdf("Osteoblast_Cx.pdf", width = 5, height = 5)
pcx <- DimPlot(data, reduction = "umap", pt.size = 1, label = F) + 
  coord_fixed() + 
  ggtitle("OS17 in Culture") + NoAxes()
pcx +  scale_color_npg(alpha = 1)
dev.off()

#cell cycle distribution of clusters

cid <- sort(unique(data@active.ident))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("OB_pie.pdf", width = 15, height = 15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
for(i in 1:2){
  temp <- WhichCells(data, idents = cid[[i]]) 
  LT.cells <- subset(data, cells = temp)
  df <- table(LT.cells$Phase)
  bp <- ggplot(as.data.frame(df), aes(x="", y=Freq, fill=Var1))+
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0) + scale_fill_brewer(palette="Blues")+
    theme_minimal() + NoAxes() +
    ggtitle(paste("Cluster", cid[[i]]))
  print(pie, vp=vplayout(ceiling(i/2), i))
}
dev.off()


###################Heatmap of OB markers
ob.markers <- FindAllMarkers(data, only.pos = FALSE, min.pct = 0.25)

ob.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

pdf("OBHeatmap.pdf", width = 7, height = 6)
p <- DoHeatmap(data, features = top10$gene, size = 5) + NoLegend()
p + theme(text = element_text(size=16, color = "black"))
dev.off()

####################Pathways associated with OB clusters
B.list <- list(OB = data)

em.hm.list <- list()
for (i in 1:length(B.list)) {
  em.hm.list[[i]] <- DGEA(B.list[[i]]) 
}

temp1 <- em.hm.list[[1]]
temp1 <- -log10(temp1)
#remove "_" by removing special characters
c <- rownames(temp1)
c <- gsub("SIG ", "", c)
c <- gsub("NABA ", "", c)
c <- gsub("SA ", "", c)
c <- gsub("_", " ", c)
rownames(temp1) <- c
#remove negative 0 values due to number of decimal points displayed
temp1 <- abs(temp1)

heatmap <- temp1 %>%
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

pdf("OB_Pathway.pdf", height = 8, width = 8)
plot(heatmap)
dev.off()

#############################Generate phentoype module
cluster.0 <- FindMarkers(data, ident.1 = "0", only.pos = TRUE, min.pct = 0.25)
cluster.1 <- FindMarkers(data, ident.1 = "1", only.pos = TRUE, min.pct = 0.25)
cluster.2 <- FindMarkers(data, ident.1 = "2", only.pos = TRUE, min.pct = 0.25)

#Extract gene names
cluster0 <- list(rownames(cluster.0))
cluster1 <- list(rownames(cluster.1))
cluster2 <- list(rownames(cluster.2))

write.xlsx(cluster.0, file = "C:/Users/rssxr002/Downloads/cluster.0.xlsx")
write.xlsx(cluster.1, file = "C:/Users/rssxr002/Downloads/cluster.1.xlsx")
write.xlsx(cluster.2, file = "C:/Users/rssxr002/Downloads/cluster.2.xlsx")

##############Annotation of phenotype module
#Ran msigdb gene set enrichment on http://www.gsea-msigdb.org/gsea/index.jsp
#Visualization of p values usign barplot

SanjanaPlots2 <- read.delim("C:/Users/rssxr002/Downloads/SanjanaPlots2.txt")

ggplot(SanjanaPlots2, aes(x = -1 * Order,
                         y = FDR.log,
                         fill = FDR.log > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Cluster, ncol = 1) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  geom_text(aes(y = Label_y * 15, label = Pathway)) +
  theme_bw() +
  ylab("") 

#Addmodule score to OS cell culture/flank models
load("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/OS.listv3.CCR.RData")

# OB <- OS.list[[1]]
# OS17.cx <- OS.list[[2]]
# OS17.Tibia <-  OS.list[[3]]
# OS17.Lung <-OS.list[[4]]
# t143B.cx <- OS.list[[5]]
# t143B.Tibia <-OS.list[[6]]
# t143B.Lung <- OS.list[[7]]
# OS2.Flank <- OS.list[[8]]
# OS2.Tibia <- OS.list[[9]]
# OS2.Lung <- OS.list[[10]]
# OS7.Flank <-  OS.list[[11]]
# OS7.Tibia <- OS.list[[12]]
# OS7.Lung <- OS.list[[13]]

for (i in 1:13) {
  OS.list[[i]] <- AddModuleScore(OS.list[[i]], features = cluster0, ctrl = 5, name = 'cluster0')
  OS.list[[i]] <- AddModuleScore(OS.list[[i]], features = cluster1, ctrl = 5, name = 'cluster1')
  OS.list[[i]] <- AddModuleScore(OS.list[[i]], features = cluster2, ctrl = 5, name = 'cluster2')
}

plot1 <- DotPlot(OS.list[[2]], features = c("cluster01", "cluster11", "cluster21")) + coord_flip() 
plot2 <- DotPlot(OS.list[[3]], features = c("cluster01", "cluster11", "cluster21")) + coord_flip()
plot3 <- DotPlot(OS.list[[4]], features = c("cluster01", "cluster11", "cluster21")) + coord_flip()

# DotPlot(OS.list[[1]], features = c("cluster01", "cluster11", "cluster21")) + coord_flip() 

# plot1 <- RidgePlot(OS.list[[11]], features = c("cluster01", "cluster11", "cluster21"))
# plot2 <- RidgePlot(OS.list[[12]], features = c("cluster01", "cluster11", "cluster21"))
# plot3 <- RidgePlot(OS.list[[13]], features = c("cluster01", "cluster11", "cluster21"))

  CombinePlots(
    plots = list(plot1, plot2, plot3),
    legend = 'none', nrow = 1)
  
#####################Optimize cluster resolution 
  data <- OS.list[[2]] 
  DimPlot(data, reduction = "umap", pt.size = 1, label = T) + 
    coord_fixed() + 
    scale_color_npg(alpha = 0.7)
  
  #Set optimum resolution for clustering
  data.nC <- nRes(data, 
                  res = seq(from = 0.1, to = 0.3, by = 0.05))
  plot <- pSil(data.nC, 0.2)
  OS.list[[4]] <- OS.list[[2]] %>% FindClusters(resolution = 0.2)
  
  # OS.list[[1]] res 0.15
  # OS.list[[2]] res 0.15
  # OS.list[[3]] res 0.25
  # OS.list[[4]] res 0.2
  
  # OS.list[[5]] res 0.25
  # OS.list[[6]] res 0.25
  # OS.list[[7]] res 0.25
  # OS.list[[8]] res 0.25
  # OS.list[[9]] res 0.25
  # OS.list[[10]] res 0.25
  # OS.list[[11]] res 0.25
  # OS.list[[12]] res 0.25
  # OS.list[[13]] res 0.25
  
  
  OS.list[[1]] <- OS.list[[1]] %>% FindClusters(resolution = 0.15)
  OS.list[[2]] <- OS.list[[2]] %>% FindClusters(resolution = 0.15)
  OS.list[[3]] <- OS.list[[3]] %>% FindClusters(resolution = 0.25)
  OS.list[[4]] <- OS.list[[4]] %>% FindClusters(resolution = 0.2)
  # OS.list[[5]] <- OS.list[[5]] %>% FindClusters(resolution = 0.2)
  # OS.list[[6]] <- OS.list[[6]] %>% FindClusters(resolution = 0.2)
  # OS.list[[7]] <- OS.list[[7]] %>% FindClusters(resolution = 0.2)
  # OS.list[[8]] <- OS.list[[8╒]] %>% FindClusters(resolution = 0.2)
  # OS.list[[9]] <- OS.list[[9]] %>% FindClusters(resolution = 0.2)
  # OS.list[[10]] <- OS.list[[10]] %>% FindClusters(resolution = 0.2)
  # OS.list[[11]] <- OS.list[[11]] %>% FindClusters(resolution = 0.2)  
  # OS.list[[12]] <- OS.list[[12]] %>% FindClusters(resolution = 0.2)
  # OS.list[[13]] <- OS.list[[13]] %>% FindClusters(resolution = 0.2)