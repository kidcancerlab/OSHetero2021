#Subset lung cells from merged datasets; identify cell, DEGs compared to rest of the lung to identify Upstream regulators

#Subset lung cells from OS17.TL
OS17.Lung <- subset(OS17.TL, subset = src == "Lung")
t143b.Lung <- subset(t143b, subset = src == "Lung")
OS2.Lung <- subset(OS2, subset = src == "Lung")
OS7.Lung <- subset(OS7, subset = src == "Lung")
  
# Use colors from the npg palatte in ggsci
list <- list(OS17.Lung, t143b.Lung, OS2.Lung, OS7.Lung)
names(list) <- c("OS17.Lung", "t143b.Lung", "OS2.Lung", "OS7.Lung")

for(i in 1:length(list)) {
  p1 <- DimPlot(list[[i]], reduction = "umap", pt.size = 1, label = T) + 
    coord_fixed() + 
    ggtitle("Lung from merged TL") + 
    scale_color_npg(alpha = 1)
  print(p1)
}

#DEGs in defined clusters
OS17.Lung.markers2 <- FindMarkers(OS17.Lung, 
                           ident.1 = "2", 
                           ident.2 = NULL, 
                           min.pct = 0.25)

write.xlsx(OS17.Lung.markers2, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure2/OS17.Lung.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

t143b.Lung.markers1 <- FindMarkers(t143b.Lung, 
                                 ident.1 = "1", 
                                 ident.2 = NULL, 
                                 min.pct = 0.25)

t143b.Lung.markers2 <- FindMarkers(t143b.Lung, 
                                  ident.1 = "2", 
                                  ident.2 = NULL, 
                                  min.pct = 0.25)

t143b.Lung.markers3 <- FindMarkers(t143b.Lung, 
                                  ident.1 = "3", 
                                  ident.2 = NULL, 
                                  min.pct = 0.25)
write.xlsx(t143b.Lung.markers1, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure2/t143b.Lung.markers1.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(t143b.Lung.markers2, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure2/t143b.Lung.markers2.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(t143b.Lung.markers3, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure2/t143b.Lung.markers3.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

OS2.Lung.markers1 <- FindMarkers(OS2.Lung, 
                                 ident.1 = "1", 
                                 ident.2 = NULL, 
                                 min.pct = 0.25)

write.xlsx(OS2.Lung.markers1, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure2/OS2.Lung.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

OS7.Lung.markers0 <- FindMarkers(OS7.Lung, 
                                ident.1 = "0", 
                                ident.2 = NULL, 
                                min.pct = 0.25)

OS7.Lung.markers3 <- FindMarkers(OS7.Lung, 
                                ident.1 = "3", 
                                ident.2 = NULL, 
                                min.pct = 0.25)

write.xlsx(OS7.Lung.markers0, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure2/OS7.Lung.markers0.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(OS7.Lung.markers3, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure2/OS7.Lung.markers3.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
#Perform IPA - Upsstream analysis
#Heatmap of IPA comparison file 
df <- read.delim("Fig3_UpReg.txt", header = TRUE)
rownames(df) <- df[,1]
df[,1] <- NULL
df[is.na(df)]=0
library(pheatmap)   
df_transpose <- as.data.frame(t(df))
df <- as.matrix(df_transpose)

subj <- rownames(df)
aka2 = data.frame(ID = factor(c(rep("OS17", 1),
                                rep("t143B", 3),
                                rep("OS2", 1),
                                rep("OS7", 2))))
rownames(aka2)<-subj
aka3 = list(ID = c(OS17 = "#E64B35FF", t143B = "#4DBBD5FF", OS2 = "#00A087FF", OS7 = "#3C5488FF"))

pheatmap(t((df)),
         color = rev(brewer.pal(7,"RdBu")), 
         annotation_col = aka2, 
         annotation_colors = aka3,
         annotation_legend = TRUE,
         show_colnames = T, show_rownames = T, cluster_rows = F, 
         cluster_cols = T, legend = TRUE,
         clustering_distance_rows = "euclidean", border_color = "white")
