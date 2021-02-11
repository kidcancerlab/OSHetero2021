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


cx <- full_join(temp1, temp2, by = "rn")
pd <- full_join(temp3, temp4, by = "rn")
cx.pd <- full_join(cx, pd, by = "rn")

rownames(cx.pd)=cx.pd[,1]
cx.pd=cx.pd[,-1]
cx.pd[is.na(cx.pd)]=1
cx.pd <- cx.pd[!grepl('NA', rownames(cx.pd)), ] #remove rows with rownames 'NA'

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

cols <- colorRampPalette(c("snow1", "red"), # distances 3 to max(distmat) colored from green to black
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

col <- colorRampPalette(c("snow1", "steelblue4"), # distances 3 to max(distmat) colored from green to black
                             100)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)
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
Heatmap(t(scale(df)), 
        col = color, 
        name = "-log10(p)", 
        top_annotation = column_ha)
dev.off()
