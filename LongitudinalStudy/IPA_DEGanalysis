#IPA analysis
# Decided to proceed with A2, B1, (B2, B3), C1, D0 and D3
# Find differentially expressed features selected cluster and remaining cells
#OS17
A0.markers <- FindMarkers(OS17.TL, ident.1 = "0", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
A1.markers <- FindMarkers(OS17.TL, ident.1 = "1", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
A2.markers <- FindMarkers(OS17.TL, ident.1 = "2", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
A3.markers <- FindMarkers(OS17.TL, ident.1 = "3", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)

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
