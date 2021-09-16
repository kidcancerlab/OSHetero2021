# Extract the barcode frequency lists
cx.lt.list <- processLTBC(cx.raw,
                          lt.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0016-zymo/lt.fq",
                          cid.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0016-zymo/cid.fq",
                          ret.list = T)
tib.lt.list <- processLTBC(tib.raw,
                           lt.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0018xS0028/lt.fq",
                           cid.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0018xS0028/cid.fq",
                           ret.list = T)
lung.lt.list <- processLTBC(lung.raw,
                            lt.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0024xS0029/lt.fq",
                            cid.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0024XS0029/cid.fq",
                            ret.list = T)

# Subset the three samples, but retain the current umap data
cx.sub <- subset(os17, subset = src == "Culture")
tib.sub <- subset(os17, subset = src == "Tibia")
lung.sub <- subset(os17, subset = src == "Lung")

# Generate cluster identification plots for the top 4 enriched clones from each sample
# Use colors from the npg palatte in ggsci
library("grid")
#Culture
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("grid_cx_v1.pdf", width = 10, height = 10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
for (i in 1:4) {
  p <- DimPlot(cx.tagged, 
               reduction = "umap",
               pt.size = 1,
               cells.highlight = WhichCells(cx.sub, expression = lt == cx.lt.list[[i,1]]), 
               cols.highlight = "#E64B35FF",
               sizes.highlight = 3, order = cell.ids) + 
    coord_fixed() +
    theme(legend.position = "none") +
    ggtitle(paste("Clone", cx.lt.list[[i,1]])) +
    xlim(-12,-6) + NoLegend() + NoAxes()
  print(p, vp=vplayout(ceiling(i/4), i))
}
dev.off()

#Tibia
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("grid_tib_v1.pdf", width = 10, height = 10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
for (i in 1:4) {
  p <- DimPlot(tib.tagged, 
               reduction = "umap",
               pt.size = 1,
               cells.highlight = WhichCells(tib.sub, expression = lt == tib.lt.list[[i,1]]), 
               cols.highlight = "#00A087FF",
               sizes.highlight = 3) + 
    coord_fixed() +
    theme(legend.position = "none") +
    ggtitle(paste("Clone", tib.lt.list[[i,1]])) +
    xlim(0,9) + NoLegend() + NoAxes()
  print(p, vp=vplayout(ceiling(i/4), i))
}
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

#Culture plots
pdf("Culture.pdf", width = 2.3, height = 2.3)
pcx <- DimPlot(cx.sub, reduction = "umap", pt.size = 1, label = F, group.by = "src") + 
  coord_fixed() + 
  ggtitle("Plate-derived") +
  xlim(-12,-6) + NoLegend() + NoAxes()
pcx + scale_color_manual(values = c("#E64B35FF"))
dev.off()

pdf("Tibia.pdf", width = 2.3, height = 2.3)
ptib <- DimPlot(tib.sub, reduction = "umap", pt.size = 1, label = F, group.by = "src") + 
  coord_fixed() + 
  ggtitle("Bone-derived")+
  xlim(0,9) + NoLegend() + NoAxes() 
ptib + scale_color_manual(values = c("#00A087FF"))
dev.off()

pdf("Lung.pdf", width = 2.3, height = 2.3)
plung <- DimPlot(lung.sub, reduction = "umap", pt.size = 1, label = F, group.by = "src") + 
  coord_fixed() + 
  ggtitle("Lung-derived") +
  xlim(0,9) + NoLegend() + NoAxes()
plung + scale_color_manual(values = c("#4DBBD5FF"))
dev.off()

#Together plot
pdf("All.pdf", width = 7, height = 7)
DimPlot(os17, reduction = "umap", group.by = "src", pt.size = 1, label = T) + 
  coord_fixed() +  
  scale_color_npg(alpha = 1)
dev.off()

DimPlot(os17, reduction = "umap", group.by = "src", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS17 Clusters") + 
  scale_color_npg(alpha = 0.7)

