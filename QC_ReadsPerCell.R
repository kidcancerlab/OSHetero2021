load(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/OS.listv3.CCR.RData")

name <- names(OS.list)
for (i in 1:13){
  data <- OS.list[[i]]
  data <- Seurat::GetAssayData(data) 
  p <- ggplot(tibble(readsPerCell = colSums(data)), 
              aes(x = readsPerCell)) + 
    geom_histogram(bins = 200) +
    ggtitle(paste0(name[i]))
  pdf(paste0(name[i], " Reads.pdf", sep =" "))
  plot(p)
  dev.off()
}

######################Total variance explained by PCs

##load only scaled data
# load(file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/OS.list_PC50.RData")

# # Determine percent of variation associated with each PC
# pct <- data@reductions$pca@stdev / sum(data@reductions$pca@stdev) * 100
# 
# # Calculate cumulative percents for each PC
# cum <- cumsum(pct)
# 
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
# co1 <- which(cum > 90 & pct < 5)[1]
# co1
# 
# # Determine the difference between variation of PC and subsequent PC
# co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
# co2


Cummulative_Var <- vector()
for (i in 1:length(OS.list)){
  data <- OS.list[[i]]
  pdf(paste0("ElbowPlot_", name[i], ".pdf"))
  p <- ElbowPlot(data, ndims = 50)
  print(p)
  dev.off()
  pct <- data@reductions$pca@stdev / sum(data@reductions$pca@stdev) * 100
  Cummulative_Var <- append(Cummulative_Var, sum(pct[1:20]))
}

min(Cummulative_Var)
max(Cummulative_Var)

# For each dataset, the first twenty principal components were selected based on the elbow plot 
# for percentage explained variances, representing ~55.5–60.3% of total variances. 

ElbowPlot <- list()
for(i in 1:length(OS.list)){
  data <- OS.list[[i]]
  ElbowPlot[[i]] <- ElbowPlot(data, ndims = 50)
}


pdf("ElbowPlot.pdf", width = 15, height = 15)
grid.arrange(ElbowPlot[[1]], ElbowPlot[[2]], ElbowPlot[[3]], ElbowPlot[[4]],
             ElbowPlot[[5]], ElbowPlot[[6]], ElbowPlot[[7]], ElbowPlot[[8]],
             ElbowPlot[[9]], ElbowPlot[[10]], ElbowPlot[[11]], ElbowPlot[[12]],
             ElbowPlot[[13]], nrow = 4)
dev.off()

vargenes <- vector()
for (i in 1:length(OS.list)){
  data <- OS.list[[i]]
  vargenes <- append(vargenes, length(data@assays$RNA@var.features))
}


####Patient data
library(qs)
tumor_comb <- qread("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/tar_dir/tumor_comb.qs")

Cummulative_Var <- vector()
for (i in 1:length(raw)){
  data <- raw[[i]] %>% RunPCA()
  pct <- data@reductions$pca@stdev / sum(data@reductions$pca@stdev) * 100
  Cummulative_Var <- append(Cummulative_Var, sum(pct[1:20]))
}

min(Cummulative_Var)
max(Cummulative_Var)
