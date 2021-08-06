library(rrrSingleCellUtils)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringr)
# library(future)
library(harmony)

set.seed(888)
# plan("multisession", workers = 10)
# options(future.globals.maxSize= 10 * 24000 * 1024^2)

# Make a list of sample names
s <- c("BC5", "BC6", "BC10", "BC11", "BC16",
       "BC17", "BC20", "BC21", "BC22")
qc <- c(18000, 25000, 25000, 30000, 70000, 
        40000, 70000, 50000, 50000)
path <- c("Conventional", "Conventional", "Conventional",
          "Conventional", "Conventional", "Chondroblastic",
          "Chondroblastic", "Intraosseous", "Chondroblastic")
type <- c("Primary", "Primary", "Lung Met", "Primary", "Primary",
          "Lung Met", "Primary", "Primary", "Primary")

# Download, file, and extract the files from GEO
if(!dir.exists("PrimaryTumor/GSE152048")) {
  tar_dir <- "PrimaryTumor/GSE152048"
  dir.create(tar_dir)
  geo_pre <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152048/suppl/GSE152048_"
  for(i in 1:length(s)){
    gse_path <- str_c(geo_pre, s[i], ".matrix.tar.gz")
    tar_file <- str_c(tar_dir, "/", s[i], ".tar.gz ")
    download.file(gse_path, destfile = tar_file, method = "auto")
    untar(tar_file, exdir = tar_dir)
    file.remove(tar_file)
  }
}

# Create a vector that contains normalized Seurat objects for all samples
raw <- c()
for(i in 1:9) {
  x <- tenx_load_qc(str_c("PrimaryTumor/GSE152048/",
                          s[i], "/"))
  x <- subset(x, subset = nCount_RNA < qc[i] & percent.mt <13)
  x$src <- s[i]
  x$type <- type[i]
  x$path <- path[i]
  x <- x %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.3) %>%
    RunUMAP(dims = 1:20)
  print(DimPlot(x, pt.size = 1, label = T) +
          coord_fixed() +
          theme(legend.position = "none") +
          ggtitle(str_c(s[i], " basic clustering")))
  raw <- c(raw, x)
}
rm(x)

# Merge into one Seurat object and integrate with harmony
comb <- merge(raw[[1]], y = c(raw[[2]], raw[[3]], raw[[4]], raw[[5]], 
                              raw[[6]], raw[[7]], raw[[8]], raw[[9]]),
              # merge.data = T,
              add.cell.ids = s,
              project = "GSE152048")

comb <- comb %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F)

comb <- RunHarmony(comb, group.by.vars = "src")
comb <- RunUMAP(comb, reduction = "harmony", dims = 1:30)
comb <- comb %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters()

if(!dir.exists("PrimaryTumor/SavedFiles")) {
  dir.create("PrimaryTumor/SavedFiles")
}
# Save a stopping point - merged and harmony aligned
qs::qsave(comb, "PrimaryTumor/SavedFiles/comb.qs")
# Start from this stopping point
comb <- qs::qread("PrimaryTumor/SavedFiles/comb.qs")

DimPlot(comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("GSE152048 composite") +
  theme(legend.position = "none")

DimPlot(comb, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("GSE152048 composite") 

DimPlot(comb, reduction = "umap", group.by = "type") +
  coord_fixed() +
  ggtitle("GSE152048 composite") 


# Separate out the primary tumor lesions
primary <- subset(comb, subset = type == "Primary")

# Ensure that groups are not just dominated by cell cycle effects
#primary <- kill_cc(primary)
# 
# # Set a post-cc-regression stopping point
# saveRDS(primary, "primary-ccReg.rds")
# # Start from this stopping point
# primary <- readRDS("primary-ccReg.rds")

DimPlot(primary, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("GSE152048 composite - primaries") +
  theme(legend.position = "none")

DimPlot(primary, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("GSE152048 composite - primaries") 

# Calculate and show module scores from the paper
ms <- list()
ms$Osteoblastic <- c("RUNX2", "COL1A1", "CDH11", "IBSP")
ms$Chondroblastic <- c("SOX9", "ACAN", "PTH1R")
ms$Osteoclast <- c("ACP5", "CTSK", "MMP9")
ms$Myeloid <- c("CD74", "CD14", "FCGR3A")
ms$TCell <- c("CD3E", "IL7R", "CD8A", "CD4", "NKG7")
ms$NKCell <- c("NKG7", "GNLY")
ms$NKTCell <- c("NKG7", "GNLY", "CD3E")
ms$DCCell <- c("CD1C", "FCER1A", "CLEC9A", "CCR7", "CD14", "CD163")
ms$Fibroblast <- c("DCN", "COL1A1")
ms$Pericyte <- c("RGS5", "ACTA2")
ms$MSC <- c("MME", "THY1", "CXCL12", "SFRP2")
ms$Endothelial <- c("PECAM1", "VWF")
ms$Myoblast <- c("MYL1", "MYLPF")
ms$BCell <- c("MS4A1", "CD19", "JCHAIN")

mod_names <- names(ms)
#Added min.cutoff to isolate for cells most module like
for(i in 1:length(mod_names)) {
  primary <- AddModuleScore(primary, ms[i], name = names(ms[i]))
  p <- FeaturePlot(primary, features = str_c(names(ms[i]), "1"), min.cutoff = 2, pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  print(p)
}

#look for proliferation
primary<-kill_cc(primary)
DimPlot(primary, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("G1/S/G2M cell cycle analysis - primary tumor") +
  theme(legend.position = "none")


#Select Clusters which are tumors. Notice how one group is G2M only. 
#Maybe possible that this group is tumor and due to 0inflation there are transcripts which are not active. 
#Could be cycling back to other tumor states
Idents(primary) <- primary$seurat_clusters
DimPlot(primary, reduction = "umap", label = T, repel = T) +
coord_fixed() +
  ggtitle("GSE152048 composite clusters - primaries") +
  theme(legend.position = "none")
saveRDS(primary,'PrimaryTumor/SavedFiles/primary.rds')
primary<-readRDS('PrimaryTumor/SavedFiles/primary.rds')

#Re-cluster after subsetting to tumor cells
tumor<- subset(primary, idents=c(2,3,27,16,15,6,28,5,20,29,0,24,14,21,33,35,13,26,7,9,23))
tumor<-tumor%>%
  NormalizeData()%>%
  ScaleData()%>%
  FindVariableFeatures()
tumor<-tumor%>%
  RunPCA(features=VariableFeatures(tumor))%>%
  RunUMAP(dims=1:20)%>%
  FindNeighbors()%>%
  FindClusters(resolution=0.3)

#Plot results
DimPlot(tumor, pt.size=1, reduction="umap", label=T, group.by="src")+
  coord_fixed()+
  ggtitle("Reclustering of Tumor Cells")
#Plot results
DimPlot(tumor, pt.size=1, reduction="umap", label=T)+
  coord_fixed()+
  ggtitle("Reclustering of Tumor Cells")
#kill_cc adds metadata but interferes with analyses
tumor<-kill_cc(tumor)

Idents(tumor) <- tumor$seurat_clusters
DimPlot(tumor, pt.size=1, reduction="umap", label=T)+
  coord_fixed()+
  ggtitle("Reclustering of Tumor Cells")

#Added min.cutoff to isolate for cells most module like
for(i in 1:length(mod_names)) {
  tumor <- AddModuleScore(tumor, ms[i], name = names(ms[i]))
  p <- FeaturePlot(tumor, features = str_c(names(ms[i]), "1"), pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  print(p)
}
saveRDS(tumor,'PrimaryTumor/SavedFiles/tumor.rds')
tumor<-readRDS('PrimaryTumor/SavedFiles/tumor.rds')

#Re-Subset data further because 7 and 8 look like fibroblast and periocytes modules scored high
tumor_2<-subset(tumor, idents=c(1,2,3,4,5,6,9,10,11,12,13,14))

#Reculuster new tumor
tumor_2<-tumor_2%>%
  NormalizeData()%>%
  ScaleData()%>%
  FindVariableFeatures()
tumor_2<-tumor_2%>%
  RunPCA(features=VariableFeatures(tumor_2))%>%
  RunUMAP(dims=1:20)%>%
  FindNeighbors()%>%
  FindClusters(resolution=0.3)

#Plot results
DimPlot(tumor_2, pt.size=.1, reduction="umap", label=T, group.by="src")+
  coord_fixed()+
  ggtitle("Reclustering of tumor_2 Cells")
#Plot results
DimPlot(tumor_2, pt.size=.1, reduction="umap", label=T)+
  coord_fixed()+
  ggtitle("Reclustering of tumor_2 Cells")
DimPlot(tumor_2, pt.size=.1, reduction="umap", label=T,group.by="seurat_clusters")+
  coord_fixed()+
  ggtitle("Reclustering of tumor_2 Cells")
#kill_cc adds metadata but interferes with analyses
tumor_2<-kill_cc(tumor_2)
Idents(tumor_2)<-tumor$seurat_clusters
for(i in 1:length(mod_names)) {
  tumor_2 <- AddModuleScore(tumor_2, ms[i], name = names(ms[i]))
  p <- FeaturePlot(tumor_2, features = str_c(names(ms[i]), "1"), pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  print(p)
}
saveRDS(tumor_2,'PrimaryTumor/SavedFiles/tumor_2.rds')
tumor_2<-readRDS('PrimaryTumor/SavedFiles/tumor_2.rds')

#Now that sure of tumor cells. Take the tumor subsets from "tumor" and combine across data with harmony
tumor_comb<-subset(tumor, idents=c(1,2,3,4,5,6,9,10,11,12,13,14))
tumor_comb <- tumor_comb %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F)

tumor_comb <- RunHarmony(tumor_comb, group.by.vars = "src")
tumor_comb <- RunUMAP(tumor_comb, reduction = "harmony", dims = 1:30)
tumor_comb <- tumor_comb %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution=0.3)

# Save a stopping point - merged and harmony aligned
qs::qsave(tumor_comb, "PrimaryTumor/SavedFiles/tumor_comb.qs")
# Start from this stopping point
tumor_comb <- qs::qread("PrimaryTumor/SavedFiles/tumor_comb.qs")

DimPlot(tumor_comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("Tumor Clusters post-integration") +
  theme(legend.position = "none")

DimPlot(tumor_comb, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("Tumor by sample post-integration") 

tumor_comb<-kill_cc(tumor_comb)

#Run cell type module as before
Idents(tumor_comb)<-tumor_comb$seurat_clusters
for(i in 1:length(mod_names)) {
  tumor_comb <- AddModuleScore(tumor_comb, ms[i], name = names(ms[i]))
  p <- FeaturePlot(tumor_comb, features = str_c(names(ms[i]), "1"), pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  print(p)
}

#
DimPlot(tumor_comb, reduction = "umap", group.by = "path") +
  coord_fixed() +
  ggtitle("Tumor by pathology type") 

DimPlot(tumor_comb, reduction = "umap", group.by = "type") +
  coord_fixed() +
  ggtitle("Primary or Secondary Tumor") 

#for next step negbiomial or DESeq2
#Identify markers for each cluster identified in the tumor integrated sample tumor_comb. Use negative bioimial
#Idents(tumor_comb) <- tumor_comb$seurat_clusters
#for(i in 1:length(unique(tumor_comb$seurat_clusters))) {
#  p=i-1
#  marks <- FindMarkers(tumor_comb, test.type= "negbinom", ident.1=p)
#  write.table(marks, file=str_c('tumor_combmarks_negbio',p,'.tsv'), quote=FALSE, sep='\t')
#}


#Identify markers for each cluster identified in the tumor integrated sample tumor_comb. Use neg
library(DESeq2)

#DESeq2 analysis. Make iterative loop to generate tsv file with gene list
if(!dir.exists("PrimaryTumor/TumorCombMarks")) {
  dir.create("PrimaryTumor/TumorCombMarks")
}
for(i in 1:length(unique(tumor_comb$seurat_clusters))) {
  p=i-1
  marks <- FindMarkers(tumor_comb, test.type= "DESeq2", ident.1=p)
  write.table(marks, file=str_c('PrimaryTumor/TumorCombMarks/tumor_combmarks',p,'.tsv'), quote=FALSE, sep='\t')
}

#Once I have markers need to begin GSEA Analysis
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
length(unique(tumor_comb$seurat_clusters))

if(!dir.exists("PrimaryTumor/TumorCombGSEGO")) {
  dir.create("PrimaryTumor/TumorCombGSEGO")
}
#Make iterative loop to perform GSEA analysis with GO gene lists for Cellular Component, Biological Processes or Molecular Function
#Output is a ridgeplot and walk plot for the GSEA analysis. These give directional to my results.
#Issue with building modules is that I do not know how to account for negatively expressed genes for markers
#I believe that my analysis relies on genes that are up regulated, as a result I am looking for pathways 
#indicated by GSEA ridgeplot and walkplots which are positively enriched
for(i in 1:length(unique(tumor_comb$seurat_clusters))) {
  #First step is to convert my gene lists on the tsv file into ENTREZIDs. Convert my file to tibble, add column with my gene names from rownames
  p=i-1
  gsea<-read.table(str_c("PrimaryTumor/TumorCombMarks/tumor_combmarks",p,".tsv"), header=T)
  gsea.tib<-as.tibble(gsea)
  gsea_names<- gsea.tib%>%add_column(SYMBOL=rownames(gsea))
  head(gsea_names)
  
  #Make a table for human genes on my list with ENTREZIDs matched in adjacent column (See Go classfication 5.2 on cluster Profiler handbook)
  gene.df<-bitr(gsea_names$SYMBOL, fromType = "SYMBOL",
              toType= c("ENSEMBL", "ENTREZID"),
              OrgDb =org.Hs.eg.db)
  head(gene.df)
  gene.tib<-as.tibble(gene.df)
  #InnerJoin the gene table (gene.tib) and the list of DESEQ 2 fold change data (gsea_names)
  gsea_list<-inner_join(gsea_names, gene.tib,by=c("SYMBOL"))
  head(gsea_list)

  #make vector with Log 2Fold Change and have ENTREZID as names, sort decreasing
  gsea_vector<-as.vector(gsea_list$avg_log2FC)
  names(gsea_vector)<-as.character(gsea_list$ENTREZID)
  head(gsea_vector)
  gsea_vector<-sort(gsea_vector, decreasing=TRUE)
  head(gsea_vector)
  length(gsea_vector)
  x <- unique(names(gsea_vector))
  gsea_vector_2 <- gsea_vector[x]
  gse_vector_2<-sort(gsea_vector_2, decreasing = TRUE)
  length(gsea_vector_2)
  head(gsea_vector_2)

  
  #GO enrichment and overrepresntation 
  ggo<-groupGO( gene = names(gsea_vector_2),
              OrgDb =org.Hs.eg.db,
              keyType = "ENTREZID",
              ont= "CC",
              level = 2,
              readable=TRUE)
  head(ggo)
  ego<- enrichGO(gene = names(gsea_vector_2),
               keyType = "ENTREZID",
               OrgDb = org.Hs.eg.db,
               ont = "MF",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.05,
               readable=TRUE)
  head(ego, n=5)
  #GSEA analysis with GO. One for each subontology
  ego3_bp <- gseGO(geneList     = gsea_vector_2,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  ego3_mf <- gseGO(geneList     = gsea_vector_2,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "MF",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  ego3_cc <- gseGO(geneList     = gsea_vector_2,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "CC",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  
  #plot the ridge plot and walk plot
  p1<-ridgeplot(ego3_cc)+
    ggtitle(str_c("Cluster ",p," CC"))
  p2 <- gseaplot(ego3_cc, geneSetID = 1, by = "runningScore", title = str_c("Cluster ",p," CC::",ego3_cc$Description[1]))
  p3<-ridgeplot(ego3_bp)+
    ggtitle(str_c("Cluster ",p," BP"))
  p4 <- gseaplot(ego3_bp, geneSetID = 1, by = "runningScore", title = str_c("Cluster ",p," BP::",ego3_bp$Description[1]))
  p5<-ridgeplot(ego3_mf)+
    ggtitle(str_c("Cluster ",p," MF"))
  p6 <- gseaplot(ego3_mf, geneSetID = 1, by = "runningScore", title = str_c("Cluster ",p," MF::",ego3_mf$Description[1]))
  
  graph<-cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2, nrow=3, labels=LETTERS[1:3])
  print(graph)
  
  #combine subonotology analysis into one table for each cluster. save as tsv file. This preserves gene list to exam within each GO pathway
  head(ego3_cc)
  CC<- as.tibble(ego3_cc$Description[1:10])
  CC_GOID<-CC%>%add_column(CC_GOID=ego3_cc$ID[1:10])
  CC_core_enrichment<-CC_GOID%>%add_column(CC_core_enrichment=ego3_cc$core_enrichment[1:10])
  BP<-CC_core_enrichment%>%add_column(BP=ego3_bp$Description[1:10])
  BP_GOID<-BP%>%add_column(BP_GOID=ego3_bp$ID[1:10])
  BP_core_enrichment<-BP_GOID%>%add_column(BP_core_enrichment=ego3_bp$core_enrichment[1:10])
  MF<-BP_core_enrichment%>%add_column(MF=ego3_mf$Description[1:10])
  MF_GOID<-MF%>%add_column(MF_GOID=ego3_mf$ID[1:10])
  gseGO<-MF_GOID%>%add_column(MF_core_enrichment=ego3_mf$core_enrichment[1:10])
  head(gseGO,n=10)
  #write tsv file 
  write.table(gseGO, file=str_c('PrimaryTumor/TumorCombGSEGO/gseGO_cluster_',p,'.tsv'), row.names=FALSE, col.name=TRUE,sep="\t")
 
}

#08_02_2021. Identifying Gene Ontologies which are in cluster 0. 
#Making lists in xcel from Top gene converted gene list with the SYMBOLS as my input
#Two types of lists, one with all genes in an ontology or one with just "hits" within the tumor data set
celladhesion<-read.csv("PrimaryTumor/GOPathways/Cluster0_celladhesion_GO.csv")
celladhesion_genehits<-read.csv("PrimaryTumor/GOPathways/Cluster0_celladhesion_genehitsonly_GO.csv")
externalencaps<-read.csv("PrimaryTumor/GOPathways/Cluster0_externalencapsulatingstructure_GO.csv")
externalencaps_genehits<-read.csv("PrimaryTumor/GOPathways/Cluster0_externalencapsulatingstructure_hitsonly_GO.csv")
StrucMolAct_genehits<-read.csv("PrimaryTumor/GOPathways/Cluster0_StructMolAct_GO.csv")

# Make a list of modules want to visualize with FeaturePlot. 
#Pick columns with the SYMBOL to compare to expression of all tumor cells
#Goal is to find a pathway or set of related processes conserved within sub populations of tumor
alexmod <- list()
alexmod$CellAdhesion<- c(celladhesion$Symbol)
alexmod$CellAdhesion_hits<-c(celladhesion_genehits$Gene.Name)
alexmod$ExternalEncap<-c(externalencaps$Symbol)
alexmod$ExternalEncapsulating_hits<-c(externalencaps_genehits$Gene.Symbol)
alexmod$StrucMolAct_genehits<-c(StrucMolAct_genehits$Gene.Name)

module_names <- names(alexmod)

#Plot each module with Feature Plot like before with cell typing modules
for(i in 1:length(module_names)) {
  tumor_comb <- AddModuleScore(tumor_comb, alexmod[i], name = names(alexmod[i]))
  p <- FeaturePlot(tumor_comb, features = str_c(names(alexmod[i]), "1"), pt.size = 1,
                   order = T, cols = c("lightgoldenrod", "darkred"))
  print(p)
}



#Repeat for other clusters where all samples seem to associate. GOal is to find subcluster modules which have high score across osteo types.
#Might ignore clusters 3 and 5 because they seem to be dividing cells. 
#Also compare to osteoblasts in culture? If same module subclustering could be developmental