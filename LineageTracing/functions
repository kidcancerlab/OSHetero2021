# Set up environment, activate library components
library(ggsci)
library("cowplot")
library("dplyr")
library("Matrix")
library("reticulate")
library("Seurat")
library("reshape2")
library("ggplot2")
# library("harmony")
# library("future")

# plan(strategy = "multicore", workers = 3)


## Define functions# For processLTBC (lineage tracing barcodes), the following values should be set:
#    sobject - the name of the Seurat object that contains the cell barcodes (cids) that will be matched and integrated
#    lt.loc - the file location of the lineage tracing barcodes (extracted from the fastq files)
#    cid.loc - the file location of the cell id barcodes (extracted from the fastq files)
#    histogram - TRUE will trigger the function to generate and output a histogram plot of the top 40 most frequent lineage tracing barcodes
#    col.fill - the color that you would like to use for the barchart on the histogram
#    ymax - if set, this will define the upper limit of the y axis (ie, for creating side-by-side comparisons)
#    relative - TRUE will normalize cell counts to total number of cells containing barcodes
#    title - verbiage for the title of the histogram, if triggered
#    ret.list - TRUE will trigger the function to return the list of barcode frequencies, rather than the seurat object
# The following procedure will generate the matching lineage tracing and cell ID tables:
#   In Linux bash, navigate to the scRNA Counts directory containing the cellranger output .bam
#   Load samtools with > module load samtools-1.7
#   Extract all lines containing a read that matches the flanking sequences for the LT barcode
#     > samtools view possorted_genome_bam.bam | egrep "CGCA[ACGT]{14}TGGT[ACGT]{30}TGGT" > match.lt.sam
#     (You might want to make sure this is actually running when you submit the command, as
#       I've been burned before.  It takes forever even when it's working on a big file.  Open
#       a separate Putty window, SSH into the node that should be running samtools and egrep, 
#       then run top.  You should see one process for each of these programs running.)
#   Remove any lines that do not have a high-confidence cell ID assigned:
#     > egrep "CB:Z:[ACGT]{16}" match.lt.sam > match.both.sam
#   Extract the lineage tracing barcode and flanking sequences
#     > egrep -o "CBCA[ACGT]{14}TGGT[ACGT]{30}TGGT" match.both.sam > lt.fq
#   Extract the matching cell identity barcodes
#     > egrep -o "CB:Z:[ACGT]{16}" match.both.sam > pre-cid.fq
#   Clean up the cell identity barcodes (get rid of the prefix)
#     > egrep -o "[ACGT]{16}$" pre-cid.fq > cid.fq
#or 	> cut -d":" -f3 pre-cid.fq > cid.fq
#   Input lt.fq and cid.fq into the lt.loc and cid.loc variables in the function below

processLTBC <- function(sobject, lt.loc, cid.loc, histogram = F, col.fill = "#4DBBD5FF", ymax = NA, relative = F, title = "LT Barcode Frequency", ret.list = F) {
  
  # Read in the master Cellecta barcode tables for QC purposes
  bc14 <- read.table("R:/RESRoberts/Sanjana/Lineage Tracing/Cellecta/Cellecta-bc14s.txt")
  bc30 <- read.table("R:/RESRoberts/Sanjana/Lineage Tracing/Cellecta/Cellecta-bc30s.txt")
  names(bc14) <- c("label", "forward", "reverse")
  names(bc30) <- c("label", "forward", "reverse")
  bc14f <- setNames(bc14$label, bc14$forward)
  bc30f <- setNames(bc30$label, bc30$forward)
  
  # Import the lineage tracing and matching cell ID tables extracted from the fastqs
  lt.bc <- read.table(lt.loc)
  cid.bc <- read.table(cid.loc)
  cid.bc <- paste(cid.bc[[1]],"-1", sep = "")
  bc <- data.frame(cid.bc, lt.bc)
  colnames(bc) <- c("cid", "lt")
  
  # Deduplicate redundant reads
  bc <- dplyr::distinct(bc)
  
  # Extract the bc14 and the bc30 reads from the lineage tracing sequences (identified by flanking regions)
  bc$bc14 <- substring(bc$lt, 5, 18)
  bc$bc30 <- substring(bc$lt, 23, 52)
  
  # Match extracted barcode reads against the Cellecta barcode tables
  bc$label14 <- bc14f[bc$bc14]
  bc$label30 <- bc30f[bc$bc30]
  bc$label14 <- substring(bc$label14, 6)
  bc$label30 <- substring(bc$label30, 6)
  
  # Eliminate barcodes that don't match the Cellecta barcode tables
  bc <- na.omit(bc, cols = c(label14, label30))
  
  # Concatenate the two barcodes into a single compound column
  bc$label <- paste(bc$label14, bc$label30, sep = "-")
  
  # Integrate the lineage tracing barcode into the Seurat object metadata
  bc <- setNames(as.character(bc$label), bc$cid)
  bc <- sort(bc)
  sobject$lt <- bc[sobject@assays$RNA@counts@Dimnames[[2]]]
  
  # Generate the frequency tables
  ylabel = "Number of Cells"
  bc.freq <- as.data.frame(table(sobject$lt))
  bc.freq <- bc.freq[order(-bc.freq$Freq),]
  if(isTRUE(relative)) {
    bc.freq$Freq = bc.freq$Freq/length(bc)*100
    ylabel = "Percentage of Cells"
  }
  bc.plot.data <- head(bc.freq, n=40)
  
  # Create histogram graphs if desired (default using blue color from npg from ggsci)
  if(isTRUE(histogram)) {
    p <- ggplot(bc.plot.data, aes(x = reorder(Var1, -Freq), y = Freq)) +
      geom_bar(fill = col.fill, stat = "identity") +
      ggtitle(title) +
      ylab(ylabel) +
      xlab("Lineage Barcode") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p + ylim(0, ymax))
  } 
  if (isTRUE(ret.list)) {
    return(bc.freq)
  } else {
    return(sobject)
  }
}


