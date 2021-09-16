library(tibble)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(stringr)

#ComplexHeatmap
cx.pd <- read.delim("C:/Users/rssxr002/Downloads/cxpd.tsv", row.names = 1)

#remove "HALLMARK_"
c <- rownames(cx.pd)
c <- gsub("HALLMARK_", "", c)
#remove "_" by removing special characters
c <- gsub("_", " ", c)
rownames(cx.pd) <- c
cx.pd.log <- -log10(cx.pd) #log transform
  
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

#Include percentage of cells in lung for each cluster
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

##############################Figure 3 ################################
P1 <- df %>%
        as.data.frame() %>%
        rownames_to_column(var = "Sample") %>%
        full_join(aka4 %>% rownames_to_column(var = "Sample")) %>%
        pivot_longer(c(-Sample, -Lpercent),
                     names_to = "Pathway",
                     values_to = "pval",
                     names_repair = "minimal") %>%
        mutate(Signif = pval >= (-1 * log10(0.05)),
               sample_order = str_remove(Sample, "[0-9]") %>%
                       rank() * 1000 - Lpercent,
               Sample = reorder(Sample, sample_order)) %>%
        ggplot(., aes(x = Sample, y = Pathway, fill = Signif)) + 
        geom_tile() +
        geom_text(aes(label = sprintf("%0.2f", pval))) +
        scale_fill_manual(values = c("gray", "red")) +
        scale_y_discrete(limits = rev) +
        theme(legend.position = "none") +
        ylab("") 
P2 <- aka4 %>%
        rownames_to_column(var = "Sample") %>%
        full_join(aka2 %>%
                          rownames_to_column(var = "Sample")) %>%
        mutate(sample_order = str_remove(Sample, "[0-9]") %>%
                       rank() * 1000 - Lpercent,
               Sample = reorder(Sample, sample_order),
               ID = reorder(ID, sample_order)) %>%
        ggplot(., aes(x = Sample, y = Lpercent, fill = ID)) +
        geom_bar(stat = "identity") +
        ylab("Percent\nin lung") +
        xlab("") +
        theme(legend.position = "top", ) +
        labs(fill = "")
cowplot::plot_grid(P2, P1, ncol = 1, align = "v",rel_heights = c(2, 10))


##############################Figure 1 ################################
Fig1 <- read_delim("C:/Users/rssxr002//Downloads/Fig1.txt", 
                   delim = "\t", 
                   col_names = TRUE)
colnames(Fig1) <- c("Pathway", "A0", "A1", "A2", "A3", "B0", "B1", "B2", "B3")
Fig1 %>%
        column_to_rownames(var = "Pathway") %>%
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
