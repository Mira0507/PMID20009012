library(tidyverse)
library(ggplot2)
library(DESeq2)
library(pheatmap)

# raw count data: it has to be a matrix of integers 
ct <- read_tsv("gilad_count_table.txt")
ct1 <- as.matrix(ct[, 2:ncol(ct)])
rownames(ct1) <- ct$gene
ct2 <- apply(ct1, 2, as.integer)

# metadata 
md <- read_delim("gilad_phenodata.txt", delim = " ") %>%
        mutate(gender = factor(gender))
md1 <- as.matrix(md[, "gender"])
rownames(md1) <- md$sample.id

# creating DESeq2 object des 
des <- DESeqDataSetFromMatrix(countData = ct2,
                              colData = md1, 
                              design = ~ gender)

# normalizing counts
des <- estimateSizeFactors(des)

# extracting normalized counts 
des <- counts(des, normalized = TRUE) 
des <- apply(des, 2, as.integer)


# logarithmic transformation for quality control
vsd <- vst(des, blind = TRUE)

# hierarchical clustering 
log_mat_des <- assays(vsd)
log_cor_des <- cor(log_mat_des)

# Quality control: PCA 
plotPCA(des, intgroup = "gender")





