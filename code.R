library(tidyverse)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


############################### data cleaning ################################

rct <- read_csv("rawcounts.csv")
rct1 <- rct %>% 
        select(GeneSymbol, R1L1.HSM1:R8L2.HSF3)

samples <- colnames(select(rct1, -GeneSymbol))

rct2 <- rct1 %>%
        gather(sample_id, counts, - GeneSymbol) %>%
        mutate(sp = case_when(
                str_detect(sample_id, "HSM") ~ "HS_Male",
                str_detect(sample_id, "HSF") ~ "HS_Female",
                str_detect(sample_id, "RMM") ~ "RM_Male",
                str_detect(sample_id, "RMF") ~ "RM_Female",
                str_detect(sample_id, "PTM") ~ "PT_Male",
                str_detect(sample_id, "PTF") ~ "PT_Female")) %>%
        separate(sp, c("species", "gender"), sep = "_") %>%
        select(sample_id, species, gender)

mdata <- data.frame(sample = colnames(rct1 %>% select(- GeneSymbol))) %>%
        inner_join(rct2, by = c("sample" = "sample_id"))

mdata1 <- unique(mdata)

# metadata: mdata2 
mdata2 <- as.matrix(select(mdata1, -sample))
rownames(mdata2) <- mdata1$sample

# raw count data: rct3
rct3 <- as.matrix(select(rct1, - GeneSymbol))
rownames(rct3) <- rct1$GeneSymbol




#################################### DESEQ #######################################


# creating DESeq2 object des 
des <- DESeqDataSetFromMatrix(countData = rct3,
                              colData = mdata2, 
                              design = ~ species + gender)

# normalization
des1 <- estimateSizeFactors(des)
rct_norm <- counts(des1, normalized = TRUE) %>%
        apply(2, as.integer)

des_norm <- DESeqDataSetFromMatrix(countData = rct_norm,
                                   colData = mdata2, 
                                   design = ~ species + gender)

vsd <- vst(des_norm, blind = TRUE) 


# Run DEseq 
des_deseq <- DESeq(des_norm)
res_names <- resultsNames(des_deseq)

# extracting results
res1 <- results(des_deseq,
                name = "gender_Male_vs_Female",
                alpha = 0.05)

res2 <- results(des_deseq, 
                name = "species_PT_vs_HS",
                alpha = 0.05)

res3 <- results(des_deseq,
                name = "species_RM_vs_HS",
                alpha = 0.05)

# shrinkage
res1_s <- lfcShrink(des_deseq,
                  coef = "gender_Male_vs_Female",
                  type = "apeglm") 

res2_s <- lfcShrink(des_deseq,
                    coef = "species_PT_vs_HS",
                    type = "apeglm")  
res3_s <- lfcShrink(des_deseq,
                   coef = "species_RM_vs_HS",
                   type = "apeglm")



############################## final data cleaning ###############################

# data cleaning for plotting 
m_vs_f <- cbind(gene = rct1$GeneSymbol, as.data.frame(res1_s))





# modeling dispersion
mean_counts <- rowMeans(rct_norm)
var_counts <- rowVars(rct_norm)
disp_df <- data.frame(mean_counts, var_counts)
disp_plot <- ggplot(disp_df, 
                    aes(x = mean_counts, 
                        y = var_counts)) +
        geom_point(alpha = 0.1) + 
        scale_x_log10() +
        scale_y_log10() + 
        theme_bw()





#################################### plotting #################################


# correlation heatmap 
htmap <- vsd %>%
        assay() %>%
        cor() %>%
        pheatmap(annotation = select(as.data.frame(mdata2), c(species, gender)),
                 main = "Correlation Heatmap")

# PCA
pca <- plotPCA(vsd, intgroup = c("gender", "species")) +
        theme_bw() +
        ggtitle("PCA")

