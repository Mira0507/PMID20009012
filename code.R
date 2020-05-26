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
rct3 <- apply(rct3, 2, as.integer)




#################################### DESEQ #######################################


# creating DESeq2 object des 
des <- DESeqDataSetFromMatrix(countData = rct3,
                              colData = mdata2, 
                              design = ~ species + gender)

# normalization
des1 <- estimateSizeFactors(des)
rct_norm <- counts(des1, normalized = TRUE) 
rownames(rct_norm) <- rct1$GeneSymbol
rct_norm1 <- rct_norm
rct_norm <- apply(rct_norm, 2, as.integer)

des_norm <- DESeqDataSetFromMatrix(countData = rct_norm,
                                   colData = mdata2, 
                                   design = ~ species + gender)

vsd <- vst(des_norm, blind = TRUE) 

# modeling dispersion
mean_counts <- rowMeans(rct_norm)
var_counts <- rowVars(rct_norm)
disp_df <- data.frame(mean_counts, var_counts)



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
sig_table <- function(res, cat) {
        df <- cbind(gene = rct1$GeneSymbol,
                    as.data.frame(res)) %>%
                mutate(Comparison = cat) %>%
                filter(padj < 0.05) %>% 
                arrange(desc(log2FoldChange))
        return(df)
}

# final result (log2 fold change) with significantly changed genes
m_vs_f <- sig_table(res1_s, "Male over Female")
p_vs_h <- sig_table(res2_s, "Chimpanzee over Human")
r_vs_h <- sig_table(res3_s, "Rhesus Macaque over Human")


comparison_table <- rbind(m_vs_f,
                          p_vs_h,
                          r_vs_h) 

# normalized count data table for plotting
rct_norm_df <- as.data.frame(rct_norm1) %>%
        rownames_to_column(var = "gene") %>%
        gather(sample, normalized_counts, -gene) %>%
        filter(normalized_counts > 0) %>%
        mutate(log10_normalized_counts = log10(normalized_counts)) %>%
        inner_join(mdata1, by = "sample") %>%
        mutate(sample = factor(sample, 
                               levels = c("R1L1.HSM1",
                                          "R2L3.HSM2",
                                          "R5L2.HSM1",
                                          "R3L6.HSM3",
                                          "R4L1.HSM3",
                                          "R4L8.HSM2",
                                          
                                          "R1L4.HSF1",
                                          "R8L2.HSF3",
                                          "R8L1.HSF3",
                                          "R4L2.HSF1",
                                          "R3L2.HSF2",
                                          "R2L7.HSF2",
                                          
                                          "R1L6.PTM1",
                                          "R2L8.PTM2",
                                          "R3L3.PTM1",
                                          "R6L4.PTM3",
                                          "R6L2.PTM3",
                                          "R4L6.PTM2",
                                          
                                          "R2L4.PTF2",
                                          "R3L7.PTF3",
                                          "R6L6.PTF2",
                                          "R5L3.PTF3",
                                          "R4L4.PTF1",
                                          "R1L2.PTF1",
                                          
                                          "R1L3.RMM1",
                                          "R2L6.RMM2",
                                          "R3L1.RMM3",
                                          "R3L8.RMM1",
                                          "R4L3.RMM3",
                                          "R5L4.RMM2",
                                          
                                          "R1L7.RMF1",
                                          "R2L2.RMF2",
                                          "R3L4.RMF3",
                                          "R4L7.RMF3",
                                          "R5L1.RMF1",
                                          "R5L8.RMF2")))

            











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

# dispersion plot
disp_plot <- ggplot(disp_df, 
                    aes(x = mean_counts, 
                        y = var_counts)) +
        geom_point(alpha = 0.1) + 
        scale_x_log10() +
        scale_y_log10() + 
        theme_bw() + 
        ggtitle("Dispersion Plot") + 
        xlab("Counts Mean") + 
        ylab("Counts Variance")

# MA plots

MA_m_vs_f <- plotMA(res1, 
                    ylim = c(-2, 2),
                    main = "MA Plot of Contrast between Genders without Shrinkage")
MA_m_vs_f_shr <- plotMA(res1_s, 
                        ylim = c(-2, 2),
                        main = "MA Plot of Contrast between Genders with Shrinkage")

MA_p_vs_h <- plotMA(res2, 
                    ylim = c(-8, 8),
                    main = "MA Plot of Contrast between Human and Chimpanzee without Shrinkage")
MA_p_vs_h_shr <- plotMA(res2_s, 
                        ylim = c(-8, 8),
                        main = "MA Plot of Contrast between Human and Chimpanzee with Shrinkage")

MA_r_vs_h <- plotMA(res3,
                    ylim = c(-9, 9),
                    main = "MA Plot of Contrast between Human and Rhesus Macaque without Shrinkage")
MA_r_vs_h_shr <- plotMA(res3_s,
                        ylim = c(-9, 9),
                        main = "MA Plot of Contrast between Human and Rhesus Macaque with Shrinkage")


# gene expression heatmap 

res_heat_map <-
        ggplot(rct_norm_df, 
               aes(x = sample, 
                   y = gene, 
                   fill = log10_normalized_counts,
                   color = log10_normalized_counts)) +
        geom_tile() +
        theme(axis.text.y = element_blank(),
              axis.text.x = element_text(size = 7, angle = 90)) +
        ggtitle("Gene Expression Profile") +
        xlab("Sample") + 
        ylab("Gene")
                           
