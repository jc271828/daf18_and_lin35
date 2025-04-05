
windowsFonts(sans = windowsFont("Arial"))

# load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(eulerr)
library(gtools)
library(arrangements)
library(dendextend)
library(gplots)
library(reshape2)
library(colorspace)
library(RColorBrewer)

# for arranging panels
library(cowplot)
library(egg)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)


# define a function -------------------------------------------------------

# define a new function called "hypergeo"
hypergeo <- function(target_genes, DE_genes, all_genes) {
  # white balls in urn
  m <- intersect(target_genes, all_genes) %>% length()
  # black balls in urn
  n <- length(all_genes) - m 
  # balls drawn
  k <- length(DE_genes)
  # white in drawn balls
  x <- intersect(target_genes, DE_genes) %>% length()
  phyper(x-1, m, n, k, lower.tail = FALSE) %>% return()
}


# heat map -----------------------------------------------------------------

# create color palettes
morecols <- colorRampPalette(c("cornflowerblue",
                               "black",
                               "khaki1"))

# load data
anova_DEGs_woFed <-
  read_xlsx("S1 Data.xlsx",
            sheet = "glmQLFTest",
            col_names = TRUE) %>% 
  filter(FDR < 0.05)

cpm <- 
  read_xlsx("S1 Data.xlsx",
            sheet = "CPM",
            col_names = TRUE)

# clean up data
cpm_woFed <- filter(cpm, gene_id %in% anova_DEGs_woFed$gene_id)
cpm_woFed <- cpm_woFed[, c(12:20, 25:27, 1)]

logCPM_woFed <- log2(cpm_woFed[, 1:12] + 1)
AvgLogCPM_woFed <-
  cbind(as.data.frame(rowMeans(logCPM_woFed[, c(1:3)])),
        as.data.frame(rowMeans(logCPM_woFed[, c(4:6)])),
        as.data.frame(rowMeans(logCPM_woFed[, c(7:9)])),
        as.data.frame(rowMeans(logCPM_woFed[, c(10:12)])))
colnames(AvgLogCPM_woFed) <- c("daf-16",
                               "daf-18",
                               "daf-16; daf-18",
                               "wildtype")
rownames(AvgLogCPM_woFed) <- cpm_woFed$gene_id

# gene-level pearson correlation
gene_cor_woFed <- cor(t(AvgLogCPM_woFed), method = "pearson", use = "pairwise.complete.obs")
gene_distance_woFed <- as.dist((1-gene_cor_woFed)/2)
gene_dend_woFed <- 
  hclust(gene_distance_woFed, method = "complete") %>%
  as.dendrogram()
nleaves(gene_dend_woFed)
nnodes(gene_dend_woFed)
plot(gene_dend_woFed)
gene_dend_woFed <- color_branches(gene_dend_woFed, k = 5, col = rainbow_hcl)
gene_dend_woFed <- color_labels(gene_dend_woFed, k = 5, col = rainbow_hcl)
plot(gene_dend_woFed)

Znorm_woFed <- t(scale(t(AvgLogCPM_woFed), center = TRUE, scale = TRUE))

# for genes, before and after z-score transformation, correlation doesn't change, but it DOES change for conditions
# condition-level pearson correlation
condition_cor_woFed <- cor(Znorm_woFed, method = "pearson", use = "pairwise.complete.obs")
condition_distance_woFed <- as.dist((1-condition_cor_woFed)/2)
condition_dend_woFed <- 
  hclust(condition_distance_woFed, method = "complete") %>%
  as.dendrogram()
nleaves(condition_dend_woFed)
nnodes(condition_dend_woFed)
plot(condition_dend_woFed)

tmp_woFed <- 
  as.data.frame(Znorm_woFed) %>%
  apply(c(1,2), as.numeric)

# reorder gene dendrogram based on wildtype mean to make the heatmap look nice (doesn't change dendrogram leaf relationships)
# gene_dend_woFed <- reorder(gene_dend_woFed, as.data.frame(tmp_woFed)$wildtype, agglo.FUN = mean)
col_labels_woFed <- get_leaves_branches_col(gene_dend_woFed)
col_labels_woFed <- col_labels_woFed[order(order.dendrogram(gene_dend_woFed))]

# plot heat map
p1 <-
  heatmap.2(tmp_woFed,
            main = "1416 DEGs from ANOVA-like glmLRT\nsorted in descending order of starved wildtype",
            dendrogram = "both",
            Rowv = gene_dend_woFed,
            Colv = condition_dend_woFed,
            col = morecols(50),
            trace = "none",
            density.info = "none",
            margins = c(5,5),
            labRow = FALSE,
            cexCol = 1.5,
            xlab = "condition",
            srtCol = 0,
            adjCol = c(0.5,0.5),
            key.par=list(mar = c(5,1,3,1)),
            key.xlab = "z-score",
            key.title = NA,
            RowSideColors = col_labels_woFed)


# d18 x d16 DEG Venn ------------------------------------------------------

fdr_cutoff <- 0.05

# load data
d18 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-18_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% .$gene_id
d16 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% .$gene_id
d18_minus_d16 <- setdiff(d18, d16)

set.seed(1)
sets <- list(d18 = d18, d16 = d16)
p3 <- plot(euler(sets, shape = "ellipse"), quantities = TRUE)
print(p3)
hypergeo(d16, d18, WBid_background) # 0


# d18 - d16 x Cui 2013 lin-35 targets Venn -------------------------------------

background <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16;daf-18_vs_daf-16")
WBid_background <- background$gene_id

GeneLists <- read_csv("DREAM_lin35_Becky_Tepper_Amanda_gene_lists.csv", col_names = TRUE)
cui_lin35_targets <- GeneLists %>% filter(Term %in% c("down by lin-35(n745), starved L1", "UP by lin-35(n745), starved L1"))

set.seed(1)
sets <- list(d18_minus_d16 = d18_minus_d16, cui_lin35_targets = intersect(cui_lin35_targets$wb_id, WBid_background))
p4 <- plot(euler(sets, shape = "ellipse"), quantities = TRUE)
print(p4)
hypergeo(cui_lin35_targets$wb_id, d18_minus_d16, WBid_background) # 1.088394e-15


# d18 - d16 x Cui 2013 lin-35 targets CDF -------------------------------------

temp <- background %>% filter(gene_id %in% intersect(cui_lin35_targets$wb_id, WBid_background))
colnames(temp)[1] <- "wb_id"
temp %<>% merge(cui_lin35_targets, by = "wb_id", all = F)

temp$Term %>% table()

p5 <-
  ggplot() +
  stat_ecdf(data = background, aes(logFC), linewidth = 1) +
  stat_ecdf(data = temp, aes(logFC, color = Term), linewidth = 1) +
  labs(tag = "E", x = "log2FC (daf-16; daf-18/daf-16)", y = "Cumulative proportion") +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = expansion(mult = c(0, 0)))

down_ks_test_p <-
  ks.test(background$logFC,
          filter(temp, Term == "down by lin-35(n745), starved L1") %>% .$logFC,
          alternative = "two.sided")$p.value %>%
  signif(2)
  
up_ks_test_p <- 
  ks.test(background$logFC,
          filter(temp, Term == "UP by lin-35(n745), starved L1") %>% .$logFC,
          alternative = "two.sided")$p.value %>%
  signif(2)

p5 <-
  p5 +
  ggtitle(paste0("Down in lin-35 Mut K-S test p ", down_ks_test_p, "\nUp in lin-35 Mut K-S test p ", up_ks_test_p, collapse = ""))

print(p5)


# d18 - d16 x Gal 2021 direct lin-35 targets Venn -------------------------------------

gal_lin35_pos_targets <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin35_downregulated.txt", col_names = F) %>%
  .$X1
gal_lin35_neg_targets <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin35_upregulated.txt", col_names = F) %>%
  .$X1
gal_lin35_direct_pos_targets <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin35_direct_downregulated.txt", col_names = F) %>%
  .$X1
gal_lin35_direct_neg_targets <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin35_direct_upregulated.txt", col_names = F) %>%
  .$X1

gal_lin35_targets <- union(gal_lin35_neg_targets, gal_lin35_pos_targets)
gal_lin35_direct_targets <- union(gal_lin35_direct_neg_targets, gal_lin35_direct_pos_targets)

set.seed(1)
sets <- list(d18_minus_d16 = d18_minus_d16, gal_lin35_direct_targets = intersect(gal_lin35_direct_targets, WBid_background))
p6 <- plot(euler(sets, shape = "ellipse"), quantities = TRUE)
print(p6)
hypergeo(gal_lin35_direct_targets, d18_minus_d16, WBid_background) # 8.01025e-36


# d18 - d16 x Gal 2021 direct lin-35 targets CDF -------------------------------------

gal_lin35_direct_targets <-
  data.frame(gene_id = na.omit(gal_lin35_direct_neg_targets), Term = "UP_in_lin-35_direct_targets_Gal2021") %>%
  rbind(data.frame(gene_id = na.omit(gal_lin35_direct_pos_targets), Term = "DOWN_in_lin-35_direct_targets_Gal2021"))
temp <- background %>% filter(gene_id %in% intersect(gal_lin35_direct_targets$gene_id, WBid_background))
temp %<>% merge(gal_lin35_direct_targets, by = "gene_id", all = F)

temp$Term %>% table()

p7 <-
  ggplot() +
  stat_ecdf(data = background, aes(logFC), linewidth = 2) +
  stat_ecdf(data = temp, aes(logFC, color = Term), linewidth = 2) +
  labs(tag = "G", x = "log2FC (double/daf-16)", y = "Cumulative proportion") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits = c(-1, 2))

down_ks_test_p <-
  ks.test(background$logFC,
          filter(temp, Term == "DOWN_in_lin-35_direct_targets_Gal2021") %>% .$logFC,
          alternative = "two.sided")$p.value %>%
  signif(2)

up_ks_test_p <-
  ks.test(background$logFC,
          filter(temp, Term == "UP_in_lin-35_direct_targets_Gal2021") %>% .$logFC,
          alternative = "two.sided")$p.value %>%
  signif(2)

p7 <-
  p7 +
  ggtitle(paste0("Down in lin-35 Mut K-S test p ", down_ks_test_p, "\nUp in lin-35 Mut K-S test p ", up_ks_test_p, collapse = ""))

print(p7)


# arrange panels 3-7 ------------------------------------------------------

cairo_pdf("fig2_panels3-7.pdf", width = 25, height = 25)
plot_grid(NULL, p3, p4, p5, p6, p7, ncol = 2)
dev.off()


