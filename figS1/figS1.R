
windowsFonts(sans = windowsFont("Arial"))

# load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)
library(nlme)
library(outliers)
library(edgeR)
library(arrangements)
library(reshape2)
library(writexl)
library(readxl)
library(eulerr)

# for arranging panels
library(cowplot)
library(egg)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)


# define function "hypergeo" -----------------------------------------------------------------

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


# define color code -------------------------------------------------------

color_wt <- "#000000"
color_d18 <- "red"
color_d16 <- "orange"
color_d16_d18 <- "purple"


# load data ---------------------------------------------------------------

fdr_cutoff <- 0.05

HE_per_hour <- read.csv("hatching_curve.csv", header = TRUE) %>% subset(strain != "n2_new")

raw <- read_csv("C:/Users/jingx/Duke Bio_Ea Dropbox/Baugh Lab/2_baugh_lab_users/Jingxian/rna_seq/2020.09.02_daf-18_L1/Chitrakar_6460_200831B7/analyses_20200902/final_version_of_analysis/step1/20200902_RNAseq_raw_reads_merged.csv", col_names = TRUE) %>% as.data.frame()
ws273_GeneInfo <- read_csv("C:/Users/jingx/Duke Bio_Ea Dropbox/Baugh Lab/2_baugh_lab_users/Jingxian/rna_seq/2020.09.02_daf-18_L1/Chitrakar_6460_200831B7/analyses_20200902/final_version_of_analysis/step1/ws273_gene_info.csv", col_names = TRUE) %>% as.data.frame()

d18 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-18_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% .$gene_id
d16 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% .$gene_id
d16_minus_d18 <- setdiff(d16, d18)

gal_lin35_direct_pos_targets <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin35_direct_downregulated.txt", col_names = F) %>%
  .$X1
gal_lin35_direct_neg_targets <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin35_direct_upregulated.txt", col_names = F) %>%
  .$X1

gal_lin35_direct_targets <- union(gal_lin35_direct_neg_targets, gal_lin35_direct_pos_targets)


# hatching curves ----------------------------------------------------------

# censor outliers
daf18 <- HE_per_hour %>% filter(hours_after_bleach == 16 & strain == "daf18")

# opposite = FALSE means it's testing whether the lowest value is an outlier
dixon.test(daf18$hatching_efficiency, opposite = FALSE) # p-value = 0.08463

# censor rep2
HE_per_hour %<>% filter(!(rep == "rep2" & hours_after_bleach == 16 & strain == "daf18"))

# clean up data
HE_per_hour$strain %<>%
  str_replace("n2", "wild type") %>%
  str_replace("daf16", "daf-16(mu86)") %>%
  str_replace("daf18", "daf-18(ok480)") %>%
  str_replace("double", "daf-16(mu86); daf-18(ok480)")

HE_per_hour$strain %<>% as.factor() %>% fct_relevel(c("wild type", "daf-16(mu86)", "daf-18(ok480)", "daf-16(mu86); daf-18(ok480)"))

stat <-
  HE_per_hour %>%
  group_by(strain, hours_after_bleach) %>%
  summarise(mean(hatching_efficiency), sd(hatching_efficiency))

colnames(stat)[3] <- "mean"
colnames(stat)[4] <- "sd"

stat$strain %<>%
  str_replace("n2", "wild type") %>%
  str_replace("daf16", "daf-16(mu86)") %>%
  str_replace("daf18", "daf-18(ok480)") %>%
  str_replace("double", "daf-16(mu86); daf-18(ok480)")

stat$strain %<>% as.factor() %>% fct_relevel(c("wild type", "daf-16(mu86)", "daf-18(ok480)", "daf-16(mu86); daf-18(ok480)"))

stat$linetype <- NA
stat[str_detect(stat$strain, "ok480", negate = TRUE), "linetype"] <- 1 # 1 is solid
stat[str_detect(stat$strain, "ok480"), "linetype"] <- 2 # 2 is dashed

stat$shape <- NA
stat[str_detect(stat$strain, "ok480", negate = TRUE), "shape"] <- 16 # 16 is closed circle
stat[str_detect(stat$strain, "ok480"), "shape"] <- 21 # 21 is open circle

# three ways of plotting the data
stat$strain %>% levels()

print(paste('Average number of sampled embryos is', round(mean(HE_per_hour$L1s_and_eggs), 0)))
print(paste('Standard deviation of sampled embryos is', round(sd(HE_per_hour$L1s_and_eggs), 0)))

p1.1<-
  ggplot(data = stat, mapping = aes(x = hours_after_bleach)) +
  geom_point(mapping = aes(y = mean, color = strain, shape = shape), size = 2) +
  geom_line(mapping = aes(y = mean, color = strain, linetype = linetype), linewidth = 1) +
  geom_errorbar(mapping = aes(ymin = mean - sd,
                              ymax = mean + sd,
                              width = 0.5,
                              color = strain,
                              linetype = linetype),
                linewidth = 0.5) +
  scale_shape_identity() +
  scale_linetype_identity() +
  scale_color_manual(values = c(color_wt, color_d16, color_d18, color_d16_d18)) +
  # scale_x_continuous(breaks = seq(12, 18, by = 1), limits = c(11.9, 18.1), expand = expansion(mult = c(0, 0)))+
  # scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1), expand = expansion(mult = c(0, 0)))+
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "A",
       x = "Duration of starvation (days)",
       y = "Proportion alive normalized by the first day",
       title = "hatching efficiency over time (4 reps)")

p1.1

p1.2 <- 
  ggplot(data = stat, mapping = aes(x = hours_after_bleach)) +
  geom_point(mapping = aes(y = mean, color = strain), size = 2) +
  geom_line(mapping = aes(y = mean, color = strain), linewidth = 1) +
  geom_errorbar(mapping = aes(ymin = mean - sd,
                              ymax = mean + sd,
                              width = 0.5,
                              color = strain),
                linewidth = 0.5) +
  scale_color_manual(values = c(color_wt, color_d16, color_d18, color_d16_d18)) +
  # scale_x_continuous(breaks = seq(12, 18, by = 1), limits = c(11.9, 18.1), expand = expansion(mult = c(0, 0)))+
  # scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1), expand = expansion(mult = c(0, 0)))+
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "A",
       x = "Hours post bleach",
       y = "Proportion hatched",
       title = "hatching efficiency over time (4 reps)")

p1.2

HE_per_hour$strain %>% levels()

p1.3 <- 
  ggplot(data = HE_per_hour,
         mapping = aes(x = hours_after_bleach,
                       y = hatching_efficiency,
                       color = strain)) +
  stat_summary(geom = "crossbar",
               fun = mean,
               width = 0.3,
               linewidth = 0.5) +
  stat_summary(geom = "line", mapping = aes(group = strain, color = strain), linewidth = 1) +
  geom_point(size = 3, shape = 21) +
  scale_color_manual(values = c(color_wt, color_d16, color_d18, color_d16_d18)) +
  theme_classic() +
  theme(aspect.ratio = 0.5,
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(from = 0, to = 1, by = 0.2),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(tag = "A",
       x = "Hours post bleach",
       y = "Proportion hatched",
       title = "hatching efficiency over time (4 reps)")

p1.3


# clean up RNA-seq analysis data ------------------------------------------

# tidy data
colnames(ws273_GeneInfo)[1] <- "gene_id"
raw_merged <- merge(ws273_GeneInfo, raw, by = "gene_id", sort = FALSE)

# get rid of four "outlier" samples
raw_merged %<>% select(-c(colnames(raw_merged)[c(2:4, 6:13)],
                          "N2_starv_sample12",
                          "daf-16_sample13",
                          "daf-18_sample24",
                          "double_sample25"))

# convert data frame to a longer form
raw_merged %<>%
  pivot_longer(cols = -c("gene_id", "type"),
               names_to = "group",
               values_to = "counts") %>%
  as.data.frame()

# put data in a DGEList object
raw_merged_wide <-
  raw_merged %>%
  filter(type == "protein_coding_gene") %>%
  pivot_wider(names_from = "group",
              values_from = "counts") %>%
  as.data.frame() %>%
  select(-c("type"))

df <- DGEList(counts = raw_merged_wide[, c(2:17)],
              genes = raw_merged_wide[, 1],
              group = substr(colnames(raw_merged_wide), start = 1, stop = nchar(colnames(raw_merged_wide))-9)[2:17]) 

df <- df[order(rowSums(df$counts), decreasing = TRUE),]

rownames(df$counts) <- df$genes$genes

# CPM filtering
keep <- rowSums(cpm(df) > 1) >= 3
df <- df[keep, , keep.lib.sizes = FALSE]

# TMM normalization
df <- calcNormFactors(df, method = "TMM")
cpm <- cpm(df, normalized.lib.sizes = TRUE) %>% as.data.frame()
cpm$gene_id <- rownames(cpm)


# sample-level Pearson correlation ----------------------------------------

cpm_for_pearson <- cpm
count_cormatrix <- round(cor(log2(cpm_for_pearson[, str_detect(colnames(cpm_for_pearson), "fed|gene_id", negate = T)] + 1),
                             use = "all.obs",
                             method = "pearson"),
                         digits = 3)
melted_cormatrix <- melt(count_cormatrix)

# plot Pearson correlation heatmap alphabetically
count_cormatrix_ordered <- count_cormatrix[order(rownames(count_cormatrix)), order(colnames(count_cormatrix))]

melted_cormatrix_ordered <- melt(count_cormatrix_ordered)

p_maybe <-
  ggplot(data = melted_cormatrix_ordered, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "lightgray",
                       high = "black",
                       mid = "#3f3f3f",
                       midpoint = 0.98,
                       limit = c(0.96, 1)) +
  geom_text(aes(Var2, Var1, label = value),
            color = "white",
            size = 3) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.title = element_blank(),
        aspect.ratio = 1) +
  labs(title = "Pearson correlation matrix alphabetically ordered")

p_maybe


# sample-level PCA --------------------------------------------------------

# get a longer data frame for CPM
cpm_long <-
  cpm %>% 
  pivot_longer(-gene_id , names_to = "group", values_to = "cpm") %>%
  as.data.frame()

# pca using covariance matrix
cpm_for_pca <- 
  cpm_long %>% 
  as_tibble() %>%
  pivot_wider(names_from = "group", values_from = "cpm") %>%
  as.data.frame() %>%
  select(-c("gene_id"))

# remove N2_fed
cpm_for_pca_raw <- cpm_for_pca
cpm_for_pca <- cpm_for_pca[, str_detect(colnames(cpm_for_pca), "fed", negate = T)]

# Amy's way of transforming the data
cpm_for_pca$mean <- rowMeans(cpm_for_pca)
tmp <- cpm_for_pca[, 1:12]/cpm_for_pca$mean
pca_using_cov <- prcomp(t(log2(tmp + 1)),
                        center = TRUE,
                        scale. = FALSE)
summary(pca_using_cov)
rm(list = c("tmp"))

# get the score matrix
pca.cov.df <- as.data.frame(pca_using_cov$x)
pca.cov.df <- cbind(pca.cov.df, rownames(pca.cov.df))

colnames(pca.cov.df)[13] <- "genotype"
pca.cov.df$genotype %<>% as.character()
pca.cov.df$genotype <- substr(pca.cov.df$genotype, start = 1, stop = nchar(pca.cov.df$genotype) - 9)

pca.cov.df$rep_raw <- c(2,3,3,4,
                        1,5,
                        1,1,1,2,2,2)

pca.cov.df$rep <- c("B","C","C","C",
                    "A","C",
                    "A","A","A","B","B","B")

# plot PCA
pca.cov.df$genotype %<>% as.character()
pca.cov.df$genotype %<>%
  str_replace_all("N2_starv", "wild type") %>%
  str_replace_all("daf-16", "daf-16(mu86)") %>%
  str_replace_all("daf-18", "daf-18(ok480)") %>%
  str_replace_all("double", "daf-16(mu86); daf-18(ok480)")

pca.cov.df$genotype %<>% as.factor() %>% fct_relevel(c("wild type",
                                                       "daf-16(mu86)",
                                                       "daf-18(ok480)",
                                                       "daf-16(mu86); daf-18(ok480)"))

pca.cov.df$genotype %>% levels()
pca.cov.df$rep %<>% as.character()

p2 <-
  ggplot(pca.cov.df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = genotype, shape = rep), size = 3) +
  scale_color_manual(values = c(color_wt, color_d16, color_d18, color_d16_d18)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "B",
       title = "PCA using cov matrix",
       x = "PC1 (26% of variance)",
       y = "PC2 (18% of variance)")

p2


# number of DEGs triangle -------------------------------------------------

d16_vs_N2 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% nrow()
d18_vs_N2 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-18_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% nrow()
dbl_vs_N2 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16;daf-18_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% nrow()
d18_vs_d16 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-18_vs_daf-16") %>% filter(FDR < fdr_cutoff) %>% nrow()
dbl_vs_d16 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16;daf-18_vs_daf-16") %>% filter(FDR < fdr_cutoff) %>% nrow()
dbl_vs_d18 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16;daf-18_vs_daf-18") %>% filter(FDR < fdr_cutoff) %>% nrow()

DEG_triangle_matrix <- matrix(nrow = 3, ncol = 3)
colnames(DEG_triangle_matrix) <- c("wild type", "daf-16" , "daf-18") 
rownames(DEG_triangle_matrix) <- c("daf-16" , "daf-18", "daf-16; daf-18") 

DEG_triangle_matrix["daf-16", "wild type"] <- d16_vs_N2
DEG_triangle_matrix["daf-18", "wild type"] <- d18_vs_N2
DEG_triangle_matrix["daf-16; daf-18", "wild type"] <- dbl_vs_N2
DEG_triangle_matrix["daf-18", "daf-16"] <- d18_vs_d16
DEG_triangle_matrix["daf-16; daf-18", "daf-16"] <- dbl_vs_d16
DEG_triangle_matrix["daf-16; daf-18", "daf-18"] <- dbl_vs_d18


# d16 - d18 x Cui 2013 lin-35 targets Venn -------------------------------------

background <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16;daf-18_vs_daf-16")
WBid_background <- background$gene_id

GeneLists <- read_csv("DREAM_lin35_Becky_Tepper_Amanda_gene_lists.csv", col_names = TRUE)
cui_lin35_targets <- GeneLists %>% filter(Term %in% c("down by lin-35(n745), starved L1", "UP by lin-35(n745), starved L1"))

set.seed(1)
sets <- list(d16_minus_d18 = d16_minus_d18, cui_lin35_targets = intersect(cui_lin35_targets$wb_id, WBid_background))
p3 <- plot(euler(sets, shape = "ellipse"), quantities = TRUE)
print(p3)
hypergeo(cui_lin35_targets$wb_id, d16_minus_d18, WBid_background) # 0.7683188


# d16 - d18 x Gal 2021 direct lin-35 targets Venn -------------------------------------

set.seed(1)
sets <- list(d16_minus_d18 = d16_minus_d18, gal_lin35_direct_targets = intersect(gal_lin35_direct_targets, WBid_background))
p4 <- plot(euler(sets, shape = "ellipse"), quantities = TRUE)
print(p4)
hypergeo(gal_lin35_direct_targets, d16_minus_d18, WBid_background) # 0.9515436


# arrange panels ------------------------------------------------------

p <-
  grid.arrange(grobs = lapply(list(p1.2, p2),
                              set_panel_size,
                              width = unit(2, "in"),
                              height = unit(2, "in")),
               widths = c(2.2, 2.2),
               ncol = 2)

cairo_pdf("figS1_AB.pdf", width = 15, height = 5) # width usually needs two more panels' widths than all panels' summed widths
plot_grid(p)
dev.off()

