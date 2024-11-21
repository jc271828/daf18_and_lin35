
windowsFonts(sans = windowsFont("Arial"))

# load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(eulerr)
library(UpSetR)
library(gtools)
library(arrangements)
library(dendextend)
library(gplots)
library(reshape2)
library(colorspace)
library(RColorBrewer)
library(car)
library(ARTool)

# for arranging panels
library(cowplot)
library(egg)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)


# define color code -------------------------------------------------------

color_wt <- "#000000"
color_d18 <- "red"

color_dream_bound <- "coral"
color_efl1 <- "blue"
color_efl2 <- "darkgoldenrod"
color_dpl1_neg <- "purple"
color_dpl1_pos <- "darkcyan"
color_lin52 <- "darkgoldenrod"

color_lin15b <- "purple"
color_lin36 <- "blue"

color_d18_x_efl1 <- "blue"
color_d18_x_lin15b <- "purple"


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


# load data ---------------------------------------------------------------

fdr_cutoff <- 0.05

background <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16;daf-18_vs_daf-16")
WBid_background <- background$gene_id

Gal2021_RNASeq_background <- read_xlsx("Gal2021_RNASeq.xlsx", sheet = "lin-35 vs wildtype")$`Wormbase ID` # their lin-15B and efl-1 RNA-seq backgrounds are the same as lin-35 RNA-seq

latorre_background <- read_csv("ws235_gene_names.csv") %>% .$ wb_id

goetsch_lin52_3A_background <- read_xlsx("SupplementalTableS3.xlsx", sheet = "GSE199287") %>% .$WormBase.ID

gal_lin35_direct_pos_targets <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin35_direct_downregulated.txt", col_names = F) %>%
  .$X1
gal_lin35_direct_neg_targets <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin35_direct_upregulated.txt", col_names = F) %>%
  .$X1
gal_lin35_direct_targets <- union(gal_lin35_direct_neg_targets, gal_lin35_direct_pos_targets)

d18 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-18_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% .$gene_id
d16 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% .$gene_id
d18_minus_d16 <- setdiff(d18, d16)

latorre_8DREAM_bound_L3 <- read_tsv("Latorre2015_L3_8DREAM_bound_gene_ids_20240603.txt", col_names = F) %>% .$X1
latorre_8DREAM_bound_Embryo <- read_tsv("Latorre2015_Embryo_8DREAM_bound_gene_ids_20240603.txt", col_names = F) %>% .$X1

direct_down_in_efl1_gal <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_efl1_direct_downregulated.txt", col_names = F) %>%
  .$X1 # has only 1 gene
direct_up_in_efl1_gal <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_efl1_direct_upregulated.txt", col_names = F) %>%
  .$X1

direct_down_in_lin15B_gal <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin15B_direct_downregulated.txt", col_names = F) %>%
  .$X1 # has only 15 genes
direct_up_in_lin15B_gal <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin15B_direct_upregulated.txt", col_names = F) %>%
  .$X1

direct_down_in_lin36_gal <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin36_direct_downregulated.txt", col_names = F) %>%
  .$X1 # has only 5 genes
direct_up_in_lin36_gal <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_lin36_direct_upregulated.txt", col_names = F) %>%
  .$X1

down_in_dpl1_gal <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_dpl1_downregulated.txt", col_names = F) %>%
  .$X1
up_in_dpl1_gal <-
  read_tsv("~/../../OneDrive - Duke University/lab/COMMON_FILES/Gal2021_DREAM_lin35_lin15B_lin36/gal2021_dpl1_upregulated.txt", col_names = F) %>%
  .$X1

DE_in_lin52_3A_vs_WT <-
  read_xlsx("~/../../OneDrive - Duke University/lab/COMMON_FILES/Goetsch2022_LIN52_3A_mutation/SupplementalTableS3.xlsx", sheet = "GSE199287") %>%
  filter(`lin-52.3AvWT.flag` %in% c(1, -1))

up_in_lin52_3A <- DE_in_lin52_3A_vs_WT %>% filter(`lin-52.3AvWT.flag` == 1) %>% .$WormBase.ID
down_in_lin52_3A <- DE_in_lin52_3A_vs_WT %>% filter(`lin-52.3AvWT.flag` == -1) %>% .$WormBase.ID

length(down_in_lin52_3A) # too few, excluded from CDF

fry_germline_genes <- read_tsv("fry2021_HC_germline_genes.txt", col_names = F) %>% .$X1


# discard directionality before Venn --------------------------------------

gal_efl1_direct_targets <- union(direct_down_in_efl1_gal, direct_up_in_efl1_gal) %>% na.omit() %>% as.character()
gal_dpl1_targets <- union(down_in_dpl1_gal, up_in_dpl1_gal) %>% na.omit() %>% as.character()
goetsch_lin52_3A_targets <- union(down_in_lin52_3A, up_in_lin52_3A) %>% na.omit() %>% as.character()

gal_lin15B_direct_targets <- union(direct_down_in_lin15B_gal, direct_up_in_lin15B_gal) %>% na.omit() %>% as.character()
gal_lin36_direct_targets <- union(direct_down_in_lin36_gal, direct_up_in_lin36_gal) %>% na.omit() %>% as.character()


# do UpSetR -------------------------------------------------------

# d18_minus_d16 %>% as.data.frame() %>% write_tsv("d18_minus_d16.txt", col_names = F)
# intersect(latorre_8DREAM_bound_L3$gene_id, WBid_background) %>% as.data.frame() %>% write_tsv("DREAM.txt", col_names = F)
# intersect(gal_efl1_direct_targets, WBid_background) %>% as.data.frame() %>% write_tsv("efl1_direct.txt", col_names = F)
# intersect(gal_dpl1_targets, WBid_background) %>% as.data.frame() %>% write_tsv("dpl1.txt", col_names = F)
# intersect(goetsch_lin52_3A_targets, WBid_background) %>% as.data.frame() %>% write_tsv("lin52_3A.txt", col_names = F)
# intersect(gal_lin15B_direct_targets, WBid_background) %>% as.data.frame() %>% write_tsv("lin15B_direct.txt", col_names = F)
# intersect(gal_lin36_direct_targets, WBid_background) %>% as.data.frame() %>% write_tsv("lin36_direct.txt", col_names = F)

read_tsv("d18_minus_d16.txt", col_names = F) %>% .$X1 %>% length() # 674
read_tsv("DREAM.txt", col_names = F) %>% .$X1 %>% length() # 510
read_tsv("efl1_direct.txt", col_names = F) %>% .$X1 %>% length() # 80
read_tsv("dpl1.txt", col_names = F) %>% .$X1 %>% length() # 275
read_tsv("lin52_3A.txt", col_names = F) %>% .$X1 %>% length() # 105
read_tsv("lin15B_direct.txt", col_names = F) %>% .$X1 %>% length() # 143
read_tsv("lin36_direct.txt", col_names = F) %>% .$X1 %>% length() # 250

d18_minus_d16_x_DREAM_THAP <- list(d18_minus_d16 = d18_minus_d16,
                                   DREAM = intersect(latorre_8DREAM_bound_L3, WBid_background),
                                   efl1_direct = intersect(gal_efl1_direct_targets, WBid_background),
                                   dpl1 = intersect(gal_dpl1_targets, WBid_background),
                                   lin52_3A = intersect(goetsch_lin52_3A_targets, WBid_background),
                                   lin15B_direct = intersect(gal_lin15B_direct_targets, WBid_background),
                                   lin36_direct = intersect(gal_lin36_direct_targets, WBid_background))

p2 <- upset(fromList(d18_minus_d16_x_DREAM_THAP),
            order.by = "freq",
            keep.order = T,
            nintersects = 100,
            # set order is from bottom to top
            sets = c("lin36_direct",
                     "lin15B_direct",
                     "lin52_3A",
                     "dpl1",
                     "efl1_direct",
                     "DREAM",
                     "d18_minus_d16"))

p2

# hypergeometric tests
hypergeo(latorre_8DREAM_bound_L3, d18_minus_d16, WBid_background) # 6.281374e-14
hypergeo(gal_efl1_direct_targets, d18_minus_d16, WBid_background) # 4.506408e-14
hypergeo(gal_dpl1_targets, d18_minus_d16, WBid_background) # 6.570084e-13
hypergeo(goetsch_lin52_3A_targets, d18_minus_d16, WBid_background) # 0.04626565
hypergeo(gal_lin15B_direct_targets, d18_minus_d16, WBid_background) # 0.01216288
hypergeo(gal_lin36_direct_targets, d18_minus_d16, WBid_background) # 1.181249e-43


# plot DREAM-related gene sets in one CDF ---------------------------------

# genes bound by 8 DREAM components
latorre_8DREAM_bound_L3 <- data.frame(gene_id = na.omit(latorre_8DREAM_bound_L3), Term = "DREAM targets")
temp <- background %>% filter(gene_id %in% intersect(latorre_8DREAM_bound_L3$gene_id, WBid_background))
latorre_8DREAM_bound_L3 %<>% merge(temp, ., by = "gene_id", all = F)

# efl-1 direct targets from Gal 2021
gal_efl1_direct_negative_targets <- data.frame(gene_id = na.omit(direct_up_in_efl1_gal), Term = "efl-1 direct negative targets")
temp <- background %>% filter(gene_id %in% intersect(gal_efl1_direct_negative_targets$gene_id, WBid_background))
gal_efl1_direct_negative_targets %<>% merge(temp, ., by = "gene_id", all = F)

# dpl-1 negative targets from Gal 2021
gal_dpl1_neg_targets <- data.frame(gene_id = na.omit(up_in_dpl1_gal), Term = "dpl-1 negative targets")
temp <- background %>% filter(gene_id %in% intersect(gal_dpl1_neg_targets$gene_id, WBid_background))
gal_dpl1_neg_targets %<>% merge(temp, ., by = "gene_id", all = F)

# dpl-1 positive targets from Gal 2021
gal_dpl1_positive_targets <- data.frame(gene_id = na.omit(down_in_dpl1_gal), Term = "dpl-1 positive targets")
temp <- background %>% filter(gene_id %in% intersect(gal_dpl1_positive_targets$gene_id, WBid_background))
gal_dpl1_positive_targets %<>% merge(temp, ., by = "gene_id", all = F)

# lin-52(3A) repressed targets (severing interaction between lin-35 and DREAM) from Goetsch 2022
goetsch_lin52_3A_repressed_targets <- data.frame(gene_id = na.omit(up_in_lin52_3A), Term = "lin-52 repressed targets")
temp <- background %>% filter(gene_id %in% intersect(goetsch_lin52_3A_repressed_targets$gene_id, WBid_background))
goetsch_lin52_3A_repressed_targets %<>% merge(temp, ., by = "gene_id", all = F)

# bind them all
DREAM_efl1_dpl1_lin52 <- rbind(latorre_8DREAM_bound_L3, gal_efl1_direct_negative_targets, gal_dpl1_neg_targets, gal_dpl1_positive_targets, goetsch_lin52_3A_repressed_targets)

DREAM_efl1_dpl1_lin52$Term %>% unique()
DREAM_efl1_dpl1_lin52$Term %<>% as.factor() %>% fct_relevel(c("DREAM targets", "efl-1 direct negative targets", "dpl-1 negative targets", "dpl-1 positive targets", "lin-52 repressed targets"))
DREAM_efl1_dpl1_lin52$Term %>% levels()

# number of genes in each gene set
length(WBid_background) # 15018
DREAM_efl1_dpl1_lin52$Term %>% table()

# plot
p3 <-
  ggplot() +
  stat_ecdf(data = background, aes(logFC), linewidth = 1, color = "black") +
  stat_ecdf(data = DREAM_efl1_dpl1_lin52, aes(logFC, color = Term), linewidth = 1) +
  scale_color_manual(values = c(color_dream_bound, color_efl1, color_dpl1_neg, color_dpl1_pos, color_lin52)) +
  labs(tag = "C", x = "log2FC (daf-16; daf-18/daf-16)", y = "Cumulative proportion") +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = expansion(mult = c(0, 0)))

p3

ks.test(background$logFC,
        filter(latorre_8DREAM_bound_L3, Term == "DREAM targets") %>% .$logFC,
        alternative = "two.sided")$p.value %>%
  signif(2) # 1.3e-09

ks.test(background$logFC,
        filter(gal_efl1_direct_negative_targets, Term == "efl-1 direct negative targets") %>% .$logFC,
        alternative = "two.sided")$p.value %>%
  signif(2) # 0.00065

ks.test(background$logFC,
        filter(gal_dpl1_neg_targets, Term == "dpl-1 negative targets") %>% .$logFC,
        alternative = "two.sided")$p.value %>%
  signif(2) # 1.8e-06

ks.test(background$logFC,
        filter(gal_dpl1_positive_targets, Term == "dpl-1 positive targets") %>% .$logFC,
        alternative = "two.sided")$p.value %>%
  signif(2) # 0.016

ks.test(background$logFC,
        filter(goetsch_lin52_3A_repressed_targets, Term == "lin-52 repressed targets") %>% .$logFC,
        alternative = "two.sided")$p.value %>%
  signif(2) # 0.0037


# plot THAP protein-related gene sets in one CDF --------------------------

# lin-15B direct targets from Gal 2021
gal_lin15B_direct_targets <- data.frame(gene_id = na.omit(direct_up_in_lin15B_gal), Term = "lin-15B direct negative targets")
temp <- background %>% filter(gene_id %in% intersect(gal_lin15B_direct_targets$gene_id, WBid_background))
gal_lin15B_direct_targets %<>% merge(temp, ., by = "gene_id", all = F)

# lin-36 direct negative negative targets from Gal 2021
gal_lin36_direct_neg_targets <- data.frame(gene_id = na.omit(direct_up_in_lin36_gal), Term = "lin-36 direct negative targets")
temp <- background %>% filter(gene_id %in% intersect(gal_lin36_direct_neg_targets$gene_id, WBid_background))
gal_lin36_direct_neg_targets %<>% merge(temp, ., by = "gene_id", all = F)

# bind
lin15B_lin36 <- rbind(gal_lin15B_direct_targets, gal_lin36_direct_neg_targets)

lin15B_lin36$Term %>% unique()
lin15B_lin36$Term %<>% as.factor()
lin15B_lin36$Term %>% levels()

# number of genes in each gene set
length(WBid_background) # 15018
lin15B_lin36$Term %>% table()

# plot
p4 <-
  ggplot() +
  stat_ecdf(data = background, aes(logFC), linewidth = 1, color = "black") +
  stat_ecdf(data = lin15B_lin36, aes(logFC, color = Term), linewidth = 1) +
  scale_color_manual(values = c(color_lin15b, color_lin36)) +
  labs(tag = "D", x = "log2FC (daf-16; daf-18/daf-16)", y = "Cumulative proportion") +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = expansion(mult = c(0, 0)))

p4

ks.test(background$logFC,
        filter(gal_lin15B_direct_targets, Term == "lin-15B direct negative targets") %>% .$logFC,
        alternative = "two.sided")$p.value %>%
  signif(2) # 0.0069

ks.test(background$logFC,
        filter(gal_lin36_direct_neg_targets, Term == "lin-36 direct negative targets") %>% .$logFC,
        alternative = "two.sided")$p.value %>%
  signif(2) # 7.8e-13


# DREAM x daf-18 SS -------------------------------------------------------

# load data
filename <- "DREAM_x_daf18_SS"

ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(!(strain == "LRB635" & rep == 1))

temp <- read_xlsx("DREAM_x_daf18_SS_more.xlsx") %>% filter(condition == "FX02359")

colnames(ss)
colnames(temp)

temp %<>% select(-c("GFP_alive", "prop_GFP_in_survivors"))
colnames(temp)

colnames(temp) <- colnames(ss)

ss %<>% rbind(temp)

# change strain names
ss$strain %<>%
  str_replace("N2", "wildtype") %>%
  str_replace("JA1850", "lin-36(we36)") %>%
  str_replace("MT2495", "lin-15B(n744)") %>%
  str_replace("JJ1549", "efl-1(se1)") %>%
  str_replace("FX02359", "efl-2(tm2359)") %>%
  str_replace("IC166", "daf-18(ok480)") %>%
  str_replace("LRB638", "lin-36(we36); daf-18(ok480)") %>%
  str_replace("LRB635", "daf-18(ok480); lin-15B(n744)") %>%
  str_replace("LRB634", "daf-18(ok480); efl-1(se1)")

# don't plot lin-36 here, plot it in the supp
ss <- ss[str_detect(ss$strain, "we36", negate = TRUE), ]

ss$strain %<>% as.factor()
ss$strain %>% levels()

ss$strain %<>% fct_relevel(c("wildtype", "daf-18(ok480)", "efl-1(se1)", "efl-2(tm2359)", "lin-15B(n744)", "daf-18(ok480); efl-1(se1)", "daf-18(ok480); lin-15B(n744)"))
ss$strain %>% levels()

print(paste('Average number of scored animals is', round(mean(ss$plated), 0)))
print(paste('Standard deviation of scored animals is', round(sd(ss$plated), 0)))

# more tidy-up
ss$rep %<>% as.factor()

ss[which(ss$proportion > 1), "proportion"] <- 1
ss$normalized_by_first_day <- NA
for (i in 1:nrow(ss)) {
  first_day_proportion <- filter(ss, rep == ss$rep[i] & strain == ss$strain[i]) %>% .$proportion %>% .[1]
  ss$normalized_by_first_day[i] <- ss$proportion[i]/first_day_proportion
}
ss$normalized_by_first_day[ss$normalized_by_first_day > 1] <- 1

# plot normalized proportion
ss$strain %>% levels()

ss$linetype <- NA
ss[str_detect(ss$strain, "ok480", negate = TRUE), "linetype"] <- 1 # 1 is solid
ss[str_detect(ss$strain, "ok480"), "linetype"] <- 2 # 2 is dashed

ss$shape <- NA
ss[str_detect(ss$strain, "ok480", negate = TRUE), "shape"] <- 16 # 16 is closed circle
ss[str_detect(ss$strain, "ok480"), "shape"] <- 21 # 21 is open circle

ss$strain %>% levels()

print(paste('Average number of scored animals is', round(mean(ss$plated), 0)))
print(paste('Standard deviation of scored animals is', round(sd(ss$plated), 0)))

p5 <-
  ggplot() +
  geom_point(data = ss,
             aes(x = days_after_bleach,
                 y = normalized_by_first_day,
                 color = strain,
                 shape = shape),
             size = 2) +
  scale_shape_identity() +
  geom_smooth(data = ss,
              aes(x = days_after_bleach,
                  y = normalized_by_first_day,
                  color = strain,
                  linetype =  linetype),
              fill = NA,
              method = "glm",
              method.args = list(family = "quasibinomial"),
              level = 0,
              linewidth = 1) +
  scale_linetype_identity() +
  scale_color_manual(values = c(color_wt, color_d18, color_efl1, color_efl2, color_lin15b, color_d18_x_efl1, color_d18_x_lin15b)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 35, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, by = 0.2)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "E",
       x = "Duration of starvation (days)",
       y = "Proportion alive normalized by the first day",
       title = "Everything was in virgin S-basal + 0.1% EtOH")

p5

# stats
result <- ss %>%
  select(c("strain", "rep")) %>%
  unique() %>%
  as.data.frame()

# delete the conditions where a strain has less than 5 observations
no_less_than_5_obs <- data.frame(obs = NA,
                                 half_life = NA,
                                 strain = result$strain,
                                 rep = result$rep) 
for (i in 1: dim(no_less_than_5_obs)[1]) {
  no_less_than_5_obs[i,1] <- 
    dim(filter(ss,
               strain == no_less_than_5_obs[i,3] &
                 rep == no_less_than_5_obs[i,4]))[1]
}
no_less_than_5_obs %<>% arrange(desc(obs))

# generate a model for each "condition" and calculate stats
for (i in 1:dim(result)[1]) {
  model <- glm(data = filter(ss, strain == result$strain[i] & rep == result$rep[i]),
               formula = normalized_by_first_day ~ days_after_bleach,
               family = "quasibinomial")
  A = model$coefficients[2]
  B = model$coefficients[1]
  fifty_percent_life <- -B/A
  half_life <- -log(2+exp(B), base = exp(A))
  # p < 0.05 for model_sig means there is at least one significant predictor in the model
  model_sig <- pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual, lower.tail = FALSE)
  goodness_of_fit <- 1 - model$deviance/model$null.deviance
  result$fifty_percent_life[i] <- fifty_percent_life
  result$half_life[i] <- half_life
  result$goodness_of_fit[i] <- goodness_of_fit
  result$model_sig[i] <- model_sig
}  
result

# round half_life to 2 digits
result$half_life <- round(result$half_life, 2)
result$rep %<>% as.character()
result_raw <- result

ggplot(result)+
  geom_col(mapping = aes(x = rep, y = half_life, color = strain), fill = NA, linewidth = 1)+
  theme(aspect.ratio = 1)+
  theme_classic()+
  geom_text(mapping = aes(x = rep, y = half_life + 1, label = half_life), size = 3)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  facet_wrap(.~strain)

# Bartlett's test
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.5417, can pool variance.

# do unpaired, variance-pooled t-tests on half-lives
all_strains <- result$strain %>% as.character() %>% as.factor() %>% levels()
all_pairs <- combinations(all_strains, 2)

t_test_result <- data.frame(strain1 = NA,
                            strain2 = NA,
                            t.test.p.unadjusted = NA) 

for (i in 1: dim(all_pairs)[1]) {
  tmp <- t.test(x = filter(result, strain == all_pairs[i,1])$half_life,
                y = filter(result, strain == all_pairs[i,2])$half_life,
                paired = FALSE,
                var.equal = TRUE)
  t_test_result[i,1] <- all_pairs[i,1]
  t_test_result[i,2] <- all_pairs[i,2]
  t_test_result[i,3] <- tmp$p.value
}
t_test_result %<>% arrange(t.test.p.unadjusted)
t_test_result$p.signif <- stars.pval(t_test_result$t.test.p.unadjusted)

# t_test_result$t.test.p.adjusted %<>% p.adjust(method = "BH")

# t_test_result %>% write_xlsx(paste(filename, "stats.xlsx", sep = "_"), col_names = TRUE)


# do two-factor test for daf-18 x efl-1 -----------------------------------

result_raw$strain %>% unique()
efl1_x_daf18 <- result_raw %>% filter(strain %in% c("wildtype",
                                                    "efl-1(se1)",
                                                    "daf-18(ok480)",
                                                    "daf-18(ok480); efl-1(se1)"))
efl1_x_daf18$strain %<>% as.character() %>% as.factor()
efl1_x_daf18$strain %>% levels()

efl1_x_daf18$efl1 <- NA
efl1_x_daf18$daf18 <- NA

efl1_x_daf18[str_detect(efl1_x_daf18$strain, "daf-18"), "daf18"] <- "Mut"
efl1_x_daf18[str_detect(efl1_x_daf18$strain, "efl-1"), "efl1"] <- "Mut"
efl1_x_daf18$daf18[is.na(efl1_x_daf18$daf18)] <- "WT"
efl1_x_daf18$efl1[is.na(efl1_x_daf18$efl1)] <- "WT"

# test for homogeneity of variance across groups
leveneTest(data = efl1_x_daf18, half_life ~ efl1*daf18) # Pr 0.8398, equal variance, can do ANOVA

aov(data = efl1_x_daf18, half_life ~ efl1*daf18) %>% summary() # efl1:daf18 interaction p 0.00201


# do two-factor test for daf-18 x lin-15B ---------------------------------

result_raw$strain %>% unique()
lin15B_x_daf18 <- result_raw %>% filter(strain %in% c("wildtype",
                                                      "lin-15B(n744)",
                                                      "daf-18(ok480)",
                                                      "daf-18(ok480); lin-15B(n744)"))
lin15B_x_daf18$strain %<>% as.character() %>% as.factor()
lin15B_x_daf18$strain %>% levels()

lin15B_x_daf18$lin15B <- NA
lin15B_x_daf18$daf18 <- NA

lin15B_x_daf18[str_detect(lin15B_x_daf18$strain, "daf-18"), "daf18"] <- "Mut"
lin15B_x_daf18[str_detect(lin15B_x_daf18$strain, "lin-15B"), "lin15B"] <- "Mut"
lin15B_x_daf18$daf18[is.na(lin15B_x_daf18$daf18)] <- "WT"
lin15B_x_daf18$lin15B[is.na(lin15B_x_daf18$lin15B)] <- "WT"

# test for homogeneity of variance across groups
leveneTest(data = lin15B_x_daf18, half_life ~ lin15B*daf18) # Pr 0.6578, equal variance, can do ANOVA

aov(data = lin15B_x_daf18, half_life ~ lin15B*daf18) %>% summary() # lin15B:daf18 interaction p 0.00172


# load data for balancer SS -----------------------------------------------

filename <- "DREAM_x_daf18_SS_more"

ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(!(condition %in% c("N2", "FX02359", "LRB372", "MT15109")))

ss$condition %>% unique()

ss$condition %<>%
  str_replace("LRB565", "+/qC1") %>%
  str_replace("LRB596", "+/hT2") %>%
  str_replace("LRB673", "lin-9(n942)/qC1") %>%
  str_replace("LRB683", "+/mIn1") %>%
  str_replace("MT15107", "lin-53(n3368)/hT2") %>%
  str_replace("VC1523", "dpl-1(gk685)/mIn1") %>%
  str_replace("LRB589", "+/qC1; daf-18") %>%
  str_replace("LRB654", "+/hT2; daf-18") %>%
  str_replace("LRB682", "lin-9(n942)/qC1; daf-18") %>%
  str_replace("LRB685", "lin-53(n3368)/hT2; daf-18")

ss$condition %>% unique()

# calculate proportion GFP(-) in survivors
ss$prop_no_GFP_in_survivors <- 1 - ss$prop_GFP_in_survivors

# calculate normalized by Day 1 proportion GFP(-) in survivors
ss$prop_no_GFP_in_survivors_norm_by_day1 <- NA

for (i in 1: nrow(ss)) {
  baseline <- filter(ss, condition == ss$condition[i] & rep == ss$rep[i] & day == 1) %>% .$prop_no_GFP_in_survivors
  ss[i, "prop_no_GFP_in_survivors_norm_by_day1"] <- ss$prop_no_GFP_in_survivors[i]/baseline
} 

ss_raw <- ss

ss$condition %>% unique()


# daf-18 x dpl-1 SS with mIn1 balancer ------------------------------------

ss_no_daf18_mutation <- ss %>% filter(condition %in% c("+/mIn1", "dpl-1(gk685)/mIn1"))
ss_no_daf18_mutation$condition %<>% as.factor() %>% fct_relevel(c("+/mIn1", "dpl-1(gk685)/mIn1"))

print(paste('Average number of scored animals is', round(mean(ss_no_daf18_mutation$plated), 0)))
print(paste('Standard deviation of scored animals is', round(sd(ss_no_daf18_mutation$plated), 0)))

# plot
p6 <-
  ggplot(data = ss_no_daf18_mutation,
         mapping = aes(x = day,
                       y = prop_no_GFP_in_survivors_norm_by_day1,
                       group = condition,
                       color = condition)) +
  geom_point(size = 2,
             # position = position_jitter(),
             na.rm = FALSE,
             # closed circle
             shape = 16) +
  stat_summary(geom = "line",
               linewidth = 1,
               fun = mean,
               # solid line
               linetype = 1) +
  scale_color_manual(values = c("black", "darkgoldenrod")) +
  theme_classic() +
  theme(aspect.ratio = 0.5) +
  scale_x_continuous(breaks = c(1, 4, 8, 16, 20)) +
  # scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
  #                    limits = c(0, 1),
  #                    expand = expansion(mult = c(0, 0.01))) +
  # facet_wrap(.~condition, ncol = 2) +
  labs(title = "Normalized by Day 1\nProportion GFP(-) in survivors\nvirgin S-basal + 0.1% EtOH\n4-6 reps",
       x = "Days of starvation",
       y = "Proportion GFP(-) in survivors (normalized by Day 1)")
p6

# stats

# test for homogeneity of variance for one-way ANOVA
bartlett.test(data = filter(ss_no_daf18_mutation, condition == "+/mIn1"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value < 2.2e-16

# non-parametric test to one-way ANOVA
kruskal.test(data = filter(ss_no_daf18_mutation, condition == "+/mIn1"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value = 0.002753

# test for homogeneity of variance for two-way ANOVA
ss_no_daf18_mutation$day %<>% as.factor() # Levene's test doesn't allow numeric variables as explanatory variables
leveneTest(data = ss_no_daf18_mutation, prop_no_GFP_in_survivors_norm_by_day1 ~ day*condition) # Pr 0.02698

# non-parametric test to two-way ANOVA
art(data = ss_no_daf18_mutation,
    prop_no_GFP_in_survivors_norm_by_day1 ~ condition * day) %>%
  anova() # condition:day p-value 0.03310045


# daf-18 x lin-9 SS with qC1 balancer ------------------------------------

ss_no_daf18_mutation <- ss %>% filter(condition %in% c("+/qC1", "lin-9(n942)/qC1"))
ss_no_daf18_mutation$condition %<>% as.factor() %>% fct_relevel(c("+/qC1", "lin-9(n942)/qC1"))

ss_has_daf18_mutation <- ss %>% filter(condition %in% c("+/qC1; daf-18", "lin-9(n942)/qC1; daf-18"))
ss_has_daf18_mutation$condition %<>% as.factor() %>% fct_relevel(c("+/qC1; daf-18", "lin-9(n942)/qC1; daf-18"))

print(paste('Average number of scored animals is', round(mean(rbind(ss_no_daf18_mutation, ss_has_daf18_mutation)$plated), 0)))
print(paste('Standard deviation of scored animals is', round(sd(rbind(ss_no_daf18_mutation, ss_has_daf18_mutation)$plated), 0)))

# plot
p7.1 <-
  ggplot(data = ss_no_daf18_mutation,
         mapping = aes(x = day,
                       y = prop_no_GFP_in_survivors_norm_by_day1,
                       group = condition,
                       color = condition)) +
  geom_point(size = 2,
             # position = position_jitter(),
             na.rm = FALSE,
             # closed circle
             shape = 16) +
  stat_summary(geom = "line",
               linewidth = 1,
               fun = mean,
               # solid line
               linetype = 1) +
  scale_color_manual(values = c("black", "darkgoldenrod")) +
  theme_classic() +
  theme(aspect.ratio = 0.5) +
  scale_x_continuous(breaks = c(1, 4, 8, 16, 20)) +
  # scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
  #                    limits = c(0, 1),
  #                    expand = expansion(mult = c(0, 0.01))) +
  # facet_wrap(.~condition, ncol = 2) +
  labs(title = "Normalized by Day 1\nProportion GFP(-) in survivors\nvirgin S-basal + 0.1% EtOH\n4-6 reps",
       x = "Days of starvation",
       y = "Proportion GFP(-) in survivors (normalized by Day 1)")
p7.1

p7.2 <-
  ggplot(data = ss_has_daf18_mutation,
         mapping = aes(x = day,
                       y = prop_no_GFP_in_survivors_norm_by_day1,
                       group = condition,
                       color = condition)) +
  geom_point(size = 2,
             # position = position_jitter(),
             na.rm = FALSE,
             # open circle
             shape = 21) +
  stat_summary(geom = "line",
               linewidth = 1,
               fun = mean,
               # dashed line
               linetype = 2) +
  scale_color_manual(values = c("black", "darkgoldenrod")) +
  theme_classic() +
  theme(aspect.ratio = 0.5) +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8)) +
  # scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
  #                    limits = c(0, 1),
  #                    expand = expansion(mult = c(0, 0.01))) +
  # facet_wrap(.~condition, ncol = 2) +
  labs(title = "Normalized by Day 1\nProportion GFP(-) in survivors\nvirgin S-basal + 0.1% EtOH\n4-6 reps",
       x = "Days of starvation",
       y = "Proportion GFP(-) in survivors (normalized by Day 1)")
p7.2

p7 <-
  grid.arrange(grobs = lapply(list(p7.1, p7.2),
                              set_panel_size,
                              width = unit(2, "in"),
                              height = unit(0.7, "in")),
               widths = 2.2,
               ncol = 1)

# stats without daf-18(ok480) mutation

# test for homogeneity of variance for one-way ANOVA
bartlett.test(data = filter(ss_no_daf18_mutation, condition == "+/qC1"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value < 2.2e-16

# non-parametric test to one-way ANOVA
kruskal.test(data = filter(ss_no_daf18_mutation, condition == "+/qC1"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value = 0.05899

# test for homogeneity of variance for two-way ANOVA
ss_no_daf18_mutation$day %<>% as.factor() # Levene's test doesn't allow numeric variables as explanatory variables
leveneTest(data = ss_no_daf18_mutation, prop_no_GFP_in_survivors_norm_by_day1 ~ day*condition) # Pr 0.002562

# non-parametric test to two-way ANOVA
art(data = ss_no_daf18_mutation,
    prop_no_GFP_in_survivors_norm_by_day1 ~ condition * day) %>%
  anova() # condition:day p-value 1.1954e-05

# stats with daf-18(ok480) mutation

# test for homogeneity of variance for one-way ANOVA
bartlett.test(data = filter(ss_has_daf18_mutation, condition == "+/qC1; daf-18"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value < 2.2e-16

# non-parametric test to one-way ANOVA
kruskal.test(data = filter(ss_has_daf18_mutation, condition == "+/qC1; daf-18"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value = 0.2723

# test for homogeneity of variance for two-way ANOVA
ss_has_daf18_mutation$day %<>% as.factor() # Levene's test doesn't allow numeric variables as explanatory variables
leveneTest(data = ss_has_daf18_mutation, prop_no_GFP_in_survivors_norm_by_day1 ~ day*condition) # Pr 0.01943

# non-parametric test to two-way ANOVA
art(data = ss_has_daf18_mutation,
    prop_no_GFP_in_survivors_norm_by_day1 ~ condition * day) %>%
  anova() # condition:day p-value 0.137227


# daf-18 x lin-53 SS with hT2 balancer ------------------------------------

ss_no_daf18_mutation <- ss %>% filter(condition %in% c("+/hT2", "lin-53(n3368)/hT2"))
ss_no_daf18_mutation$condition %<>% as.factor() %>% fct_relevel(c("+/hT2", "lin-53(n3368)/hT2"))

ss_has_daf18_mutation <- ss %>% filter(condition %in% c("+/hT2; daf-18", "lin-53(n3368)/hT2; daf-18"))
ss_has_daf18_mutation$condition %<>% as.factor() %>% fct_relevel(c("+/hT2; daf-18", "lin-53(n3368)/hT2; daf-18"))

print(paste('Average number of scored animals is', round(mean(rbind(ss_no_daf18_mutation, ss_has_daf18_mutation)$plated), 0)))
print(paste('Standard deviation of scored animals is', round(sd(rbind(ss_no_daf18_mutation, ss_has_daf18_mutation)$plated), 0)))

# plot
p8.1 <-
  ggplot(data = ss_no_daf18_mutation,
         mapping = aes(x = day,
                       y = prop_no_GFP_in_survivors_norm_by_day1,
                       group = condition,
                       color = condition)) +
  geom_point(size = 2,
             # position = position_jitter(),
             na.rm = FALSE,
             # closed circle
             shape = 16) +
  stat_summary(geom = "line",
               linewidth = 1,
               fun = mean,
               # solid line
               linetype = 1) +
  scale_color_manual(values = c("black", "darkgoldenrod")) +
  theme_classic() +
  theme(aspect.ratio = 0.5) +
  scale_x_continuous(breaks = c(1, 4, 8, 16, 20)) +
  # scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
  #                    limits = c(0, 1),
  #                    expand = expansion(mult = c(0, 0.01))) +
  # facet_wrap(.~condition, ncol = 2) +
  labs(title = "Normalized by Day 1\nProportion GFP(-) in survivors\nvirgin S-basal + 0.1% EtOH\n4-6 reps",
       x = "Days of starvation",
       y = "Proportion GFP(-) in survivors (normalized by Day 1)")
p8.1

p8.2 <-
  ggplot(data = ss_has_daf18_mutation,
         mapping = aes(x = day,
                       y = prop_no_GFP_in_survivors_norm_by_day1,
                       group = condition,
                       color = condition)) +
  geom_point(size = 2,
             # position = position_jitter(),
             na.rm = FALSE,
             # open circle
             shape = 21) +
  stat_summary(geom = "line",
               linewidth = 1,
               fun = mean,
               # dashed line
               linetype = 2) +
  scale_color_manual(values = c("black", "darkgoldenrod")) +
  theme_classic() +
  theme(aspect.ratio = 0.5) +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8)) +
  # scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
  #                    limits = c(0, 1),
  #                    expand = expansion(mult = c(0, 0.01))) +
  # facet_wrap(.~condition, ncol = 2) +
  labs(title = "Normalized by Day 1\nProportion GFP(-) in survivors\nvirgin S-basal + 0.1% EtOH\n4-6 reps",
       x = "Days of starvation",
       y = "Proportion GFP(-) in survivors (normalized by Day 1)")
p8.2

p8 <-
  grid.arrange(grobs = lapply(list(p8.1, p8.2),
                              set_panel_size,
                              width = unit(2, "in"),
                              height = unit(0.7, "in")),
               widths = 2.2,
               ncol = 1)

# stats without daf-18(ok480) mutation

# test for homogeneity of variance for one-way ANOVA
bartlett.test(data = filter(ss_no_daf18_mutation, condition == "+/hT2"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value < 2.2e-16

# non-parametric test to one-way ANOVA
kruskal.test(data = filter(ss_no_daf18_mutation, condition == "+/hT2"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value = 0.3906

# test for homogeneity of variance for two-way ANOVA
ss_no_daf18_mutation$day %<>% as.factor() # Levene's test doesn't allow numeric variables as explanatory variables
leveneTest(data = ss_no_daf18_mutation, prop_no_GFP_in_survivors_norm_by_day1 ~ day*condition) # Pr 0.01124

# non-parametric test to two-way ANOVA
art(data = ss_no_daf18_mutation,
    prop_no_GFP_in_survivors_norm_by_day1 ~ condition * day) %>%
  anova() # condition:day p-value 0.02290075

# stats with daf-18(ok480) mutation

# test for homogeneity of variance for one-way ANOVA
bartlett.test(data = filter(ss_has_daf18_mutation, condition == "+/hT2; daf-18"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value < 2.2e-16

# non-parametric test to one-way ANOVA
kruskal.test(data = filter(ss_has_daf18_mutation, condition == "+/hT2; daf-18"), prop_no_GFP_in_survivors_norm_by_day1 ~ day) # p-value = 0.8601

# test for homogeneity of variance for two-way ANOVA
ss_has_daf18_mutation$day %<>% as.factor() # Levene's test doesn't allow numeric variables as explanatory variables
leveneTest(data = ss_has_daf18_mutation, prop_no_GFP_in_survivors_norm_by_day1 ~ day*condition) # Pr 0.04521

# non-parametric test to two-way ANOVA
art(data = ss_has_daf18_mutation,
    prop_no_GFP_in_survivors_norm_by_day1 ~ condition * day) %>%
  anova() # condition:day p-value 0.0910533


# CDF of Fry 2021's high-confidence germline genes in daf-16; daf-18 vs daf-16--------------------------------------

p9 <-
  ggplot() +
  stat_ecdf(data = background, aes(logFC), linewidth = 1, color = "black") +
  stat_ecdf(data = fry_germline_genes, aes(logFC, color = Term), linewidth = 1) +
  scale_color_manual(values = "red") +
  labs(tag = "I", x = "log2FC (daf-16; daf-18/daf-16)", y = "Cumulative proportion") +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = expansion(mult = c(0, 0)))

p9

# stats
ks.test(background$logFC,
        fry_germline_genes$logFC,
        alternative = "two.sided")$p.value %>%
  signif(2) # 0 (< 2.22e−16)


#  some hypergeometric tests ----------------------------------------------

# discard directionality
gal_efl1_direct_targets <- union(direct_down_in_efl1_gal, direct_up_in_efl1_gal) %>% na.omit() %>% as.character()
gal_dpl1_targets <- union(down_in_dpl1_gal, up_in_dpl1_gal) %>% na.omit() %>% as.character()
goetsch_lin52_3A_targets <- union(down_in_lin52_3A, up_in_lin52_3A) %>% na.omit() %>% as.character()

gal_lin15B_direct_targets <- union(direct_down_in_lin15B_gal, direct_up_in_lin15B_gal) %>% na.omit() %>% as.character()
gal_lin36_direct_targets <- union(direct_down_in_lin36_gal, direct_up_in_lin36_gal) %>% na.omit() %>% as.character()

# hypergeometric tests
hypergeo(fry_germline_genes, d18_minus_d16, WBid_background) # 8.415059e-07

hypergeo(fry_germline_genes, gal_lin35_direct_targets, Gal2021_RNASeq_background) # 8.893569e-14

hypergeo(fry_germline_genes, latorre_8DREAM_bound_L3$gene_id, latorre_background) # 9.352733e-06

hypergeo(fry_germline_genes, gal_efl1_direct_targets, Gal2021_RNASeq_background) # 0.1335914
hypergeo(fry_germline_genes, gal_dpl1_targets, Gal2021_RNASeq_background) # 4.757746e-08

hypergeo(fry_germline_genes, goetsch_lin52_3A_targets, goetsch_lin52_3A_background) # 1.631594e-05

hypergeo(fry_germline_genes, gal_lin15B_direct_targets, Gal2021_RNASeq_background) # 5.461041e-14
hypergeo(fry_germline_genes, gal_lin36_direct_targets, Gal2021_RNASeq_background) # 0.004511102


# arrange panels ----------------------------------------------------------

p2

p3

p4

p5

p6

plot(p7)

plot(p8)

p9

# m <- matrix(c(NA, 1, NA, 2, NA, 3, NA, 4, NA, 5, 6, NA), ncol = 3, byrow = TRUE)
# 
# p <-
#   grid.arrange(grobs = lapply(list(p2.L3, p4, p6, p8, p10, p11),
#                               set_panel_size,
#                               width = unit(4, "in"),
#                               height = unit(4, "in")),
#                widths = c(5, 5, 5),
#                layout_matrix = m,
#                ncol = 3)
# 
# cairo_pdf("fig6.pdf", width = 25, height = 25) # width usually needs two more panels' widths than all panels' summed widths
# plot_grid(p)
# dev.off()
