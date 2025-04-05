
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
color_lin36 <- "purple"
color_d18_x_lin36 <- "purple"


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

d18 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-18_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% .$gene_id
d16 <- read_xlsx("RNASeq_20200902_d18_d16_dbl_N2Starv_N2Fed_analysis_report.xlsx", sheet = "daf-16_vs_N2_starv") %>% filter(FDR < fdr_cutoff) %>% .$gene_id

d18_minus_d16 <- setdiff(d18, d16)
d16_minus_d18 <- setdiff(d16, d18)

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

gal_efl1_direct_targets <- union(direct_down_in_efl1_gal, direct_up_in_efl1_gal)
gal_dpl1_targets <- union(down_in_dpl1_gal, up_in_dpl1_gal)
goetsch_lin52_3A_targets <- union(down_in_lin52_3A, up_in_lin52_3A)

gal_lin15B_direct_targets <- union(direct_down_in_lin15B_gal, direct_up_in_lin15B_gal)
gal_lin36_direct_targets <- union(direct_down_in_lin36_gal, direct_up_in_lin36_gal)


# do UpSetR -------------------------------------------------------

# d16_minus_d18 %>% as.data.frame() %>% write_tsv("d16_minus_d18.txt", col_names = F)

read_tsv("d16_minus_d18.txt", col_names = F) %>% .$X1 %>% length() # 273
read_tsv("../fig6/DREAM.txt", col_names = F) %>% .$X1 %>% length() # 510
read_tsv("../fig6/efl1_direct.txt", col_names = F) %>% .$X1 %>% length() # 80
read_tsv("../fig6/dpl1.txt", col_names = F) %>% .$X1 %>% length() # 275
read_tsv("../fig6/lin52_3A.txt", col_names = F) %>% .$X1 %>% length() # 105
read_tsv("../fig6/lin15B_direct.txt", col_names = F) %>% .$X1 %>% length() # 143
read_tsv("../fig6/lin36_direct.txt", col_names = F) %>% .$X1 %>% length() # 250

d16_minus_d18_x_DREAM_THAP <- list(d16_minus_d18 = d16_minus_d18,
                                   DREAM = intersect(latorre_8DREAM_bound_L3, WBid_background),
                                   efl1_direct = intersect(gal_efl1_direct_targets, WBid_background),
                                   dpl1 = intersect(gal_dpl1_targets, WBid_background),
                                   lin52_3A = intersect(goetsch_lin52_3A_targets, WBid_background),
                                   lin15B_direct = intersect(gal_lin15B_direct_targets, WBid_background),
                                   lin36_direct = intersect(gal_lin36_direct_targets, WBid_background))

p1 <- upset(fromList(d16_minus_d18_x_DREAM_THAP),
            order.by = "freq",
            keep.order = T,
            # set order is from bottom to top
            sets = c("lin36_direct",
                     "lin15B_direct",
                     "lin52_3A",
                     "dpl1",
                     "efl1_direct",
                     "DREAM",
                     "d16_minus_d18"),
            nintersects = 100)

p1

# hypergeometric tests
hypergeo(latorre_8DREAM_bound_L3, d16_minus_d18, WBid_background) # 0.9957198
hypergeo(gal_efl1_direct_targets, d16_minus_d18, WBid_background) # 1
hypergeo(gal_dpl1_targets, d16_minus_d18, WBid_background) # 0.7399902
hypergeo(goetsch_lin52_3A_targets, d16_minus_d18, WBid_background) # 0.8552915
hypergeo(gal_lin15B_direct_targets, d16_minus_d18, WBid_background) # 0.9283524
hypergeo(gal_lin36_direct_targets, d16_minus_d18, WBid_background) # 1


# lin-36 doesn't affect SS ------------------------------------------------

# load data
filename <- "DREAM_x_daf18_SS"

ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(strain %in% c("JA1850", "N2", "IC166", "LRB638"))

# change strain names
ss$strain %<>%
  str_replace("N2", "wild type") %>%
  str_replace("JA1850", "lin-36(we36)") %>%
  str_replace("IC166", "daf-18(ok480)") %>%
  str_replace("LRB638", "lin-36(we36); daf-18(ok480)")

ss$strain %<>% as.factor()
ss$strain %>% levels()

ss$strain %<>% fct_relevel(c("wild type", "lin-36(we36)", "daf-18(ok480)", "lin-36(we36); daf-18(ok480)"))
ss$strain %>% levels()

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

# plot
p2 <-
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
                  linetype = linetype),
              fill = NA,
              method = "glm",
              method.args = list(family = "quasibinomial"),
              level = 0,
              linewidth = 1) +
  scale_linetype_identity() +
  scale_color_manual(values = c(color_wt, color_lin36, color_d18, color_d18_x_lin36)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 35, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, by = 0.2)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "B",
       x = "Duration of starvation (days)",
       y = "Proportion alive normalized by the first day",
       title = "Everything was in virgin S-basal + 0.1% EtOH")

p2

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
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.7698, can pool variance.

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


# Venn x Fry 2021 high-confidence germline genes --------------------------

set.seed(1)
sets <- list(d16_minus_d18 = d16_minus_d18, fry_germline_genes = fry_germline_genes)
p3 <- plot(euler(sets, shape = "ellipse"), quantities = TRUE)
print(p3)
hypergeo(fry_germline_genes, d16_minus_d18, WBid_background) # 0.9486767


# arrange panels ----------------------------------------------------------

p1
p2

# p <-
#   grid.arrange(grobs = lapply(list(p1, p2),
#                               set_panel_size,
#                               width = unit(2, "in"),
#                               height = unit(2, "in")),
#                widths = c(2.2, 2.2),
#                ncol = 2,
#                nrow = 1)
# 
# cairo_pdf("figS3.pdf", width = 15, height = 5) # width usually needs two more panels' widths than all panels' summed widths
# plot_grid(p)
# dev.off()
