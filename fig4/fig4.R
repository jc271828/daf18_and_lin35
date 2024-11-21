
windowsFonts(sans = windowsFont("Arial"))

# load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(gtools)
library(arrangements)
library(gplots)
library(reshape2)
library(RColorBrewer)

# for arranging panels
library(cowplot)
library(egg)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)


# define color code -------------------------------------------------------

color_TIR1_lin35degron <- "#566B31"
color_TIR1_lin35degron_d18 <- "#D9A528"

color_wt <- "#000000"
color_d18 <- "#51B848"


# lin-35 WB +/- daf-18 ------------------------------------------------------

# load and clean up data

filename <- "LIN35GFP_WB"

WB <- read_xlsx(paste(filename, ".xlsx", sep = ""))
WB <- WB[4:18, ]

colnames(WB) <- WB[1, ] %>% as.character()
WB <- WB[-1, ]

WB %<>% na.omit() %>% select(c("samples", "rep", "LIN35GFP_IntDen", "GFP_to_tubulin_IntDen"))

WB$samples %<>% as.character()
WB$rep %<>% as.numeric()
WB$LIN35GFP_IntDen %<>% as.numeric()
WB$GFP_to_tubulin_IntDen %<>% as.numeric()

# stats 

# AWR58 is degron::GFP::lin-35; rpl-28p::TIR1
# LRB486 is degron::GFP::lin-35; rpl-28p::TIR1; daf-18(ok480)

# compare plus auxin with background to see if there's still detectable LIN-35 after auxin addition

# p value 0.9698. can pool variance.
bartlett.test(LIN35GFP_IntDen ~ samples, data = filter(WB, samples %in% c("AWR58_IAA", "background")))

# one-sided paired vairance-pooled t-test p-value 0.0454 (alternative = "greater"), there's still some LIN-35 after auxin addition
t.test(x = filter(WB, samples == "AWR58_IAA")$LIN35GFP_IntDen,
       y = filter(WB, samples == "background")$LIN35GFP_IntDen,
       var.equal = TRUE,
       paired = TRUE,
       alternative = "greater")

# compare with and without ok480

# p value 0.8995. can pool variance.
bartlett.test(GFP_to_tubulin_IntDen ~ samples, data = filter(WB, samples %in% c("AWR58_EtOH", "AWR58_ok480_EtOH")))

# p-value comparing plus auxin with background is 0.06815 (alternative = "two.sided"), almost significant -- less LIN-35 in daf-18 mutant
t.test(x = filter(WB, samples == "AWR58_EtOH")$GFP_to_tubulin_IntDen,
       y = filter(WB, samples == "AWR58_ok480_EtOH")$GFP_to_tubulin_IntDen,
       var.equal = TRUE,
       paired = TRUE,
       alternative = "two.sided")

# compare with and without auxin

# p value 0.2038. can pool variance.
bartlett.test(GFP_to_tubulin_IntDen ~ samples, data = filter(WB, samples %in% c("AWR58_IAA", "AWR58_EtOH")))

# p-value comparing plus auxin with background is 0.03445 (alternative = "two.sided")
t.test(x = filter(WB, samples == "AWR58_IAA")$GFP_to_tubulin_IntDen,
       y = filter(WB, samples == "AWR58_EtOH")$GFP_to_tubulin_IntDen,
       var.equal = TRUE,
       paired = TRUE,
       alternative = "two.sided")


# plot
WB$samples %<>% as.factor() %>% fct_relevel(c("background", "AWR58_IAA", "AWR58_EtOH", "AWR58_ok480_EtOH"))

ggplot(data = filter(WB, samples %in% c("background", "AWR58_IAA")),
       mapping = aes(x = samples,
                     y = LIN35GFP_IntDen,
                     color = samples)) +
  geom_point(size = 5) +
  stat_summary(geom = "crossbar",
               fun = mean,
               linewidth = 1,
               width = 0.3) +
  stat_summary(geom = "line", mapping = aes(group = rep), color = "black", linewidth = 1) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_y_continuous(limits = c(0, 2.5),
                     breaks = seq(from = 0, to = 2.5, by = 0.5),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "3 reps")

p2 <- 
  ggplot(data = filter(WB, samples != "background"),
       mapping = aes(x = samples,
                     y = GFP_to_tubulin_IntDen,
                     color = samples)) +
  geom_point(size = 2) +
  stat_summary(geom = "crossbar",
               fun = mean,
               width = 0.3,
               linewidth = 1) +
  stat_summary(geom = "line", mapping = aes(group = rep), color = "black", linewidth = 1) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 0.8),
                     breaks = seq(from = 0, to = 0.8, by = 0.2),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "3 reps", y = "LIN-35 protein abundance (normalized to alpha-tubulin)", tag = "B")

p2


# IP/WB of degron::GFP::LIN-35 --------------------------------------------

filename <- "LIN35GFP_IP_WB"

df <- read_xlsx(paste(filename, ".xlsx", sep = ""))

df <- df[1:16, 1:7]

df_for_ttest <-
  df %>%
  select(rep, `daf-18`, bg_removed_intensity) %>%
  group_by(rep, `daf-18`) %>%
  summarise(lower_to_upper_ratio = bg_removed_intensity[2]/bg_removed_intensity[1])

df_for_ttest$lower_to_upper_ratio # matches what was calculated using Excel

# Bartlett's test
bartlett.test(data = df_for_ttest, lower_to_upper_ratio ~ `daf-18`) # p-value = 0.2418, can pool variance

# paired variance-pooled t-test
t.test(x = df_for_ttest$lower_to_upper_ratio[df_for_ttest$`daf-18` == "WT"],
       y = df_for_ttest$lower_to_upper_ratio[df_for_ttest$`daf-18` == "ok480"],
       paired = TRUE,
       var.equal = TRUE) # p-value = 0.04508

# plot
df_for_ttest$`daf-18` %<>% as.factor() %>% fct_relevel(c("WT", "ok480"))

p3.2 <- 
  ggplot(data = filter(df_for_ttest),
         mapping = aes(x = `daf-18`,
                       y = lower_to_upper_ratio,
                       color = `daf-18`)) +
  geom_point(size = 2) +
  stat_summary(geom = "crossbar",
               fun = mean,
               width = 0.3,
               linewidth = 1) +
  stat_summary(geom = "line", mapping = aes(group = rep), color = "black", linewidth = 1) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  # scale_y_continuous(limits = c(0, 0.8),
  #                    breaks = seq(from = 0, to = 0.8, by = 0.2),
  #                    expand = expansion(mult = c(0, 0))) +
  labs(title = "3 reps", y = "LIN-35 lower band to upper band intensity ratio in LIN-35 IP lysates", tag = "C")

p3.2


# IP/MS peptide intensity ~ coordinates overlaid with Calpain cleavage prediction scores -----------------------------------

filename <- "LIN35GFP_MS"

raw <- read_xlsx(paste(filename, ".xlsx", sep = ""))
raw <- raw[str_detect(raw$Accession, "GFP_Lin35"), ]

raw$All.Proteins %>% unique() # G5EDT1 is LIN-35; P42212 is GFP so everything is indeed unique to GFP_Lin35

# A is AWR58 degron::GFP::lin-35; rpl-28p::TIR1
# B is LRB486 degron::GFP::lin-35; rpl-28p::TIR1; daf-18(ok480)

raw %<>% select(c("Sequence", "Modifications", "101164 A Normalized", "101166 A Normalized", "101168 A Normalized", "101165 B Normalized", "101167 B Normalized", "101169 B Normalized"))

seq <- read_tsv("degron_GFP_lin35_Reinke_Lab_translated_based_on_new_MS_result.txt", col_names = TRUE)
seq$`>LoxP, linker, and introns in GFP (annotated in SUNY Biotech GFP sequences file) are not calculated for translation (total 3744 bp and 1246 aa and 142.8 kDa)` %<>% str_remove_all("\\*")
seq$`>LoxP, linker, and introns in GFP (annotated in SUNY Biotech GFP sequences file) are not calculated for translation (total 3744 bp and 1246 aa and 142.8 kDa)` %>% nchar()
colnames(seq) <- "sequence"

raw$start <- NA
raw$stop <- NA

for (i in 1: nrow(raw)) {
  raw[i, "start"] <- str_locate(seq$sequence, raw$Sequence[i]) %>% .[1, 1]
  raw[i, "stop"] <- str_locate(seq$sequence, raw$Sequence[i]) %>% .[1, 2]
}

raw$A_mean <- rowMeans(raw[, 3:5])
raw$B_mean <- rowMeans(raw[, 6:8])

# plot
raw$BtoA_log2FC <- log2(raw$B_mean/raw$A_mean)

str_locate(seq$sequence, "LKS")
K541inLKS <- "red" # K541 in the GFP-tagged LIN-35 is AA 284 + 541 = 825 (based on new MS result)

p4.1 <-
  ggplot(raw)+
  geom_segment(mapping = aes(x = start, xend = stop, y = BtoA_log2FC, yend = BtoA_log2FC),
               lineend = "round",
               linejoin = "round",
               linewidth = 1) +
  annotate(geom = "segment",
           x = 0,
           xend = 1300,
           y = 0,
           yend = 0,
           color = "red",
           linetype = "dashed",
           linewidth = 1,
           lineend = "round",
           linejoin = "round") +
  annotate(geom = "segment",
           x = 825,
           xend = 825,
           y = -5.5,
           yend = 5.5,
           color = "red",
           linetype = "dashed",
           linewidth = 1,
           lineend = "round",
           linejoin = "round") +
  annotate(geom = "text", label = "K541 in LKS", color = "red", x = 880, y = 3) +
  annotate(geom = "text", label = "log2FC = 0", color = "red", x = 1350, y = 0) +
  annotate(geom = "point", x = 825, y = 0, color = "red", size = 3) +
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"))+
  scale_x_continuous(breaks = seq(from = 0, to = 1250, by = 100),
                     expand = expansion(mult = c(0, 0.15)))+
  scale_y_continuous(limits = c(-6, 6), expand = expansion(mult = c(0, 0))) +
  labs(title = "log2FC of peptide intensity in ok480 vs. WT\n(normalized across all samples) ~ peptide position",
       x = "Position of peptides",
       y = "log2FC(ok480/WT)",
       tag = "D")

p4.1

# load Calpain cleavage site predcition score file (predicted using https://calpcleav.szialab.org/ and degron::GFP::LIN-35 1246 AA sequence)
cleavage_score <- read_tsv("LIN35_calpain_cleavage_prediction_score_20240604.txt")
colnames(cleavage_score) <- c("cleavage_site", "score")
cleavage_score$cleavage_site %<>% floor()

p4.2 <-
  ggplot()+
  geom_line(data = cleavage_score, mapping = aes(x = cleavage_site, y = score), color = "black", linewidth = 0.2) +
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"))+
  scale_x_continuous(breaks = seq(from = 0, to = 1250, by = 100),
                     expand = expansion(mult = c(0, 0.15)))+
  scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0)), breaks = seq(0, 30, by = 5)) +
  annotate(geom = "text", label = "K541 in LKS", color = "red", x = 880, y = 20) +
  annotate(geom = "segment",
           x = 825,
           xend = 825,
           y = 15,
           yend = 30,
           color = "red",
           linetype = "dashed",
           linewidth = 1,
           lineend = "round",
           linejoin = "round") +
  labs(title = "Predicted Calpain cleavage site score ~ peptide position",
       x = "Position of peptides",
       y = "Calpain cleavage prediction scores",
       tag = "D.2")

p4.2


# plot RB1 Calpain cleavage site scores -----------------------------------

cleavage_score <- read_tsv("RB1_calpain_cleavage_prediction_score_20240604.txt")
colnames(cleavage_score) <- c("cleavage_site", "score")
cleavage_score$cleavage_site %<>% floor()

ggplot()+
  geom_line(data = cleavage_score, mapping = aes(x = cleavage_site, y = score), color = "black", linewidth = 0.2) +
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"))+
  scale_x_continuous(breaks = seq(from = 0, to = 928, by = 100),
                     expand = expansion(mult = c(0, 0.15)))+
  scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0)), breaks = seq(0, 30, by = 5)) +
  annotate(geom = "text", label = "K810 in LKS", color = "red", x = 870, y = 20) +
  annotate(geom = "segment",
           x = 810,
           xend = 810,
           y = 20,
           yend = 30,
           color = "red",
           linetype = "dashed",
           linewidth = 1,
           lineend = "round",
           linejoin = "round") +
  labs(title = "Predicted Calpain cleavage site score ~ peptide position",
       x = "Position of peptides",
       y = "Calpain cleavage prediction scores")


# lin-35 overexpression --------------------------------------------------

# load data
filename <- "lin35Ex_SS"

ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(!(rep == 3 & strain == "tm690_ok480"))

ss %<>% filter(strain %in% c("LRB610", "LRB610_GFP_positive", "LRB611", "LRB611_GFP_positive"))

ss %>% filter(strain %in% c("LRB610", "LRB610_GFP_positive")) %>% .$rep %>% unique()

# calculate plated, alive and proportion of GFP(-) populations for LRB610
ss %>% filter(strain %in% c("LRB610", "LRB610_GFP_positive")) %>% .$rep %>% unique()
for(i in 1:3) { # 1:num reps
  
  # LRB610
  LRB610 <- filter(ss, rep == i & strain == "LRB610")
  colnames(LRB610)[4:5] %<>% paste("_all", sep = "")
  
  LRB610_GFP <- filter(ss, rep == i & strain == "LRB610_GFP_positive")
  colnames(LRB610_GFP)[4:5] %<>% paste("_GFP", sep = "")
  
  LRB610_NO_GFP <- merge(LRB610[, 3:5], LRB610_GFP[, 3:5], by = "days_after_bleach", all = T)
  LRB610_NO_GFP %<>% mutate(plated = plated_all - plated_GFP, alive = alive_all - alive_GFP) %>% mutate(proportion = alive/plated)
  
  LRB610_NO_GFP$strain <- "LRB610_GFP_negative"
  LRB610_NO_GFP$rep <- i
  
  LRB610_NO_GFP <- LRB610_NO_GFP[, c(9:10, 1, 6:8)]
  
  ss %<>% rbind(LRB610_NO_GFP)
  
}

# calculate plated, alive and proportion of GFP(-) populations for LRB611
ss %>% filter(strain %in% c("LRB611", "LRB611_GFP_positive")) %>% .$rep %>% unique()
for(i in 1:4) { # 1:num reps
  
  # LRB611
  LRB611 <- filter(ss, rep == i & strain == "LRB611")
  colnames(LRB611)[4:5] %<>% paste("_all", sep = "")
  
  LRB611_GFP <- filter(ss, rep == i & strain == "LRB611_GFP_positive")
  colnames(LRB611_GFP)[4:5] %<>% paste("_GFP", sep = "")
  
  LRB611_NO_GFP <- merge(LRB611[, 3:5], LRB611_GFP[, 3:5], by = "days_after_bleach", all = T)
  LRB611_NO_GFP %<>% mutate(plated = plated_all - plated_GFP, alive = alive_all - alive_GFP) %>% mutate(proportion = alive/plated)
  
  LRB611_NO_GFP$strain <- "LRB611_GFP_negative"
  LRB611_NO_GFP$rep <- i
  
  LRB611_NO_GFP <- LRB611_NO_GFP[, c(9:10, 1, 6:8)]
  
  ss %<>% rbind(LRB611_NO_GFP)
  
}

# change strain names
ss$strain %<>%
  str_replace("LRB610", "kuEx119[lin-35+; sur-5p::GFP]") %>%
  str_replace("LRB611", "daf-18(ok480); kuEx119[lin-35+; sur-5p::GFP]")

# filter out GFP(+) and GFP(-) combined populations
ss %<>% filter(!(strain %in% c("kuEx119[lin-35+; sur-5p::GFP]", "daf-18(ok480); kuEx119[lin-35+; sur-5p::GFP]")))

# factorize strain
ss$strain %<>% as.factor()
ss$strain %>% levels()

ss$strain %<>% fct_relevel(c("kuEx119[lin-35+; sur-5p::GFP]_GFP_negative",
                             "kuEx119[lin-35+; sur-5p::GFP]_GFP_positive",
                             "daf-18(ok480); kuEx119[lin-35+; sur-5p::GFP]_GFP_negative",
                             "daf-18(ok480); kuEx119[lin-35+; sur-5p::GFP]_GFP_positive"))
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

# add a linetype variable
ss$linetype <- NA
ss[str_detect(ss$strain, "positive"), "linetype"] <- 1 # 1 is solid
ss[str_detect(ss$strain, "negative"), "linetype"] <- 2 # 2 is dashed

# plot normalized proportion
ss$strain %>% levels()

p5 <-
  ggplot() +
  geom_point(data = ss,
             aes(x = days_after_bleach,
                 y = normalized_by_first_day,
                 color = strain),
             size = 2) +
  geom_smooth(data = ss,
              aes(x = days_after_bleach,
                  y = normalized_by_first_day,
                  color = strain,
                  linetype = linetype),
              fill = NA,
              method = "glm",
              method.args = list(family = "quasibinomial"),
              level = 0,
              linewidth = 3) +
  scale_linetype_identity() +
  scale_color_manual(values = c(color_wt, color_wt, color_d18, color_d18)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 25, by = 5)) +
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
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.08668, can pool variance.

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

t_test_result %>% write_xlsx(paste(filename, "stats.xlsx", sep = "_"), col_names = TRUE)


# arrange panels ----------------------------------------------------------

m <- matrix(c(NA, 1, NA, 2:4), ncol = 3, byrow = TRUE)

p <-
  grid.arrange(grobs = lapply(list(p2, p3.2, p4, p5),
                              set_panel_size,
                              width = unit(4, "in"),
                              height = unit(4, "in")),
               widths = c(5, 5, 5),
               layout_matrix = m,
               ncol = 3)

cairo_pdf("fig4.pdf", width = 25, height = 25) # width usually needs two more panels' widths than all panels' summed widths
plot_grid(p)
dev.off()


p <-
  grid.arrange(grobs = lapply(list(p2, p3.2),
                              set_panel_size,
                              width = unit(2, "in"),
                              height = unit(2, "in")),
               widths = c(2.2, 2.2),
               ncol = 2)

cairo_pdf("fig4_lineplots.pdf", width = 15, height = 5) # width usually needs two more panels' widths than all panels' summed widths
plot_grid(p)
dev.off()
