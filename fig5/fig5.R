
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

color_wt <- "#000000"
color_d18 <- "red"
color_clp1 <- "darkgoldenrod"
color_clp1_d18 <- "red"

color_degron_LIN35 <- "#000000"
color_degron_LIN35_d18 <- "red"


# dot plot of LIN-35 MS peptide intensity before vs after cleavage --------

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

raw$BtoA_log2FC <- log2(raw$B_mean/raw$A_mean)

# both before and after K825 peptide log2FCs fail the normality test
shapiro.test(filter(raw, stop <= 825)$BtoA_log2FC)
shapiro.test(filter(raw, start > 825)$BtoA_log2FC)

# use the non-parametric alternative to one-sample t-test (tests median instead of mean against a given mu value)
wilcox.test(filter(raw, stop <= 825)$BtoA_log2FC, mu = 0, alternative = "two.sided") # p-value = 0.005691
wilcox.test(filter(raw, start > 825)$BtoA_log2FC, mu = 0, alternative = "two.sided") # p-value = 0.6136

# dotplot
raw_dotplot <- raw %>% select(c("start", "stop", "BtoA_log2FC"))
raw_dotplot$peptides_before_K825 <- NA

raw_dotplot$peptides_before_K825[raw_dotplot$stop <= 825] <- TRUE
raw_dotplot$peptides_before_K825[raw_dotplot$start > 825] <- FALSE

which(is.na(raw_dotplot$peptides_before_K825)) # no peptide spans K825

raw_dotplot$peptides_before_K825 %<>% as.character() %>% as.factor() %>% fct_relevel(c("TRUE", "FALSE"))

p2 <-
  ggplot(data = raw_dotplot,
         mapping = aes(x = peptides_before_K825,
                       y = BtoA_log2FC,
                       color = peptides_before_K825)) +
  geom_dotplot(fill = NA,
               binaxis = "y",
               stackdir = "center",
               dotsize = 2,
               binwidth = 0.05) +
  stat_summary(geom = "crossbar",
               color = "black",
               fun = mean,
               linewidth = 1,
               width = 0.3) +
  geom_segment(mapping = aes(x = 0, xend = 3, y = 0, yend = 0),
               color = "red",
               linetype = "dashed",
               linewidth = 1,
               lineend = "round",
               linejoin = "round") +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5), expand = expansion(mult = c(0, 0))) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black")) +
  labs(tag = "B")

p2


# degron::GFP::LIN-35 WB in +/- clp-1 +/- daf-18 --------------------------

# AWR58 is degron::GFP::lin-35; rpl-28p::TIR1
# LRB486 is degron::GFP::lin-35; rpl-28p::TIR1; daf-18(ok480)
# LRB598 is degron::GFP::lin-35; rpl-28p::TIR1; clp-1(tm858); daf-18(ok480)
# LRB599 is degron::GFP::lin-35; rpl-28p::TIR1; clp-1(tm858)

# load data
filename <- "LIN35_GFP_WB_in_daf18_clp1"

WB <- read_xlsx(paste(filename, ".xlsx", sep = ""))

WB %<>% na.omit() %>% .[, c(1:2, 11)]

WB$DAF18 <- NA
WB$DAF18[WB$sample %in% c("AWR58", "LRB599")] <- 1
WB$DAF18[WB$sample %in% c("LRB486", "LRB598")] <- 0

WB$CLP1 <- NA 
WB$CLP1[WB$sample %in% c("AWR58", "LRB486")] <- 1
WB$CLP1[WB$sample %in% c("LRB599", "LRB598")] <- 0

WB$sample %<>% as.factor() %>% fct_relevel(c("AWR58", "LRB486", "LRB599", "LRB598"))
WB$rep %<>% as.character()
WB$rep %<>% as.factor() %>% fct_relevel(unique(WB$rep))

p3 <- 
  ggplot(data = WB,
         mapping = aes(x = sample,
                       y = GFP_tubulin_ratio,
                       color = sample)) +
  geom_point(size = 5) +
  stat_summary(geom = "crossbar",
               fun = mean,
               width = 0.3,
               linewidth = 1) +
  stat_summary(geom = "line", mapping = aes(group = rep), color = "black", linewidth = 1) +
  theme_classic() +
  theme(aspect.ratio = 0.3,
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 0.7),
                     breaks = seq(from = 0, to = 0.7, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "3 reps", y = "LIN-35 protein abundance (normalized to alpha-tubulin)", tag = "C")

p3

# ggsave(filename = "fig5_lineplot.pdf", device = cairo_pdf, width = 15, height = 5)

# stats
# p-value = 0.004598, cannot pool variance
bartlett.test(GFP_tubulin_ratio ~ sample, data = WB)

# paired t-test between AWR58 and LRB486 p-value = 0.01947
t.test(x = filter(WB, sample == "AWR58")$GFP_tubulin_ratio,
       y = filter(WB, sample == "LRB486")$GFP_tubulin_ratio,
       var.equal = FALSE,
       paired = TRUE)

# paired t-test between LRB599 and LRB598 p-value = 0.6426
t.test(x = filter(WB, sample == "LRB599")$GFP_tubulin_ratio,
       y = filter(WB, sample == "LRB598")$GFP_tubulin_ratio,
       var.equal = FALSE,
       paired = TRUE)

# paired t-test between AWR58 and LRB599 p-value = 0.08304
t.test(x = filter(WB, sample == "AWR58")$GFP_tubulin_ratio,
       y = filter(WB, sample == "LRB599")$GFP_tubulin_ratio,
       var.equal = FALSE,
       paired = TRUE)


# clp-1(tm858) x daf-18 SS ------------------------------------------------

# load data
filename <- "clp1_x_daf18_SS"

ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(!(strain == "RB2084" & rep %in% c(1, 3))) %>% filter(!(strain == "RB1509" & rep == 4))

ss %<>% filter(strain %in% c("N2", "tm858", "IC166", "tm858_ok480"))

# change strain names
ss$strain %<>%
  str_replace("N2", "wildtype") %>%
  str_replace("IC166", "daf-18(ok480)")

ss$strain[ss$strain == "tm858_ok480"] <- "clp-1(tm858); daf-18(ok480)"
ss$strain[ss$strain == "tm858"] <- "clp-1(tm858)"

ss$strain %<>% as.factor()
ss$strain %>% levels()

ss$strain %<>% fct_relevel(c("wildtype", "daf-18(ok480)", "clp-1(tm858)", "clp-1(tm858); daf-18(ok480)"))
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
ss[str_detect(ss$strain, "tm858", negate = TRUE), "linetype"] <- 2 # 2 is dashed
ss[str_detect(ss$strain, "tm858"), "linetype"] <- 1 # 1 is solid
ss[str_detect(ss$strain, "wildtype"), "linetype"] <- 1 # 1 is solid

ss$shape <- NA
ss[str_detect(ss$strain, "tm858", negate = TRUE), "shape"] <- 21 # 21 is open circle
ss[str_detect(ss$strain, "tm858"), "shape"] <- 16 # 16 is closed circle
ss[str_detect(ss$strain, "wildtype"), "shape"] <- 16 # 16 is closed circle

ss$strain %>% levels()

print(paste('Average number of scored animals is', round(mean(ss$total), 0)))
print(paste('Standard deviation of scored animals is', round(sd(ss$total), 0)))

p6 <-
  ggplot() +
  geom_point(data = ss,
             aes(x = days_after_bleach,
                 y = normalized_by_first_day,
                 color = strain,
                 shape = shape),
             size = 1) +
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
  scale_color_manual(values = c(color_wt, color_d18, color_clp1, color_clp1_d18)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 25, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, by = 0.2)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "F",
       x = "Duration of starvation (days)",
       y = "Proportion alive normalized by the first day",
       title = "Everything was in virgin S-basal + 0.1% EtOH")

p6

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
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.17, can pool variance.

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


# K541A rescues daf-18 SS sensitivity -------------------------------------

# load data
filename <- "lin35K541A_daf18ok480_SS"

ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(!(strain == "LRB652" & rep == 2)) # don't have enough observations

ss %<>% filter(!(strain %in% c("AWR41", "PHX8778") & rep %in% c(1, 2))) # don't have enough observations

ss %<>% filter(!(strain == "PHX8778" & rep == 4)) # dies off weirdly fast

ss %<>% filter(!(strain == "LRB652" & rep %in% c(1, 4))) # goodness of fit isn't good

# change strain names
ss$strain %<>%
  str_replace("AWR41", "degron::GFP::lin-35") %>%
  str_replace("PHX8778", "degron::GFP::lin-35(K541A)") %>%
  str_replace("LRB659", "degron::GFP::lin-35; daf-18(ok480)") %>%
  str_replace("LRB652", "degron::GFP::lin-35(K541A); daf-18(ok480)")

ss$strain %<>% as.factor()
ss$strain %>% levels()

ss$strain %<>% fct_relevel(c("degron::GFP::lin-35", "degron::GFP::lin-35(K541A)", "degron::GFP::lin-35; daf-18(ok480)", "degron::GFP::lin-35(K541A); daf-18(ok480)"))
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
ss[str_detect(ss$strain, "K541A", negate = TRUE), "linetype"] <- 2 # 2 is dashed
ss[str_detect(ss$strain, "K541A"), "linetype"] <- 1 # 1 is solid

ss$shape <- NA
ss[str_detect(ss$strain, "K541A", negate = TRUE), "shape"] <- 21 # 21 is open circle
ss[str_detect(ss$strain, "K541A"), "shape"] <- 16 # 16 is closed circle

ss$strain %>% levels()

print(paste('Average number of scored animals is', round(mean(ss$plated), 0)))
print(paste('Standard deviation of scored animals is', round(sd(ss$plated), 0)))

p7 <-
  ggplot() +
  geom_point(data = ss,
             aes(x = days_after_bleach,
                 y = normalized_by_first_day,
                 color = strain,
                 shape = shape),
             size = 1) +
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
  scale_color_manual(values = c(color_degron_LIN35, color_degron_LIN35, color_degron_LIN35_d18, color_degron_LIN35_d18)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 25, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, by = 0.2)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "G",
       x = "Duration of starvation (days)",
       y = "Proportion alive normalized by the first day",
       title = "Everything was in virgin S-basal + 0.1% EtOH")

p7

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
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.3414, can pool variance.

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


# arrange panels ----------------------------------------------------------

m <- matrix(c(NA, 1, 2, NA, NA, 3, 4, NA, NA), ncol = 3, byrow = TRUE)

p <-
  grid.arrange(grobs = lapply(list(p2, p3, p6, p7),
                              set_panel_size,
                              width = unit(4, "in"),
                              height = unit(4, "in")),
               widths = c(5, 5, 5),
               layout_matrix = m,
               ncol = 3)

cairo_pdf("fig5.pdf", width = 25, height = 25) # width usually needs two more panels' widths than all panels' summed widths
plot_grid(p)
dev.off()

p <-
  grid.arrange(grobs = lapply(list(p6, p7),
                              set_panel_size,
                              width = unit(2, "in"),
                              height = unit(2, "in")),
               widths = c(2.2, 2.2),
               ncol = 2)

cairo_pdf("fig5_SS.pdf", width = 15, height = 5) # width usually needs two more panels' widths than all panels' summed widths
plot_grid(p)
dev.off()
