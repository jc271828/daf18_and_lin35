
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
color_d16 <- "orange"
color_lin35 <- "blue"
color_lin35_d16 <- "purple"

color_d18 <- "red"
color_TIR1 <- "#008080" # teal
color_TIR1_lin35degron <- "darkgoldenrod"
color_TIR1_lin35degron_d18 <- "purple"


# daf-16 x lin-35 SS ------------------------------------------------------

# load and clean up data

filename <- "daf16_lin35_ss"

ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(rep %in% c(1:3) & days_after_bleach > 0)

ss %<>% filter(strain %in% c("N2", "CF1038", "LRB447", "LRB461"))

ss$strain %<>%
  str_replace("N2", "wildtype") %>%
  str_replace("LRB461", "lin-35(n745) daf-16(mu86)") %>%
  str_replace("LRB447", "lin-35(n745)") %>%
  str_replace("CF1038", "daf-16(mu86)")

ss$strain %<>% as.factor()
ss$strain %>% levels()

ss$strain %<>% fct_relevel(c("wildtype", "daf-16(mu86)", "lin-35(n745)", "lin-35(n745) daf-16(mu86)"))
ss$strain %>% levels()

ss[which(ss$proportion > 1), "proportion"] <- 1
ss$normalized_by_first_day <- NA
for (i in 1:nrow(ss)) {
  first_day_proportion <- filter(ss, rep == ss$rep[i] & strain == ss$strain[i]) %>% .$proportion %>% .[1]
  ss$normalized_by_first_day[i] <- ss$proportion[i]/first_day_proportion
}
ss$normalized_by_first_day[ss$normalized_by_first_day > 1] <- 1

ss$strain %>% levels()

print(paste('Average number of scored animals is', round(mean(ss$total), 0)))
print(paste('Standard deviation of scored animals is', round(sd(ss$total), 0)))

# plot normalized proportion
p1 <-
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
  scale_color_manual(values = c(color_wt, color_d16, color_lin35, color_lin35_d16)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 25, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, by = 0.2)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "A",
       x = "Duration of starvation (days)",
       y = "Proportion alive normalized by the first day",
       title = "Everything was in virgin S-basal + 0.1% EtOH")

p1

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
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.2033, can pool variance.

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

# save stats 
# t_test_result %>% write_xlsx(paste(filename, "stats.xlsx", sep = "_"), col_names = TRUE)

# do two-factor test for daf-16 x lin-35
lin35_x_daf16 <- result_raw
lin35_x_daf16$strain %<>% as.character() %>% as.factor()
lin35_x_daf16$strain %>% levels()

lin35_x_daf16$lin35 <- NA
lin35_x_daf16$daf16 <- NA

lin35_x_daf16[str_detect(lin35_x_daf16$strain, "daf-16"), "daf16"] <- "Mut"
lin35_x_daf16[str_detect(lin35_x_daf16$strain, "lin-35"), "lin35"] <- "Mut"
lin35_x_daf16$lin35[is.na(lin35_x_daf16$lin35)] <- "WT"
lin35_x_daf16$daf16[is.na(lin35_x_daf16$daf16)] <- "WT"

# test for homogeneity of variance across groups
leveneTest(data = lin35_x_daf16, half_life ~ lin35*daf16) # Pr 0.3899, equal variance, can do ANOVA

bartlett.test(data = lin35_x_daf16, half_life ~ interaction(lin35, daf16)) # Pr 0.2033, equal variance, can do ANOVA

aov(data = lin35_x_daf16, half_life ~ lin35*daf16) %>% summary() # lin35:daf16 interaction 0.068


# daf-18 x lin-35 AID SS --------------------------------------------------

filename <- "daf18_x_lin35_AID_ss"

# reps 5 to 8 are from 20190614_ss.xlsx
ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(strain != "AWR58_Aux")

# also include lin-35 null (n745) data for comparison
lin35 <- read_xlsx("lin35_ss.xlsx", col_names = TRUE) %>% filter(strain %in%c("LRB447", "AWR58_Aux"))
lin35$rep %<>% paste("-", ., sep = "")

ss %<>% rbind(lin35)
ss$rep %<>% as.numeric()

ss$strain %<>%
  str_replace("IC166", "daf-18(ok480) + EtOH") %>%
  str_replace("LRB447", "lin-35(n745) + EtOH") %>%
  str_replace("rpl28_TIR1", "rpl-28p::TIR1 + IAA") %>%
  str_replace("trpl_EtOH", "degron::lin-35; rpl-28p::TIR1; daf-18(ok480) + EtOH") %>%
  str_replace("trpl_Aux", "degron::lin-35; rpl-28p::TIR1; daf-18(ok480) + IAA") %>%
  str_replace("AWR58_EtOH", "degron::lin-35; rpl-28p::TIR1 + EtOH") %>%
  str_replace("AWR58_Aux", "degron::lin-35; rpl-28p::TIR1 + IAA") %>%
  str_replace("N2", "wildtype + EtOH")

ss$strain %<>% as.factor()
ss$strain %>% levels()

ss$strain %<>%
  fct_relevel(c("wildtype + EtOH",
                "daf-18(ok480) + EtOH",
                "lin-35(n745) + EtOH",
                "rpl-28p::TIR1 + IAA",
                "degron::lin-35; rpl-28p::TIR1 + EtOH",
                "degron::lin-35; rpl-28p::TIR1 + IAA",
                "degron::lin-35; rpl-28p::TIR1; daf-18(ok480) + EtOH",
                "degron::lin-35; rpl-28p::TIR1; daf-18(ok480) + IAA"))

ss[which(ss$proportion > 1), "proportion"] <- 1
ss$normalized_by_first_day <- NA
for (i in 1:nrow(ss)) {
  first_day_proportion <- filter(ss, rep == ss$rep[i] & strain == ss$strain[i]) %>% .$proportion %>% .[1]
  ss$normalized_by_first_day[i] <- ss$proportion[i]/first_day_proportion
}
ss$normalized_by_first_day[ss$normalized_by_first_day > 1] <- 1

ss$strain %>% levels()

ss$linetype <- NA
ss[str_detect(ss$strain, "EtOH"), "linetype"] <- 1 # 1 is solid
ss[str_detect(ss$strain, "IAA"), "linetype"] <- 2 # 2 is dashed

ss$shape <- NA
ss[str_detect(ss$strain, "EtOH"), "shape"] <- 16 # 1 is closed circle
ss[str_detect(ss$strain, "IAA"), "shape"] <- 21 # 2 is open circle

ss$strain %>% levels()

print(paste('Average number of scored animals is', round(mean(ss$total), 0)))
print(paste('Standard deviation of scored animals is', round(sd(ss$total), 0)))

# plot normalized proportion
p2 <-
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
  scale_color_manual(values = c(color_wt, color_d18, color_lin35, color_TIR1, color_TIR1_lin35degron, color_TIR1_lin35degron, color_TIR1_lin35degron_d18, color_TIR1_lin35degron_d18)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 25, by = 5)) +
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
       title = "Everything was in virgin S-basal + either 0.15% EtOH or 200 uM IAA")

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
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.004373, cannot pool variance.

# do unpaired, variance-not-pooled t-tests on half-lives
all_strains <- result$strain %>% as.character() %>% as.factor() %>% levels()
all_pairs <- combinations(all_strains, 2)

t_test_result <- data.frame(strain1 = NA,
                            strain2 = NA,
                            t.test.p.unadjusted = NA) 

for (i in 1: dim(all_pairs)[1]) {
  tmp <- t.test(x = filter(result, strain == all_pairs[i,1])$half_life,
                y = filter(result, strain == all_pairs[i,2])$half_life,
                paired = FALSE,
                var.equal = FALSE)
  t_test_result[i,1] <- all_pairs[i,1]
  t_test_result[i,2] <- all_pairs[i,2]
  t_test_result[i,3] <- tmp$p.value
}
t_test_result %<>% arrange(t.test.p.unadjusted)
t_test_result$p.signif <- stars.pval(t_test_result$t.test.p.unadjusted)

# t_test_result$t.test.p.adjusted %<>% p.adjust(method = "BH")

# t_test_result %>% write_xlsx(paste(filename, "stats.xlsx", sep = "_"), col_names = TRUE)

# do two-factor test for daf-18 x lin-35
result_raw$strain %>% unique()
lin35_x_daf18 <- result_raw %>% filter(strain %in% c("degron::lin-35; rpl-28p::TIR1; daf-18(ok480) + EtOH",
                                                     "degron::lin-35; rpl-28p::TIR1; daf-18(ok480) + IAA",
                                                     "degron::lin-35; rpl-28p::TIR1 + EtOH",
                                                     "degron::lin-35; rpl-28p::TIR1 + IAA"))
lin35_x_daf18$strain %<>% as.character() %>% as.factor()
lin35_x_daf18$strain %>% levels()

lin35_x_daf18$lin35 <- NA
lin35_x_daf18$daf18 <- NA

lin35_x_daf18[str_detect(lin35_x_daf18$strain, "daf-18"), "daf18"] <- "Mut"
lin35_x_daf18[str_detect(lin35_x_daf18$strain, "IAA"), "lin35"] <- "degraded"
lin35_x_daf18$lin35[is.na(lin35_x_daf18$lin35)] <- "WT"
lin35_x_daf18$daf18[is.na(lin35_x_daf18$daf18)] <- "WT"

# test for homogeneity of variance across groups
leveneTest(data = lin35_x_daf18, half_life ~ lin35*daf18) # Pr 0.227, equal variance, can do ANOVA

bartlett.test(data = lin35_x_daf18, half_life ~ interaction(lin35, daf18)) # Pr 0.1701, equal variance, can do ANOVA

aov(data = lin35_x_daf18, half_life ~ lin35*daf18) %>% summary() # lin35:daf18 interaction p 0.00314


# arrange panels ----------------------------------------------------------

p <-
  grid.arrange(grobs = lapply(list(p1, p2),
                              set_panel_size,
                              width = unit(2, "in"),
                              height = unit(2, "in")),
               widths = c(2.2, 2.2),
               ncol = 2)

cairo_pdf("fig3.pdf", width = 15, height = 5) # width usually needs two more panels' widths than all panels' summed widths
plot_grid(p)
dev.off()


