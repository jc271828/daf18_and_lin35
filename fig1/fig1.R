
windowsFonts(sans = windowsFont("Arial"))

# load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(gtools)
library(arrangements)
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
color_akt1 <- "blue"
color_d18_akt1 <- "blue"
color_age1_d18 <- "purple"

color_akt1_gof <- "darkgoldenrod"
color_pdk1_gof <- "#008080" # teal
color_akt1_pdk1_gof <- "darkblue"

color_d16 <- "orange"
color_d16_d18 <- "orange"


# age-1; daf-18 SS ---------------------------------------------------------------

filename <- "age1_daf18_ss"

ss <- read_xlsx(paste(filename, ".xlsx", sep = ""))

ss %<>% filter(!(strain == "LRB441" & rep == 1) & !(strain == "N2" & rep == 2))

ss$strain %<>%
  str_replace("N2", "wild type") %<>%
  str_replace("BQ1", "akt-1(mg306)") %>%
  str_replace("IC166", "daf-18(ok480)") %>%
  str_replace("LRB430", "daf-18(ok480); akt-1(mg306)") %>%
  str_replace("LRB441", "age-1(m333); daf-18(ok480)")

ss$strain %<>% as.factor()
ss$strain %>% levels()
ss$strain %<>% fct_relevel(c("wild type",
                             "daf-18(ok480)",
                             "akt-1(mg306)",
                             "daf-18(ok480); akt-1(mg306)",
                             "age-1(m333); daf-18(ok480)"))

ss[which(ss$proportion > 1), "proportion"] <- 1
ss$normalized_by_first_day <- NA
for (i in 1:nrow(ss)) {
  first_day_proportion <- filter(ss, rep == ss$rep[i] & strain == ss$strain[i]) %>% .$proportion %>% .[1]
  ss$normalized_by_first_day[i] <- ss$proportion[i]/first_day_proportion
}
ss$normalized_by_first_day[ss$normalized_by_first_day > 1] <- 1

ss$strain %>% levels()

ss$linetype <- NA
ss[str_detect(ss$strain, "ok480", negate = TRUE), "linetype"] <- 1 # 1 is solid
ss[str_detect(ss$strain, "ok480"), "linetype"] <- 2 # 2 is dashed

ss$shape <- NA
ss[str_detect(ss$strain, "ok480", negate = TRUE), "shape"] <- 16 # 16 is closed circle
ss[str_detect(ss$strain, "ok480"), "shape"] <- 21 # 21 is open circle

ss$strain %>% levels()

ss_A <- ss

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
  scale_color_manual(values = c(color_wt, color_d18, color_akt1, color_d18_akt1, color_age1_d18)) +
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
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.2771, can pool variance.

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


# akt-1; pdk-1 GoF SS ---------------------------------------------------------------

filename <- "akt1_pdk1_gof_ss"

ss <- read_csv(paste(filename, ".csv", sep = ""))

ss %<>% filter(strain != "CF1038")

colnames(ss)[6] <- "proportion"

ss$strain %<>%
  str_replace("N2", "wild type") %<>%
  str_replace("IC166", "daf-18(ok480)") %>%
  str_replace("LRB429", "akt-1(mg144GoF); pdk-1(mg142GoF)") %>%
  str_replace("GR1318", "pdk-1(mg142GoF)") %>%
  str_replace("GR1310", "akt-1(mg144GoF)")

ss$strain %<>% as.factor()
ss$strain %>% levels()
ss$strain %<>% fct_relevel(c("wild type",
                             "daf-18(ok480)",
                             "akt-1(mg144GoF)",
                             "pdk-1(mg142GoF)",
                             "akt-1(mg144GoF); pdk-1(mg142GoF)"))

ss[which(ss$proportion > 1), "proportion"] <- 1
ss$normalized_by_first_day <- NA
for (i in 1:nrow(ss)) {
  first_day_proportion <- filter(ss, rep == ss$rep[i] & strain == ss$strain[i]) %>% .$proportion %>% .[1]
  ss$normalized_by_first_day[i] <- ss$proportion[i]/first_day_proportion
}
ss$normalized_by_first_day[ss$normalized_by_first_day > 1] <- 1

ss$strain %>% levels()

ss$linetype <- NA
ss[str_detect(ss$strain, "ok480", negate = TRUE), "linetype"] <- 1 # 1 is solid
ss[str_detect(ss$strain, "ok480"), "linetype"] <- 2 # 2 is dashed

ss$shape <- NA
ss[str_detect(ss$strain, "ok480", negate = TRUE), "shape"] <- 16 # 16 is closed circle
ss[str_detect(ss$strain, "ok480"), "shape"] <- 21 # 21 is open circle

ss$strain %>% levels()

ss_B <- ss

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
  scale_color_manual(values = c(color_wt, color_d18, color_akt1_gof, color_pdk1_gof, color_akt1_pdk1_gof)) +
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
  geom_col(mapping = aes(x = rep, y = half_life, color = strain), fill = NA, size = 1)+
  theme(aspect.ratio = 1)+
  theme_classic()+
  geom_text(mapping = aes(x = rep, y = half_life + 1, label = half_life), size = 3)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  facet_wrap(.~strain)

# Bartlett's test
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.4813, can pool variance.

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


# daf-16; daf-18 SS ---------------------------------------------------------------

filename <- "daf16_daf18_ss"

ss <- read_csv(paste(filename, ".csv", sep = ""))

ss %<>% filter(scoring_method == "indirect")
ss %<>% filter(!(strain == "daf18" & rep == 1))

colnames(ss)[8] <- "proportion"

ss$strain %<>%
  str_replace("n2", "wild type") %<>%
  str_replace("daf18", "daf-18(ok480)") %>%
  str_replace("daf16", "daf-16(mu86)") %>%
  str_replace("double", "daf-16(mu86); daf-18(ok480)")

ss$strain %<>% as.factor()
ss$strain %>% levels()
ss$strain %<>% fct_relevel(c("wild type",
                             "daf-16(mu86)",
                             "daf-18(ok480)",
                             "daf-16(mu86); daf-18(ok480)"))

ss[which(ss$proportion > 1), "proportion"] <- 1
ss$normalized_by_first_day <- NA
for (i in 1:nrow(ss)) {
  first_day_proportion <- filter(ss, rep == ss$rep[i] & strain == ss$strain[i]) %>% .$proportion %>% .[1]
  ss$normalized_by_first_day[i] <- ss$proportion[i]/first_day_proportion
}
ss$normalized_by_first_day[ss$normalized_by_first_day > 1] <- 1

ss$strain %>% levels()

ss$linetype <- NA
ss[str_detect(ss$strain, "ok480", negate = TRUE), "linetype"] <- 1 # 1 is solid
ss[str_detect(ss$strain, "ok480"), "linetype"] <- 2 # 2 is dashed

ss$shape <- NA
ss[str_detect(ss$strain, "ok480", negate = TRUE), "shape"] <- 16 # 16 is closed circle
ss[str_detect(ss$strain, "ok480"), "shape"] <- 21 # 21 is open circle

ss$strain %>% levels()

ss_C <- ss

print(paste('Average number of scored animals is', round(mean(ss$total), 0)))
print(paste('Standard deviation of scored animals is', round(sd(ss$total), 0)))

# plot normalized proportion
p3 <-
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
  scale_color_manual(values = c(color_wt, color_d16, color_d18, color_d16_d18)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 25, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, by = 0.2)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  labs(tag = "C",
       x = "Duration of starvation (days)",
       y = "Proportion alive normalized by the first day",
       title = "Everything was in virgin S-basal + 0.1% EtOH")

p3

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
  geom_col(mapping = aes(x = rep, y = half_life, color = strain), fill = NA, size = 1)+
  theme(aspect.ratio = 1)+
  theme_classic()+
  geom_text(mapping = aes(x = rep, y = half_life + 1, label = half_life), size = 3)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  facet_wrap(.~strain)

# Bartlett's test
bartlett.test(half_life ~ strain, data = result) #  p-value = 0.251, can pool variance.

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

# do two-factor test for daf-16 x daf-18
daf18_x_daf16 <- result_raw
daf18_x_daf16$strain %<>% as.character() %>% as.factor()
daf18_x_daf16$strain %>% levels()

daf18_x_daf16$daf18 <- NA
daf18_x_daf16$daf16 <- NA

daf18_x_daf16[str_detect(daf18_x_daf16$strain, "daf-16"), "daf16"] <- "Mut"
daf18_x_daf16[str_detect(daf18_x_daf16$strain, "daf-18"), "daf18"] <- "Mut"
daf18_x_daf16$daf18[is.na(daf18_x_daf16$daf18)] <- "WT"
daf18_x_daf16$daf16[is.na(daf18_x_daf16$daf16)] <- "WT"

# test for homogeneity of variance across groups
leveneTest(data = daf18_x_daf16, half_life ~ daf18*daf16) # Pr 0.8086, equal variance, can do ANOVA

bartlett.test(data = daf18_x_daf16, half_life ~ interaction(daf18, daf16)) # Pr 0.251, equal variance, can do ANOVA

# ANOVA assumes normal distribution and equal variance
aov(data = daf18_x_daf16, half_life ~ daf18*daf16) %>% summary() # daf18:daf16 interaction 6.98e-06


# number of scored animals ------------------------------------------------

print(paste('Average number of scored animals is', round(mean(c(ss_A$total, ss_B$total, ss_C$total)), 0)))
print(paste('Standard deviation of scored animals is', round(sd(c(ss_A$total, ss_B$total, ss_C$total)), 0)))



# arrange panels ----------------------------------------------------------

p <-
  grid.arrange(grobs = lapply(list(p1, p2, p3),
                              set_panel_size,
                              width = unit(2, "in"),
                              height = unit(2, "in")),
               widths = c(2.2, 2.2, 2.2),
               ncol = 3)

cairo_pdf("fig1.pdf", width = 15, height = 5) # width usually needs two more panels' widths than all panels' summed widths
plot_grid(p)
dev.off()




