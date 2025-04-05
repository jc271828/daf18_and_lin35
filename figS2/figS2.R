library(tidyverse)
library(magrittr)
library(readxl)

lin35_fed_rep123 <-
  read_xlsx("LIN35GFP_fed_WB_quantification.xlsx", col_names = TRUE) %>%
  select(c("samples", "rep", "GFP_to_tubulin_ratio")) %>%
  filter(samples != "background")

lin35_fed_rep123$samples %<>% str_replace("AWR58_EtOH", "WT")
lin35_fed_rep123$samples %<>% str_replace("LRB486_EtOH", "ok480")

lin35_fed_rep123$samples %<>% as.character() %>% as.factor() %>% fct_relevel(levels = c("WT", "ok480"))
lin35_fed_rep123$GFP_to_tubulin_ratio %<>% as.numeric()
lin35_fed_rep123$rep %<>% as.character()

# stats
# Bartlett's test p-value is 0.9108. Can pool variance
bartlett.test(formula = GFP_to_tubulin_ratio ~ samples, data = lin35_fed_rep123)
# paired t-test p-value is 0.9372
t.test(x = filter(lin35_fed_rep123, samples == "WT")$GFP_to_tubulin_ratio,
       y = filter(lin35_fed_rep123, samples == "ok480")$GFP_to_tubulin_ratio,
       paired = TRUE,
       var.equal = TRUE)

# plot
ggplot(data = lin35_fed_rep123) +
  geom_point(mapping = aes(x = samples, y = GFP_to_tubulin_ratio, color = samples), size = 3) +
  geom_line(mapping = aes(x = samples, y = GFP_to_tubulin_ratio, group = rep)) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2)) +
  labs(title = "3 reps. GFP::LIN-35 fed samples. Collected 16 hr after bleach.\nt-test p-value 0.9372",
       y = "GFP_to_tubulin_ratio (A.U.)")

  