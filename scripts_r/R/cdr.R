library(tidyverse)

theme_set(theme_bw())


dirs <- dir("out")
dirs <- dirs[startsWith(dirs, "dy_")]
#
dirs <- "dy_20%_5"

summ <- bind_rows(lapply(dirs, function(d) {
  read.csv(here::here("out", d, "Runs_CDR.csv")) %>% 
    mutate(Baseline = d) %>% 
    extract(Baseline, c("PrComorb", "OR_prev_TB"), "dy_(\\d+\\%)_(\\S+)") %>% 
    mutate(
      PrComorb = factor(PrComorb, c("5%", "10%", "20%"))
    )
})) %>% 
  select(-X) %>% 
  pivot_longer(c(DetACF, DetPCF, DurUT, ACF_Reached)) %>% 
  group_by(R_ACF, ACF, PrComorb, OR_prev_TB, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  mutate(
    ACF = factor(ACF, c("Mod", "High"))
  )



g1 <- summ %>% 
  filter(PrComorb == "20%" & OR_prev_TB == "5") %>% 
  filter(startsWith(name, "Det")) %>% 
  ggplot() +
  geom_bar(aes(x = R_ACF, y = M, fill = name), stat = "identity", position = "stack") +
  scale_y_continuous("Case detection, %", labels = scales::percent, limits = 0:1) +
  scale_x_continuous("ACF attendance, per year", breaks = seq(0, 15, 3)) +
  scale_fill_discrete("Source", labels = c(DetACF="Active case-finding", DetPCF = "Routine service")) +
  facet_grid(. ~ ACF, 
             labeller = labeller(.cols = label_both), scales = "free_y") +
  theme(legend.position = "bottom")


g2 <- summ %>% 
  filter(PrComorb == "20%" & OR_prev_TB == "5") %>% 
  filter(name == "DurUT") %>% 
  ggplot() +
  geom_line(aes(x = R_ACF, y = M)) +
  geom_pointrange(aes(x = R_ACF, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("Duration with untreated TB, year") +
  scale_x_continuous("ACF attendance, per year", breaks = seq(0, 15, 3)) +
  facet_grid(. ~ ACF, 
             labeller = labeller(.cols = label_both), scales = "free_y") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")



g3 <- summ %>% 
  filter(PrComorb == "20%" & OR_prev_TB == "5") %>% 
  filter(name == "ACF_Reached") %>% 
  ggplot() +
  geom_line(aes(x = R_ACF, y = M)) +
  geom_pointrange(aes(x = R_ACF, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("ACF coverage, % per active TB", labels = scales::percent) +
  scale_x_continuous("ACF attendance, per year", breaks = seq(0, 15, 3)) +
  facet_grid(. ~ ACF, 
             labeller = labeller(.cols = label_both), scales = "free_y") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")



ggsave(g1, filename = here::here("docs", "CDR", "CaseDetection.png"), width = 7, height = 5)
ggsave(g2, filename = here::here("docs", "CDR", "DelayReduction.png"), width = 7, height = 5)
ggsave(g3, filename = here::here("docs", "CDR", "Coverage.png"), width = 7, height = 5)

g1
g2
g3

