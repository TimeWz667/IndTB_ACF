library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


cost <- read_csv(here::here("data", "cost.csv"))
cost <- as.list(setNames(cost$M, cost$Item))


ds <- c("dy_hi", "dy_lo")

stats <- bind_rows(lapply(ds, function(d) {
  read_csv(here::here("out", d, "Sim_VulACF_followup_stats.csv"))[-1] %>% 
    mutate(Population = d)
}))


stats0 <- stats %>% 
  filter(Scenario == "Baseline") %>% 
  select(Key, Population, Inc0 = IncR, Mor0 = MorR)


stats <- stats %>% 
  filter(Scenario != "Baseline") %>% 
  mutate(
    N_Screened = ACF_Vul_Screened,
    N_Confirmed = ACF_Vul_Confirmed,
    N_Fl_DS = ACF_Vul_DS_Fl,
    N_Fl_DR = ACF_Vul_DR_Fl,
    N_Sl_DR = ACF_Vul_DR_Sl,
    C_Screened = N_Screened * cost$Vul + cost$CXR,
    C_Confirmed = N_Confirmed * cost$Xpert,
    C_Tx = (N_Fl_DS + N_Fl_DR) * cost$Tx_Fl + (N_Fl_DR + N_Sl_DR) * cost$Tx_Sl,
    C_Total = C_Screened + C_Confirmed + C_Tx,
  ) %>% 
  select(Key, Scenario, Population, Coverage, Pop0, Inc1 = IncR, Mor1 = MorR,
         N_Screened, N_Confirmed, starts_with("C_")) %>% 
  left_join(stats0) %>% 
  mutate(
    AvtInc = (Inc0 - Inc1) / Inc0,
    AvtMor = (Mor0 - Mor1) / Mor0,
    across(starts_with("N_"), function(x) x / Pop0),
    across(starts_with("C_"), function(x) x / Pop0)
  )




g_fu <- stats %>% 
  group_by(Scenario, Population) %>% 
  summarise(
    M = mean(AvtInc),
    L = quantile(AvtInc, 0.025),
    U = quantile(AvtInc, 0.975)
  ) %>% 
  ggplot() +
  geom_histogram(aes(x = Scenario, y = M, fill = Scenario), colour="black", stat = "identity") +
  geom_linerange(aes(x = Scenario, ymin = L, ymax = U)) +
  scale_x_discrete("") +
  scale_y_continuous("Averted cases, %", labels = scales::percent) + 
  scale_fill_brewer("Test frequency\n/Follow-up period", palette = 3, 
                    labels = c(Vul_0_0="No follow-up", Vul_2_6="6mo / 2yr", Vul_3_3="3mo / 3yr", Vul_3_6="6mo / 3yr")) +
  facet_grid(Population~., labeller = labeller(Population=c(dy_lo="Low comorbidity", dy_hi = "High comorbidity"))) +
  theme(axis.text.x = element_blank())


ggsave(g_fu, filename = here::here("docs", "g_fu.png"), width = 6, height = 6)

  
  