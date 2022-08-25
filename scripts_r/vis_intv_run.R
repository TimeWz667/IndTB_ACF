library(tidyverse)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


ds <- dir("out")
ds <- ds[startsWith(ds, "dy_")]


sims17 <- read_csv(here::here("out", d, "Runs_IntvE.csv")) %>% 
  mutate(Year = Time)

stats17 <- sims17 %>% 
  filter(Scenario != "Baseline") %>% 
  select(Year, MorR, IncR, Scenario) %>% 
  extract(Scenario, c("ACF", "Coverage", "Type"), 
          "(\\w+), Budget=(\\S+), (Risk|Universal)", convert = T) %>% 
  filter(Type == "Risk") %>% 
  group_by(Year, Coverage) %>% 
  summarise(across(c(IncR, MorR), list(
    M = mean,
    L = function(x) quantile(x, 0.025),
    U = function(x) quantile(x, 0.975)
  ))) %>% 
  pivot_longer(-c(Year, Coverage)) %>% 
  separate(name, c("Index", "Stat"), "_") %>% 
  pivot_wider(names_from = Stat, values_from = value)




sims <- bind_rows(lapply(ds, function(d) {
  load(here::here("out", d, "Acc_Intv.rdata"))
  
  ss <- sims_acc %>% left_join(pars %>% select(Key, p_comorb, rr_risk_comorb))
  
  ss0 <- ss %>% 
    filter(Scenario == "Baseline") %>% 
    select(Key, Inc0 = Inc, Mor0 = Mor)
  ss1 <- ss %>% filter(Scenario != "Baseline")
  
  ss1 %>% 
    extract(Scenario, c("ACF", "Coverage", "Type"), 
            "(\\w+), Budget=(\\S+), (Risk|Universal)", convert = T) %>% 
    left_join(ss0) %>% 
    mutate(gp=d)
}))


g_cost <- sims %>% 
  mutate(
    AvtInc = (Inc0 - Inc) / Inc0,
    AvtMor = (Mor0 - Mor) / Mor0,
    p_comorb = scales::percent(p_comorb, accuracy = .01),
    gp = case_when(
      gp == "dy_free" & Type == "Universal" ~ "All",
      gp != "dy_free" & Type == "Risk" ~ p_comorb,
      T ~ "Other"
    )
  ) %>% 
  filter(gp != "Other") %>% 
  group_by(gp, Coverage) %>% 
  summarise(
    across(starts_with("Avt"), mean),
    Cost_TB = mean(Reached_TB) * 17.53,
    Cost_NonTB = mean(Reached_NonTB) * 17.53,
    Cost = Cost_TB + Cost_NonTB
  ) %>% 
  ungroup() %>% 
  filter(Cost < 21e6) %>% 
  ggplot() +
  geom_line(aes(x = Cost, y = AvtInc, colour = gp)) +
  scale_color_discrete("ACF strategy", 
                       labels = c("0.70%"="In vulnerable, 0.7%", 
                                  "17.85%"="In vulnerable, 17.9%", 
                                  "All"="Untargeted")) +
  scale_y_continuous("Percentage case averted, 2023-2030", labels = scales::percent) + 
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", 
                     labels = scales::number_format(scale=1e-6, accuracy = 1))



g_ppv <- sims %>% 
  mutate(
    p_comorb = scales::percent(p_comorb, accuracy = .01),
    gp = case_when(
      gp == "dy_free" & Type == "Universal" ~ "All",
      gp != "dy_free" & Type == "Risk" ~ p_comorb,
      T ~ "Other"
    )
  ) %>% 
  filter(gp != "Other") %>% 
  group_by(gp, Coverage) %>% 
  summarise(
    PPV = mean(PPV_Acf),
    Cost_TB = mean(Reached_TB) * 17.53,
    Cost_NonTB = mean(Reached_NonTB) * 17.53,
    Cost = Cost_TB + Cost_NonTB
  ) %>% 
  ungroup() %>% 
  filter(Cost < 21e6 & Cost > 0) %>% 
  ggplot() +
  geom_line(aes(x = Cost, y = PPV, colour = gp)) +
  scale_color_discrete("ACF strategy", 
                       labels = c("0.70%"="In vulnerable, 0.7%", 
                                  "17.85%"="In vulnerable, 17.9%", 
                                  "All"="Untargeted")) +
  scale_y_continuous("Positive predictive value among ACF cases, 2030", labels = scales::percent, limits = c(0, 1)) + 
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", 
                     labels = scales::number_format(scale=1e-6, accuracy = 1))


g_spent <- sims %>% 
  mutate(
    p_comorb = scales::percent(p_comorb, accuracy = .01),
    gp = case_when(
      gp == "dy_free" & Type == "Universal" ~ "All",
      gp != "dy_free" & Type == "Risk" ~ p_comorb,
      T ~ "Other"
    )
  ) %>% 
  filter(gp != "Other") %>% 
  group_by(gp, Coverage) %>% 
  summarise(
    Cost_TB = mean(Reached_TB) * 17.53,
    Cost_NonTB = mean(Reached_NonTB) * 17.53,
    Cost = Cost_TB + Cost_NonTB,
    PrTB = Cost_TB / Cost
  ) %>% 
  ungroup() %>% 
  filter(Cost < 21e6 & Cost > 0) %>% 
  ggplot() +
  geom_line(aes(x = Cost, y = PrTB, colour = gp)) +
  scale_color_discrete("ACF strategy", 
                       labels = c("0.70%"="In vulnerable, 0.7%", 
                                  "17.85%"="In vulnerable, 17.9%", 
                                  "All"="Untargeted")) +
  scale_y_continuous("ACF expenditure on True TB, 2030", labels = scales::percent) + 
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", 
                     labels = scales::number_format(scale=1e-6, accuracy = 1)) +
  expand_limits(y = 0)


covs <- unique(stats17$Coverage)[c(1, 6, 13, 21)]
g_trend <- stats17 %>% 
  filter(Coverage %in% covs) %>% 
  mutate(
    Coverage = case_when(
      Coverage == 0 ~ "Baseline",
      Coverage <= covs[2] ~ "Chennai ACF coverage",
      Coverage <= covs[3] ~ "Coverage X2",
      T ~ "Coverage X4"
    )
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U, fill = Coverage), alpha = 0.1) + 
  geom_line(aes(x = Year, y = M, colour = Coverage)) + 
  facet_wrap(.~Index, scale = "free_y", labeller = labeller(Index = c(IncR="Incidence", MorR="Mortality"))) +
  scale_y_continuous("rate per 100 000", labels = scales::number_format(scale = 1e5)) +
  scale_x_continuous("Year", breaks = c(2022, 2025, 2030)) +
  expand_limits(y = 0)




ggsave(g_cost, filename = here::here("docs", "figs", "g_cost.png"), width = 7, height = 4.5)
ggsave(g_ppv, filename = here::here("docs", "figs", "g_ppv.png"), width = 7, height = 4.5)
ggsave(g_spent, filename = here::here("docs", "figs", "g_spent.png"), width = 7, height = 4.5)
ggsave(g_trend, filename = here::here("docs", "figs", "g_trend.png"), width = 9, height = 4.5)
