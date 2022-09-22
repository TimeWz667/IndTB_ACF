library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


cost <- read_csv(here::here("data", "cost.csv"))
cost <- as.list(setNames(cost$M, cost$Item))



folder <- "dy_hi"
sims <- read_csv(here::here("out", folder, "Sim_BgACF_IncrementalAverted.csv"))[-1]



g_combine <- sims %>% 
  pivot_longer(-Key, names_to = "Index") %>% 
  group_by(Index) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  ) %>% 
  ungroup() %>% 
  mutate(
    Index = case_when(
      Index == "Avt_D2D" ~ "D2D screening",
      Index == "Avt_MDU" ~ "MDU screening",
      Index == "Avt_MDU+D2D" ~ "D2D + MDU"
    ),
    Index = factor(Index, c("D2D screening", "MDU screening", "D2D + MDU"))
  ) %>% 
  ggplot() + 
  geom_histogram(aes(x = Index, y = M), stat = "identity") + 
  geom_linerange(aes(x = Index, ymin = L, ymax = U)) +
  scale_y_continuous("Averted incidence, %", labels= scales::percent)


sims %>% 
  mutate(
    p_mdu = Avt_MDU / `Avt_MDU+D2D`,
    p_d2d = Avt_D2D / `Avt_MDU+D2D`
  ) %>% 
  summarise(
    p_mdu = mean(p_mdu),
    p_d2d = mean(p_d2d)
  )


sims <- read_csv(here::here("out", folder, "Sim_BgACF_ScaleUp.csv"))[-1]


g_trend <- sims %>% 
  select(Scale, Time, IncR, MorR) %>%  
  pivot_longer(-c(Scale, Time)) %>% 
  group_by(Scale, Time, name) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.05),
    U = quantile(value, 0.95)
  ) %>%
  ungroup() %>% 
  mutate(
    Scale = factor(Scale, c("Baseline", "1x", "2x", "4x")),
    Time = Time - 0.5
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scale), alpha = 0.1) + 
  geom_line(aes(x = Time, y = M, colour = Scale)) + 
  scale_x_continuous("Year", breaks = seq(2022, 2030, 4)) +
  scale_y_continuous("per 100 000 population", labels = scales::number_format(scale = 1e5)) +
  scale_color_discrete("", labels = c(Baseline="Baseline", 
                                      `1x`="Chennai ACF coverage", `2x`="Coverage x2", `4x`="Coverage x4")) +
  expand_limits(y = 0) +
  facet_wrap(.~name, labeller = labeller(name = c(IncR = "Annual incidence rate", 
                                                  MorR = "Annual mortality rate")), scale = "free_y") + 
  guides(fill = guide_none()) +
  theme(legend.position = c(0, 0), legend.justification = c(-0.1, -0.1))

g_trend


stats <- sims %>% 
  mutate(
    Cost = ACF_MDU_Screened * cost$CXR + ACF_D2D_Screened * cost$Sym, 
    Cost = Cost + (ACF_MDU_Confirmed + ACF_D2D_Confirmed) * cost$Xpert,
    Cost = Cost + (ACF_D2D_DS_Fl + ACF_D2D_DR_Fl) * cost$Tx_Fl,
    Cost = Cost + (ACF_MDU_DS_Fl + ACF_MDU_DR_Fl) * cost$Tx_Fl,
    Cost = Cost + (ACF_MDU_DR_Sl + ACF_D2D_DR_Sl) * cost$Tx_Sl
    # Cost = Cost + Yield_ACF_MU * C_Tx_Fl + Yield_ACF_D2D * C_Tx_Fl
  ) %>% 
  select(Time, Pop, IncR, MorR, Cost, Key, Scale) %>% 
  group_by(Key, Scale) %>% 
  summarise(
    Cost = sum(Cost * Pop) / Pop[1],
    IncN = sum(IncR * Pop) / Pop[1],
    MorN = sum(MorR * Pop) / Pop[1]
  )


g_ce <- stats %>% 
  filter(Scale != "Baseline") %>% 
  left_join(stats %>% filter(Scale == "Baseline") %>% select(Key, IncN0 = IncN)) %>% 
  mutate(Cost = Cost * 1e6, Avt = (IncN0 - IncN) * 1e6) %>% 
  group_by(Scale) %>% 
  summarise(
    across(c(Cost, Avt), list(
      M = mean,
      L = function(x) quantile(x, 0.05),
      U = function(x) quantile(x, 0.95)
    ))
  ) %>% 
  ungroup() %>% 
  mutate(
    Scale = case_when(
      Scale == "1x" ~ "Channei ACF",
      Scale == "2x" ~ "Channei ACF x2",
      Scale == "4x" ~ "Channei ACF x4"
    )
  ) %>% 
  ggplot() +
  geom_pointrange(aes(x = Avt_M, y = Cost_M, ymin = Cost_L, ymax = Cost_U)) +
  geom_linerange(aes(xmin = Avt_L, xmax = Avt_U, y = Cost_M)) +
  geom_segment(aes(x = 0, xend = Avt_M, y = 0, yend = Cost_M), linetype = "dashed") +
  geom_text(aes(x = Avt_M, y = Cost_M, label = Scale, hjust = - 0.2, vjust = 1.5)) +
  scale_y_continuous("Total ACF cost, in millions of 2019 USD", labels = scales::number_format(scale = 1e-6)) +
  scale_x_continuous("Incident case averted, 2023-2030") +
  expand_limits(x = 0, y = 0)


stats %>% 
  filter(Scale != "Baseline") %>% 
  left_join(stats %>% filter(Scale == "Baseline") %>% select(Key, IncN0 = IncN)) %>% 
  mutate(Cost = Cost * 1e6, Avt = (IncN0 - IncN) * 1e6) %>% 
  group_by(Scale) %>% 
  summarise(ce = mean(Cost) / mean(Avt))





ss <- sims %>% 
  mutate(
    N_Screen = ACF_D2D_Screened + ACF_MDU_Screened
    # N_Screen = ACF_D2D_Confirmed + ACF_MDU_Confirmed
  ) %>% 
  select(Time, Pop, IncR, MorR, N_Screen, Key, Scale) %>% 
  group_by(Key, Scale) %>% 
  summarise(
    IncN = sum(IncR * Pop) / Pop[1],
    MorN = sum(MorR * Pop) / Pop[1],
    N_Screen = sum(N_Screen * Pop) / Pop[1]
  )

tab_yield <- ss %>% 
  filter(Scale != "Baseline") %>% 
  left_join(ss %>% filter(Scale == "Baseline") %>% select(Key, IncN0 = IncN, MorN0 = MorN)) %>% 
  mutate(
    AvtInc = 1 - IncN / IncN0,
    AvtMor = 1 - MorN / MorN0,
    Yield = (IncN0 - IncN) / N_Screen
  ) %>% 
  group_by(Scale) %>% 
  summarise(
    across(c(AvtInc, AvtMor, Yield), list(
      M = median,
      L = function(x) quantile(x, 0.025),
      U = function(x) quantile(x, 0.975)
    ))
  ) %>% 
  mutate(
    AvtInc = sprintf("%s (%s-%s)", scales::percent(AvtInc_M), scales::percent(AvtInc_L), scales::percent(AvtInc_U)),
    AvtMor = sprintf("%s (%s-%s)", scales::percent(AvtMor_M), scales::percent(AvtMor_L), scales::percent(AvtMor_U)),
    Yield = sprintf("%.3f (%.3f-%.3f)", Yield_M, Yield_L, Yield_U)
  ) %>% 
  select(Scale, AvtInc, AvtMor, Yield)



ggsave(g_combine, filename = here::here("docs", "g_bg_acf_combine.png"), width = 6, height = 4.5)

ggsave(g_trend, filename = here::here("docs", "g_bg_acf_trend.png"), width = 8, height = 4.5)

ggsave(g_ce, filename = here::here("docs", "g_bg_acf_ce.png"), width = 6, height = 4.5)

write_csv(tab_yield, here::here("docs", "tab_yield_bg_acf.csv"))

