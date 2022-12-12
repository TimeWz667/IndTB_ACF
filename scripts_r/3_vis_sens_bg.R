library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


cost <- read_csv(here::here("data", "cost.csv"))
cost_d2d <- as.list(setNames(cost$D2D, cost$Item))
cost_mdu <- as.list(setNames(cost$MDU, cost$Item))


folder <- "main"
sims <- read_csv(here::here("out", folder, "Sim_BgACF_IncrementalAverted.csv"))[-1]


pop.size <- 3e6



g_combine <- sims %>% 
  pivot_longer(-Key, names_to = "Index") %>% 
  group_by(Index) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
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
  scale_x_discrete("") +
  scale_y_continuous("Averted incidence, %", labels= scales::percent)

g_combine


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


sims %>% 
  group_by(Key, Scale) %>% 
  summarise(
    PPV_MDU = sum(ACF_MDU_TP) / sum(ACF_MDU_Yield),
    PPV_D2D = sum(ACF_D2D_TP) / sum(ACF_D2D_Yield)
  ) %>% 
  group_by(Scale) %>% 
  summarise(across(starts_with("PPV"), mean))
  

sims_trend <- sims %>% 
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
  )


g_trend <- sims_trend %>% 
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


g_inc <- sims_trend %>% 
  filter(name == "IncR") %>% 
  ggplot() +
  geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scale), alpha = 0.1) + 
  geom_line(aes(x = Time, y = M, colour = Scale)) + 
  scale_x_continuous("Year", breaks = seq(2022, 2030, 4)) +
  scale_y_continuous("per 100 000 population", labels = scales::number_format(scale = 1e5)) +
  scale_color_discrete("", labels = c(Baseline="Baseline", 
                                      `1x`="Chennai ACF coverage", `2x`="Coverage x2", `4x`="Coverage x4")) +
  expand_limits(y = 0) +
  guides(fill = guide_none()) +
  labs(subtitle = '(A) Annual incidence') +
  theme(legend.position = 'None')


g_mor <- sims_trend %>% 
  filter(name == "MorR") %>% 
  ggplot() +
  geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scale), alpha = 0.1) + 
  geom_line(aes(x = Time, y = M, colour = Scale)) + 
  scale_x_continuous("Year", breaks = seq(2022, 2030, 4)) +
  scale_y_continuous("per 100 000 population", labels = scales::number_format(scale = 1e5)) +
  scale_color_discrete("", labels = c(Baseline="Baseline", 
                                      `1x`="Chennai ACF coverage", `2x`="Coverage x2", `4x`="Coverage x4")) +
  expand_limits(y = 0) +
  guides(fill = guide_none()) +
  labs(subtitle = '(B) Annual Mortality') +
  theme(legend.position = c(0, 0), legend.justification = c(-0.1, -0.1))


stats <- sims %>% 
  mutate(
    Cost = ACF_MDU_CXR * cost_mdu$CXR + ACF_D2D_CXR * cost_d2d$CXR, 
    Cost = Cost + ACF_MDU_Sym * cost_mdu$Sym + ACF_D2D_Sym * cost_d2d$Sym, 
    Cost = Cost + ACF_MDU_Xpert * cost_mdu$Xpert + ACF_D2D_Xpert * cost_d2d$Xpert, 
    Cost = Cost + ACF_MDU_DS_Fl * cost_mdu$Tx_Fl + ACF_D2D_DS_Fl * cost_d2d$Tx_Fl, 
    Cost = Cost + ACF_MDU_DR_Fl * cost_mdu$Tx_Fl + ACF_D2D_DR_Fl * cost_d2d$Tx_Fl, 
    Cost = Cost + ACF_MDU_DR_Sl * cost_mdu$Tx_Sl + ACF_D2D_DR_Sl * cost_d2d$Tx_Sl
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
  mutate(Cost = Cost * pop.size, Avt = (IncN0 - IncN) / IncN0) %>% 
  group_by(Scale) %>% 
  summarise(
    across(c(Cost, Avt), list(
      M = mean,
      L = function(x) quantile(x, 0.25),
      U = function(x) quantile(x, 0.75)
    ))
  ) %>% 
  ungroup() %>% 
  mutate(
    Scale = case_when(
      Scale == "1x" ~ "Chennai ACF",
      Scale == "2x" ~ "Chennai ACF x2",
      Scale == "4x" ~ "Chennai ACF x4"
    )
  ) %>% 
  ggplot() +
  geom_pointrange(aes(x = Cost_M, y = Avt_M, ymin = Avt_L, ymax = Avt_U)) +
  geom_linerange(aes(xmin = Cost_L, xmax = Cost_U, y = Avt_M)) +
  geom_segment(aes(x = 0, yend = Avt_M, y = 0, xend = Cost_M), linetype = "dashed") +
  geom_text(aes(y = Avt_M, x = Cost_M, label = Scale, hjust = - 0.2, vjust = 1.5)) +
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", labels = scales::number_format(scale = 1e-6, accuracy = 1)) +
  scale_y_continuous("Incident case averted, %, 2023-2030", labels = scales::percent) +
  expand_limits(x = c(0, 30e6), y = 0)


g_ce


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


g_bind <- ggarrange(ggarrange(g_inc, g_mor, nrow=2), 
                    g_ce + labs(subtitle = "(C)"), nrow=1, ncol=2)


dggsave(g_combine, filename = here::here("docs", "figs", "g_bg_acf_combine.png"), width = 6, height = 4.5)

ggsave(g_trend, filename = here::here("docs", "figs", "g_bg_acf_trend.png"), width = 8, height = 4.5)

ggsave(g_ce, filename = here::here("docs", "figs", "g_bg_acf_ce.png"), width = 6, height = 4.5)

ggsave(g_bind, filename = here::here("docs", "figs", "g_bg_acf_bind.png"), width = 10, height = 8)

write_csv(tab_yield, here::here("docs", "tabs", "tab_yield_bg_acf.csv"))

