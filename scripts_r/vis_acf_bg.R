library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


cost <- read_csv(here::here("data", "cost.csv"))
cost <- setNames(cost$M, cost$Item)



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
      Index == "Avt_D2D" ~ "DTD screening",
      Index == "Avt_MU" ~ "MDU screening",
      Index == "Avt_MU+D2D" ~ "Both"
    ),
    Index = factor(Index, c("MDU screening", "DTD screening", "Both"))
  ) %>% 
  ggplot() + 
  geom_histogram(aes(x = Index, y = M), stat = "identity") + 
  geom_linerange(aes(x = Index, ymin = L, ymax = U)) +
  scale_y_continuous("Averted incidence, %", labels= scales::percent)



sims <- read_csv(here::here("out", folder, "Sim_BgACF_ScaleUp.csv"))[-1]


sims_cost <- sims %>% 
  select(Key, Scale) %>% 
  distinct() %>% 
  mutate(
    C_CXR = cost_fn$CXR(n()),
    C_Sym = cost_fn$Sym(n()),
    C_Xpert = cost_fn$Xpert(n()),
    C_Tx_Fl = cost_fn$Tx_Fl(n()),
    C_Tx_Sl = cost_fn$Tx_Sl(n())
  )


sims <- sims %>% 
  left_join(sims_cost)


g_trend <- sims %>% 
  group_by(Scale, Time) %>% 
  summarise(
    M = mean(IncR),
    L = quantile(IncR, 0.025),
    U = quantile(IncR, 0.975)
  ) %>%
  ungroup() %>% 
  mutate(
    Scale = factor(Scale, c("Baseline", "1x", "2x", "4x"))
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scale), alpha = 0.1) + 
  geom_line(aes(x = Time, y = M, colour = Scale)) + 
  scale_y_continuous("Incidence, per 100 000", labels = scales::number_format(scale = 1e5)) +
  scale_color_discrete("", labels = c(Baseline="Baseline", 
                                      `1x`="Chennai ACF coverage", `2x`="Coverage x2", `4x`="Coverage x4")) +
  expand_limits(y = 0) +
  guides(fill = guide_none()) +
  theme(legend.position = c(0, 0), legend.justification = c(-0.1, -0.1))

g_trend


stats <- sims %>% 
  left_join(sims_cost) %>% 
  mutate(
    Cost = Reach_ACF_MU1 * C_CXR + Reach_ACF_D2D2 * C_Sym + (Reach_ACF_MU2 + Reach_ACF_D2D2) * C_Xpert
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
      L = function(x) quantile(x, 0.025),
      U = function(x) quantile(x, 0.975)
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



ggsave(g_combine, filename = here::here("out", folder, "g_bg_acf_combine.pdf"), width = 6, height = 4.5)

ggsave(g_trend, filename = here::here("out", folder, "g_bg_acf_trend.pdf"), width = 6, height = 4.5)

ggsave(g_ce, filename = here::here("out", folder, "g_bg_acf_ce.pdf"), width = 6, height = 4.5)



