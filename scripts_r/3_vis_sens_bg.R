library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


cost <- read_csv(here::here("data", "cost.csv"))
cost_d2d <- as.list(setNames(cost$D2D, cost$Item))
cost_mdu <- as.list(setNames(cost$MDU, cost$Item))

cost_all <- as.list(setNames(cost$Vul, cost$Item))


folder <- "main"
sims <- read_csv(here::here("out", folder, "Sim_BgACF_IncrementalAverted.csv"))[-1]


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
    PPV_MDU = sum(ACF_TP) / sum(ACF_Yield),
    PPV_D2D = sum(ACF_Alt_TP) / sum(ACF_Alt_Yield)
  ) %>% 
  group_by(Scale) %>% 
  summarise(across(starts_with("PPV"), mean))
  

sims_trend <- sims %>% 
  select(Scale, Year, IncR, MorR) %>%  
  pivot_longer(-c(Scale, Year)) %>% 
  group_by(Scale, Year, name) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.05),
    U = quantile(value, 0.95)
  ) %>%
  ungroup() %>% 
  mutate(
    Scale = factor(Scale, c("Baseline", "1x", "2x", "4x")),
  )


g_trend <- sims_trend %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U, fill = Scale), alpha = 0.1) + 
  geom_line(aes(x = Year, y = M, colour = Scale)) + 
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
  geom_ribbon(aes(x = Year, ymin = L, ymax = U, fill = Scale), alpha = 0.1) + 
  geom_line(aes(x = Year, y = M, colour = Scale)) + 
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
  geom_ribbon(aes(x = Year, ymin = L, ymax = U, fill = Scale), alpha = 0.1) + 
  geom_line(aes(x = Year, y = M, colour = Scale)) + 
  scale_x_continuous("Year", breaks = seq(2022, 2030, 4)) +
  scale_y_continuous("per 100 000 population", labels = scales::number_format(scale = 1e5)) +
  scale_color_discrete("", labels = c(Baseline="Baseline", 
                                      `1x`="Chennai ACF coverage", `2x`="Coverage x2", `4x`="Coverage x4")) +
  expand_limits(y = 0) +
  guides(fill = guide_none()) +
  labs(subtitle = '(B) Annual Mortality') +
  theme(legend.position = c(0, 0), legend.justification = c(-0.1, -0.1))


stats <- local({
  cost <- cost_all
  
  stats <- read_csv(here::here("out", "main", "Sim_BgACF_ScaleUp_Stats.csv"))[-1]
  
  
  stats0 <- stats %>% 
    filter(Scale == "Baseline") %>% 
    select(Key, Inc0 = Inc, Mor0 = Mor)
  
  
  stats <- stats %>% 
    mutate(
      across(starts_with("ACF_"), function(x) ifelse(is.na(x), 0, x)),
      N_Footfall = ACF_Footfall + ACF_Alt_Footfall,
      N_Screened = ACF_Screened + ACF_Alt_Screened,
      N_Confirmed = ACF_Uti_xpert + ACF_Alt_Uti_xpert,
      N_Vul = ACF_Uti_vul + ACF_Alt_Uti_vul,
      N_Sym = ACF_Uti_sym + ACF_Alt_Uti_sym,
      N_VulSym = ACF_Uti_vs + ACF_Alt_Uti_vs,
      N_CXR = ACF_Uti_cxr + ACF_Alt_Uti_cxr,
      N_Xpert = ACF_Uti_xpert + ACF_Alt_Uti_xpert,
      N_Fl = ACF_Fl + ACF_Alt_Fl,
      N_Sl = ACF_Sl + ACF_Alt_Sl,
      N_TPT = ACF_TPT + ACF_Alt_TPT,
      C_Screened = N_Sym * cost$Sym + N_Vul * cost$Vul + N_VulSym * cost$VSC + N_CXR * cost$CXR,
      C_Confirmed = N_Xpert * cost$Xpert,
      C_Tx = (N_Fl + N_TPT) * cost$Tx_Fl + N_Sl * cost$Tx_Sl,
      C_Total = C_Screened + C_Confirmed + C_Tx,
    ) %>% 
    select(Key, Scale, Inc1 = Inc, Mor1 = Mor,
           N_Footfall, N_Screened, N_Confirmed, starts_with("C_")) %>% 
    left_join(stats0) %>% 
    mutate(
      AvtInc = (Inc0 - Inc1) / Inc0,
      AvtMor = (Mor0 - Mor1) / Mor0,
      Yield = (Inc0 - Inc1) / N_Screened
    )
})  


g_ce <- stats %>% 
  group_by(Scale) %>% 
  rename(Avt = AvtInc, Cost = C_Total) %>% 
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
      Scale == "1x" ~ "Chennai ACF",
      Scale == "2x" ~ "x2",
      Scale == "4x" ~ "x4"
    )
  ) %>% 
  ggplot() +
  geom_pointrange(aes(x = Cost_M, y = Avt_M, ymin = Avt_L, ymax = Avt_U)) +
  geom_linerange(aes(xmin = Cost_L, xmax = Cost_U, y = Avt_M)) +
  geom_segment(aes(x = 0, yend = Avt_M, y = 0, xend = Cost_M), linetype = "dashed") +
  geom_text(aes(y = Avt_M, x = Cost_M, label = Scale, hjust = - 0.2, vjust = 1.5)) +
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", labels = scales::number_format(scale = 1e-6, accuracy = 1)) +
  scale_y_continuous("Incident case averted, %, 2023-2030", labels = scales::percent) +
  expand_limits(x = c(0, 12e6), y = 0)


g_ce


stats %>% 
  group_by(Scale) %>% 
  mutate(Avt = AvtInc * Inc0) %>% 
  rename(Cost = C_Total) %>% 
  summarise(ce = sum(Cost) / sum(Avt)) %>% 
  pull(ce)


tab_yield <- stats %>% 
  mutate(Yield = ifelse(is.na(Yield), 0, Yield)) %>% 
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
    Yield = sprintf("%.1f (%.1f-%.1f)", Yield_M * 1e3, Yield_L * 1e3, Yield_U * 1e3)
  ) %>% 
  select(Scale, AvtInc, AvtMor, Yield)


g_bind <- ggarrange(ggarrange(g_inc, g_mor, nrow=2), 
                    g_ce + labs(subtitle = "(C)"), nrow=1, ncol=2)


ggsave(g_combine, filename = here::here("docs", "figs", "g_bg_acf_combine.png"), width = 6, height = 4.5)

ggsave(g_trend, filename = here::here("docs", "figs", "g_bg_acf_trend.png"), width = 8, height = 4.5)

ggsave(g_ce, filename = here::here("docs", "figs", "g_bg_acf_ce.png"), width = 6, height = 4.5)

ggsave(g_bind, filename = here::here("docs", "figs", "g_bg_acf_bind.png"), width = 10, height = 8)

write_csv(tab_yield, here::here("docs", "tabs", "tab_yield_bg_acf.csv"))

