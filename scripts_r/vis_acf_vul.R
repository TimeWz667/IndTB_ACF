library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


pop.size <- 3e6


cost <- read_csv(here::here("data", "cost.csv"))
cost_plain <- as.list(setNames(cost$Plain, cost$Item))
cost <- as.list(setNames(cost$Vul, cost$Item))


ds <- c("dy_hi", "dy_lo")

stats <- bind_rows(lapply(ds, function(d) {
  read_csv(here::here("out", d, "Sim_VulACF_budget_stats.csv"))[-1] %>% 
    mutate(Population = d)
}))


stats0 <- stats %>% 
  filter(Coverage == 0) %>% 
  select(Key, Type, Population, Inc0 = IncR, Mor0 = MorR)



stats <- stats %>% 
  mutate(
    N_Footfall = pmax(ACF_Plain_Footfall, ACF_Vul_Footfall),
    N_Screened = pmax(ACF_Plain_Screened, ACF_Vul_Screened),
    N_Confirmed = pmax(ACF_Plain_Confirmed, ACF_Vul_Confirmed),
    N_Vul = pmax(ACF_Plain_Vul, ACF_Vul_Vul),
    N_Sym = pmax(ACF_Plain_Sym, ACF_Vul_Sym),
    N_CXR = pmax(ACF_Plain_CXR, ACF_Vul_CXR),
    N_Xpert = pmax(ACF_Plain_Xpert, ACF_Vul_Xpert),
    N_Fl_DS = pmax(ACF_Plain_DS_Fl, ACF_Vul_DS_Fl),
    N_Fl_DR = pmax(ACF_Plain_DR_Fl, ACF_Vul_DR_Fl),
    N_Sl_DR = pmax(ACF_Plain_DR_Sl, ACF_Vul_DR_Sl),
    C_Screened = N_Sym * cost$Sym + N_Vul * cost$Vul + N_CXR * cost$CXR,
    C_Confirmed = ACF_Vul_Xpert * cost$Xpert + ACF_Plain_Xpert * cost_plain$Xpert,
    C_Tx = (N_Fl_DS + N_Fl_DR) * cost$Tx_Fl + (N_Fl_DR + N_Sl_DR) * cost$Tx_Sl,
    C_Total = C_Screened + C_Confirmed + C_Tx,
  ) %>% 
  select(Key, Type, Population, Coverage, Pop0, Inc1 = IncR, Mor1 = MorR,
         N_Footfall, N_Screened, N_Confirmed, starts_with("C_")) %>% 
  left_join(stats0) %>% 
  mutate(
    AvtInc = (Inc0 - Inc1) / Inc0,
    AvtMor = (Mor0 - Mor1) / Mor0,
    across(starts_with("N_"), function(x) x / Pop0 * pop.size),
    across(starts_with("C_"), function(x) x / Pop0 * pop.size)
  ) %>% 
  mutate(
    Gp = case_when(
      Type == "VulACF" & Population == "dy_lo" ~ "Vul_lo",
      Type == "VulACF" & Population == "dy_hi" ~ "Vul_hi",
      Population == "dy_hi" ~ "cxr",
      T ~ "None"
    )
  ) %>% 
  filter(Gp != "None")



g_avt1 <- stats %>% 
  group_by(Gp, Coverage, Population, Type) %>% 
  summarise(
    M = mean(AvtInc),
    L = quantile(AvtInc, 0.025),
    U = quantile(AvtInc, 0.975)
  ) %>% 
  ungroup() %>% 
  ggplot() + 
  geom_ribbon(aes(x = Coverage, ymin = L, ymax = U, fill = Gp), alpha = 0.1) +
  geom_line(aes(x = Coverage, y = M, colour = Gp)) +
  scale_y_continuous("Averted cases, %", labels = scales::percent) + 
  scale_x_continuous("Annual ACF screened, percentage population", 
                     labels = scales::percent) +
  scale_color_discrete("Scenario", labels=c(cxr="Untargeted screening",
                                            Vul_hi="Vulnerability-led, low threshold",
                                            Vul_lo="Vulnerability-led, high threshold"
                                            )) +
  guides(fill = guide_none()) +
  theme(legend.position = "right")



g_vul_imp <- stats %>% 
  group_by(Gp, Coverage, Population, Type) %>% 
  summarise(
    M = mean(AvtInc),
    L = quantile(AvtInc, 0.025),
    U = quantile(AvtInc, 0.975)
  ) %>% 
  ungroup() %>% 
  ggplot() + 
  geom_line(aes(x = Coverage, y = M, colour = Gp)) +
  scale_y_continuous("Averted cases, %", labels = scales::percent) + 
  scale_x_continuous("Annual ACF screened, percentage population", 
                     labels = scales::percent) +
  scale_color_discrete("Scenario", labels=c(cxr="Untargeted screening",
                                            Vul_hi="Vulnerability-led, low threshold",
                                            Vul_lo="Vulnerability-led, high threshold"
  )) +
  guides(fill = guide_none()) +
  theme(legend.position = c(1, 0), legend.justification = c(1 + 0.05, -0.05))



g_vul_ci <- stats %>% 
  group_by(Gp, Coverage, Population, Type) %>% 
  summarise(across(c(AvtInc, C_Total), mean)) %>% 
  ungroup() %>% 
  ggplot() + 
  geom_line(aes(x = C_Total, y = AvtInc, colour = Gp)) +
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", breaks=seq(0, 40, 10) * 1e6, labels = scales::number_format(scale = 1e-6)) + 
  # scale_y_continuous("Incident case averted, %, 2023-2030", labels = scales::percent) +
  scale_y_continuous("Averted cases, %", labels = scales::percent) + 
  scale_color_discrete("Scenario", labels=c(cxr="Untargeted screening",
                                            Vul_hi="Vulnerability-led, low threshold",
                                            Vul_lo="Vulnerability-led, high threshold"
  )) +
  expand_limits(x = 0, y = 0) +
  theme(legend.position = c(1, 0), legend.justification = c(1 + 0.05, -0.05))



g_avt1
g_vul_imp
g_vul_ci


ggsave(g_avt, filename = here::here("docs", "g_avt.png"), width = 6, height = 5)
ggsave(g_vul_imp, filename = here::here("docs", "g_vul_imp.png"), width = 6, height = 5)
ggsave(g_vul_ci, filename = here::here("docs", "g_vul_ci.png"), width = 6, height = 5)
# ggsave(g_yield, filename = here::here("docs", "g_yield.png"), width = 6, height = 3)

