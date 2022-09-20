library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


cost <- read_csv(here::here("data", "cost.csv"))
cost <- as.list(setNames(cost$M, cost$Item))


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
    N_Fl_DS = pmax(ACF_Plain_DS_Fl, ACF_Vul_DS_Fl),
    N_Fl_DR = pmax(ACF_Plain_DR_Fl, ACF_Vul_DR_Fl),
    N_Sl_DR = pmax(ACF_Plain_DR_Sl, ACF_Vul_DR_Sl),
    C_Screened = N_Screened * ifelse(Type == "VulACF", cost$Vul + cost$CXR, cost$CXR),
    C_Confirmed = N_Confirmed * cost$Xpert,
    C_Tx = (N_Fl_DS + N_Fl_DR) * cost$Tx_Fl + (N_Fl_DR + N_Sl_DR) * cost$Tx_Sl,
    C_Total = C_Screened + C_Confirmed + C_Tx,
  ) %>% 
  select(Key, Type, Population, Coverage, Pop0, Inc1 = IncR, Mor1 = MorR,
         N_Footfall, N_Screened, N_Confirmed, starts_with("C_")) %>% 
  left_join(stats0) %>% 
  mutate(
    AvtInc = (Inc0 - Inc1) / Inc0,
    AvtMor = (Mor0 - Mor1) / Mor0,
    across(starts_with("N_"), function(x) x / Pop0),
    across(starts_with("C_"), function(x) x / Pop0)
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



g_avt <- stats %>% 
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
  scale_x_continuous("People screened, per year-person", 
                     labels = scales::percent) +
  scale_color_discrete("Scenario", labels=c(cxr="Universal screening",
                                            Vul_hi="Vulnerability-led, high comorbidity",
                                            Vul_lo="Vulnerability-led, low comorbidity"
                                            )) +
  guides(fill = guide_none()) +
  theme(legend.position = "right")



g_avt


# 
# stats %>% 
#   mutate(kg = paste0(Gp, Key)) %>% 
#   ggplot() +
#   geom_line(aes(x = C_Total, y = AvtInc, colour = Gp, group = kg)) +
#   scale_color_discrete("Scenario", labels=c(cxr="Universal screening",
#                                             Vul_hi="Vulnerability-led, high comorbidity",
#                                             Vul_lo="Vulnerability-led, low comorbidity"
#   )) +
#   scale_y_continuous("Averted cases, %", labels = scales::percent) + 
#   scale_x_continuous("Total ACF cost, per year-person", 
#                      labels = scales::percent)
# 
# 
# stats %>% 
#   mutate(kg = paste0(Gp, Key)) %>% 
#   ggplot() +
#   geom_line(aes(x = N_Footfall, y = AvtInc, colour = Gp, group = kg)) +
#   scale_color_discrete("Scenario", labels=c(cxr="Universal screening",
#                                             Vul_hi="Vulnerability-led, high comorbidity",
#                                             Vul_lo="Vulnerability-led, low comorbidity"
#   )) +
#   scale_y_continuous("Averted cases, %", labels = scales::percent) + 
#   scale_x_continuous("People reached, per year-person", 
#                      labels = scales::percent) +




# g_yield <- sims %>% 
#   mutate(
#     p_hi = sprintf("Pr(High risk)=%s", scales::percent(P_Comorb, accuracy=.1)),
#     OR_ComorbTB = factor(OR_ComorbTB),
#     yield = - (N_Reached_M * 17.53) / dAvt_M 
#   ) %>% 
#   ggplot() +
#   geom_ribbon(aes(x = yield, ymin = N_Reached_L, ymax = N_Reached_U), alpha = 0.1) +
#   geom_line(aes(x = yield, y = N_Reached_M)) +
#   scale_x_continuous("Cost per case averted, $") + 
#   scale_y_continuous("Number of people screened per year", 
#                      labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
#   facet_grid(.~p_hi, scales = 'free_x') +
#   expand_limits(y = 0) +
#   labs(caption = 'Population size: 1 million')
# 
# 
# g_yield


ggsave(g_avt, filename = here::here("docs", "g_avt.png"), width = 6, height = 5)
# ggsave(g_yield, filename = here::here("docs", "g_yield.png"), width = 6, height = 3)

