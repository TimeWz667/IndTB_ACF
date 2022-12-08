library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


labs_display <- c(
  ut="MDU screening",
  utsym="Symptom screening",
  vul="Vulnerability-led MDU"
) 


cost <- read_csv(here::here("data", "cost.csv"))
cost <- as.list(setNames(cost$Vul, cost$Item))


ds <- "main" #, "dy_lo")

stats <- local({
  stats <- read_csv(here::here("out", ds, "Sim_VulACF_budget_stats.csv"))[-1]
  
  
  stats0 <- stats %>% 
    filter(Coverage == 0) %>% 
    select(Key, Alg, Inc0 = Inc, Mor0 = Mor)
  
  
  stats <- stats %>% 
    mutate(
      across(starts_with("ACF_"), function(x) ifelse(is.na(x), 0, x)),
      N_Footfall = ACF_Footfall,
      N_Screened = ACF_Screened,
      N_Confirmed = ACF_Uti_xpert,
      N_Vul = ACF_Uti_vul,
      N_Sym = ACF_Uti_sym,
      N_VulSym = ACF_Uti_vs,
      N_CXR = ACF_Uti_cxr,
      N_Xpert = ACF_Uti_xpert,
      N_Fl = ACF_Fl,
      N_Sl = ACF_Sl,
      N_TPT = ACF_TPT,
      C_Screened = N_Sym * cost$Sym + N_Vul * cost$Vul + N_VulSym * cost$VSC + N_CXR * cost$CXR,
      C_Confirmed = N_Xpert * cost$Xpert,
      C_Tx = (N_Fl + N_TPT) * cost$Tx_Fl + N_Sl * cost$Tx_Sl,
      C_Total = C_Screened + C_Confirmed + C_Tx,
    ) %>% 
    select(Key, Alg, Coverage, Inc1 = Inc, Mor1 = Mor,
           N_Footfall, N_Screened, N_Confirmed, starts_with("C_")) %>% 
    left_join(stats0) %>% 
    mutate(
      AvtInc = (Inc0 - Inc1) / Inc0,
      AvtMor = (Mor0 - Mor1) / Mor0
    ) %>% 
    mutate(
      Gp = case_when(
        Alg == "VSC" ~ "vul",
        Alg == "SyCx" ~ "ut",
        Alg == "Sy" ~ "utsym",
        T ~ "None"
      )
    ) %>% 
    filter(Gp != "None")
})


g_vul_covimp <- stats %>% 
  group_by(Gp, Coverage, Alg) %>% 
  summarise(
    M = mean(AvtInc),
  ) %>% 
  ungroup() %>% 
  ggplot() + 
  # geom_point(aes(x = Coverage, y = M, colour = Gp)) +
  geom_line(aes(x = Coverage, y = M, colour = Gp)) +
  scale_y_continuous("Averted cases, %", labels = scales::percent) + 
  scale_x_continuous("Annual ACF screened, percentage population", 
                     labels = scales::percent) +
  scale_color_discrete("Scenario", labels=labs_display) +
  guides(fill = guide_none()) +
  theme(legend.position = c(0, 1), legend.justification = c(- 0.05, 1 +0.05))


g_vul_covimp


g_vul_covimp_ci <- stats %>% 
  group_by(Gp, Coverage, Alg) %>% 
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
  scale_color_discrete("Scenario", labels=labs_display) +
  guides(fill = guide_none()) +
  theme(legend.position = c(0, 1), legend.justification = c(-0.1, 1.1))


g_vul_covimp_ci



g_vul_costimp <- stats %>% 
  group_by(Gp, Coverage, Alg) %>% 
  summarise(across(c(AvtInc, C_Total), mean)) %>% 
  ungroup() %>% 
  ggplot() + 
  geom_line(aes(x = C_Total, y = AvtInc, colour = Gp)) +
  #geom_point(aes(x = C_Total, y = AvtInc, colour = Gp)) +
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", 
                     breaks=c(seq(0, 50, 25), seq(100, 1000, 100)) * 1e6, 
                     # limits = c(0, 40e6),
                     labels = scales::number_format(scale = 1e-6)) + 
  # scale_y_continuous("Incident case averted, %, 2023-2030", labels = scales::percent) +
  scale_y_continuous("Averted cases, %", labels = scales::percent
                     # limits = c(0, 0.15)
                     ) + 
  scale_color_discrete("Scenario", labels=labs_display) +
  expand_limits(x = 0, y = 0) +
  theme(legend.position = c(1, 0), legend.justification = c(1 + 0.05, -0.05))


ref <- stats %>% 
  group_by(Gp, Coverage, Alg) %>% 
  summarise(across(c(AvtInc, C_Total), mean)) %>% 
  ungroup() %>% 
  filter(Coverage %in% c(0.1, 0.2, 0.5)) %>% 
  mutate(
    Coverage = scales::percent(Coverage)
  )


g_vul_costimp_ref <- g_vul_costimp + 
  geom_point(data = ref, aes(x = C_Total, y = AvtInc, colour = Gp, shape = Coverage)) +
  scale_shape("Coverage: % Screened per year") +
  expand_limits(x = c(0, 80e6))


g_vul_costimp
g_vul_costimp_ref

c0 <- stats %>% 
  group_by(Alg, Gp) %>% 
  filter(Coverage == max(Coverage)) %>% 
  summarise(C0 = mean(C_Total))

g_vul_costimp_ci <- stats %>% 
  select(Key, Alg, Gp, Coverage, C_Total, AvtInc) %>% 
  group_by(Key, Alg, Gp) %>% 
  left_join(c0) %>% 
  summarise(
    C_Total_out = seq(0, min(C0), length.out=30), 
    AvtInc = spline(x = C_Total, y = AvtInc, xout = C_Total_out)$y
  ) %>% 
  rename(C_Total = C_Total_out) %>% 
  group_by(Alg, Gp, C_Total) %>% 
  summarise(
    M = mean(AvtInc),
    L = quantile(AvtInc, 0.025),
    U = quantile(AvtInc, 0.975)
  ) %>% 
  ungroup() %>% 
  ggplot() + 
  geom_ribbon(aes(x = C_Total, ymin = L, ymax = U, fill = Gp), alpha = 0.1) + 
  geom_point(data = ref, aes(x = C_Total, y = AvtInc, colour = Gp, shape = Coverage)) +
  geom_line(aes(x = C_Total, y = M, colour = Gp)) +
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", 
                     breaks=c(0, 10, 30, 50) * 1e6, 
                     # limits = c(0, 40e6),
                     labels = scales::number_format(scale = 1e-6)) + 
  scale_y_continuous("Averted cases, %", labels = scales::percent) + 
  scale_color_discrete("Scenario", labels=labs_display) +
  guides(fill = guide_none()) +
  theme(legend.position = c(1, 0), legend.justification = c(1.05, -0.05), legend.box = "horizontal")
  
  
g_vul_costimp_ci


g_vul_covimp
g_vul_covimp_ci
g_vul_costimp
g_vul_costimp_ci

g_vul_bind <- ggpubr::ggarrange(g_vul_covimp + theme(legend.position = "none"), g_vul_costimp, nrow = 1)


ggsave(g_vul_covimp, filename = here::here("docs", "figs", "g_vul_covimp.png"), width = 8, height = 5)
ggsave(g_vul_covimp_ci, filename = here::here("docs", "figs", "g_vul_covimp_ci.png"), width = 6, height = 5)
ggsave(g_vul_costimp, filename = here::here("docs", "figs", "g_vul_costimp.png"), width = 6, height = 5)
ggsave(g_vul_costimp_ref, filename = here::here("docs", "figs", "g_vul_costimp_ref.png"), width = 6, height = 5)
ggsave(g_vul_costimp_ci, filename = here::here("docs", "figs", "g_vul_costimp_ci.png"), width = 6, height = 5)
ggsave(g_vul_bind, filename = here::here("docs", "figs", "g_vul_bind.png"), width = 10, height = 5)
# ggsave(g_yield, filename = here::here("docs", "g_yield.png"), width = 6, height = 3)

