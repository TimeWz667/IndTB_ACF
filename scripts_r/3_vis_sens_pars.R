library(tidyverse)
library(latex2exp)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))



sens <- local({
  sens <- read_csv(here::here("out", "sens", "Sens_Pars_stats.csv"))[, -1]
  
  sens0 <- sens %>% 
    filter(Change == 0) %>% 
    select(Inc0=Inc, Mor0=Mor, Par, Key, Alg, Coverage)
  
  sens1 <- sens %>% 
    select(Inc1=Inc, Mor1=Mor, Par, Change, Key, Alg, Coverage)
  
  sens1 %>% left_join(sens0)
})


stats <- sens %>% 
  filter(abs(Change) == 0.05) %>% 
  group_by(Par, Change, Alg, Coverage) %>% 
  summarise(across(c(Inc1, Inc0, Mor1, Mor0), mean))


ord_inc <- stats %>% group_by(Par) %>% 
  summarise(x = diff(range(Inc1))) %>% 
  mutate(Par = reorder(Par, x)) %>% 
  arrange(Par) %>% pull(Par) %>% as.character()

ord_mor <- stats %>% group_by(Par) %>% 
  summarise(x = diff(range(Mor1))) %>% 
  mutate(Par = reorder(Par, x)) %>% 
  arrange(Par) %>% pull(Par) %>% as.character()


lab_pars <- c(
  "rr_beta_dr" = r"($\rho_{dr}$)", 
  "r_relapse_td" = r"($\rho_{dr}$)", 
  "rr_sus_ltbi"= r"($\pi$)",  
  "r_relapse_st" = r"($\gamma_{st}$)",
  "r_clear" = r"($\nu_{c}$)",      
  "r_cs_s" = r"($\delta_{0}$)", 
  "r_relapse_tc"= r"($\gamma_{tc}$)", 
  "r_onset"= r"($\sigma$)", 
  "r_react"= r"($\gamma_{re}$)",
  "p_entry_pri" = r"($p_{Pr0}$)", 
  "beta_ds" = r"($\beta_{0}$)",
  "rr_inf_asym" = r"($\rho_{A}$)",
  "p_dx_pri" = r"($p_{dx, Pub}$)",
  "r_cs_c" = r"($\delta_{1}$)",
  "p_dx_pub" = r"($p_{dx, Pub}$)"
)


g_sens_par_inc <- stats %>% 
  mutate(
    Par = factor(Par, ord_inc),
    i = as.numeric(Par),
    col = ifelse(sign(Change) > 0, "+5%", "-5%")
  ) %>% 
  ggplot() + 
  #geom_rect(aes(ymax = i + 0.45, ymin = i - 0.45, xmin = Inc1, xmax = Inc0, fill = col)) +
  geom_bar(aes(y = Par, x = Inc1 - Inc0 , fill = col), stat = 'Identity') +
  scale_x_continuous("Change in averted incidence, difference in %", 
                     labels = scales::percent, 
                     sec.axis = sec_axis(~.+as.numeric(stats[1, "Inc0"]), 
                                         labels = scales::percent, 
                                         name = "Averted incidence, %")) +
  scale_y_discrete("Parameter", labels = TeX(lab_pars)) +
  scale_fill_discrete("Change in parameter value") +
  expand_limits(x = c(-0.0018, 0.0018)) +
  theme(legend.position = c(1, 0), legend.just = c(1.05, 0))  


g_sens_par_mor <- stats %>% 
  mutate(
    Par = factor(Par, ord_mor),
    i = as.numeric(Par),
    col = ifelse(sign(Change) > 0, "+5%", "-5%")
  ) %>% 
  ggplot() + 
  #geom_rect(aes(ymax = i + 0.45, ymin = i - 0.45, xmin = Inc1, xmax = Inc0, fill = col)) +
  geom_bar(aes(y = Par, x = Mor1 - Mor0 , fill = col), stat = 'Identity') +
  scale_x_continuous("Change in averted mortality, difference in %", 
                     labels = scales::percent, 
                     sec.axis = sec_axis(~.+as.numeric(stats[1, "Mor0"]), 
                                         labels = scales::percent, 
                                         name = "Averted mortality, %")) +
  scale_y_discrete("Parameter", labels = TeX(lab_pars)) +
  scale_fill_discrete("Change in parameter value") +
  expand_limits(x = c(-0.0025, 0.0025)) +
  theme(legend.position = c(1, 0), legend.just = c(1.05, 0))  


g_sens_par_mor
g_sens_par_inc


g_sens_par <- ggpubr::ggarrange(g_sens_par_inc, g_sens_par_mor, ncol = 1)

g_sens_par

ggsave(g_sens_par, filename = here::here("docs", "figs", "g_sens_par.png"), width = 6, height = 7)

