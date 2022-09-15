library(tidyverse)


theme_set(theme_bw())


ds <- dir("out")
ds <- ds[startsWith(ds, "dy_sc")]
ds <- setNames(ds, ds)
ds <- ds[ds != "dy_sc2-1"]


scs <- bind_rows(lapply(ds, function(d) {
  ys <- read_csv(here::here("out", d, "Runs_Post.csv"))
  ys %>% 
    mutate(File = d, OR_ComorbTB = or_comorb, P_Comorb = pr_tb) %>% 
    select(File, OR_ComorbTB, P_Comorb) %>% 
    distinct()
}))


scs


ppv0 <- 0.8


sims <- bind_rows(lapply(ds, function(d) {
  post <- read_csv(here::here("out", d, "Runs_ACF_t.csv"))
  
  
  p_nontb <- post %>% 
    mutate(Year = Time - 0.5) %>% 
    filter(Year == 2020) %>% 
    extract(Scenario, "ACF", "(High|Mod|Baseline)", remove = F) %>% 
    mutate(
      Spec = ifelse(ACF == "High", 0.99, 0.98),
      N_ACF_FP_Detected = N_ACF_Detected * (1 / ppv0 - 1),
      N_ACF_NonTB_Reached = N_ACF_FP_Detected / (1 - Spec),
      P_ACF_NonTB_Reached = N_ACF_NonTB_Reached / Pop
    ) %>% 
    select(Key, Scenario, P_ACF_NonTB_Reached)
  
  
  incs = post %>% 
    group_by(Key, Scenario) %>% 
    left_join(p_nontb) %>% 
    summarise(
      Pop = sum(Pop) * 0.5,
      Pop_RiskHi = sum(Pop_RiskHi) * 0.5,
      AccInc = sum(Pop * IncR) * 0.5,
      TP_Reached = sum(N_ACF_Reached) * 0.5,
      FP_Reached = sum(Pop * P_ACF_NonTB_Reached) * 0.5,
      ACF_Reached = TP_Reached + FP_Reached,
      PrTB_ACF = TP_Reached / max(ACF_Reached, 1e-10)
    ) %>% 
    ungroup()
  
  
  incs0 <- incs %>% 
    filter(Scenario == "Baseline") %>% 
    select(Key, AccInc0 = AccInc)
  
  
  impacts <- incs %>% 
    filter(Scenario != "Baseline") %>% 
    rename(AccInc1 = AccInc) %>% 
    left_join(incs0) %>% 
    mutate(
      PrAvt = (1 - AccInc1 / AccInc0)
    ) %>% 
    extract(Scenario, c("ACF", "Coverage", "r_acf"), "(Mod|High), (Focus|Universal), R_ACF=(\\S+)") %>% 
    mutate(
      PrPop = ACF_Reached / Pop,
      PrTar = ifelse(Coverage == "Focus", ACF_Reached / Pop_RiskHi, PrPop)
    ) %>% 
    group_by(r_acf, ACF, Coverage) %>% 
    summarise(
      across(c(PrPop, PrTar, PrAvt), list(
        M = median,
        L = function(x) quantile(x, 0.25),
        U = function(x) quantile(x, 0.75)
      ))
    )
  
  impacts %>% mutate(File = d)
})) %>% 
  left_join(scs)


g_per_tar <- sims %>% 
  mutate(
    p_hi = sprintf("Pr(High risk)=%s", scales::percent(P_Comorb, accuracy=.1)),
    OR_ComorbTB = factor(OR_ComorbTB)
  ) %>% 
  ggplot() +
  geom_line(aes(x = PrTar_M, y = PrAvt_M, colour = OR_ComorbTB)) + 
  geom_point(aes(x = PrTar_M, y = PrAvt_M, colour = OR_ComorbTB)) + 
  scale_y_continuous("Incident cases averted, %", labels = scales::percent) + 
  scale_x_continuous("ACF coverage per targeted population-year", labels = scales::percent) + 
  scale_colour_discrete("Odds ratio of TB") +
  facet_wrap(p_hi~., nrow=2) +
  expand_limits(x = 0, y = 0)



g_per_cap <- sims %>% 
  mutate(
    p_hi = sprintf("Pr(High risk)=%s", scales::percent(P_Comorb, accuracy=.1)),
    OR_ComorbTB = factor(OR_ComorbTB)
  ) %>% 
  ggplot() +
  geom_line(aes(x = PrPop_M, y = PrAvt_M, colour = OR_ComorbTB)) + 
  geom_point(aes(x = PrPop_M, y = PrAvt_M, colour = OR_ComorbTB)) + 
  scale_y_continuous("Incident cases averted, %", labels = scales::percent) + 
  scale_x_continuous("ACF coverage per targeted population-year", labels = scales::percent) + 
  scale_colour_discrete("Odds ratio of TB") +
  facet_wrap(p_hi~., nrow=2) +
  expand_limits(x = 0, y = 0)



ggsave(g_per_tar, filename = here::here("docs", "scale", "g_per_tar.png"), width = 6, height = 7)
ggsave(g_per_cap, filename = here::here("docs", "scale", "g_per_cap.png"), width = 6, height = 7)


