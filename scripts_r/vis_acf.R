library(tidyverse)



theme_set(theme_bw())


ds <- c("dy_free")

d <- ds[1]
ys <- read_csv(here::here("out", d, "Sim_VulACF_budget.csv"))[-1]


stats <- ys %>% 
  group_by(Key, Coverage, p_comorb, rr) %>%
  mutate(across(everything(), function(x) x * Pop)) %>% 
  select(-Pop) %>% 
  summarise(across(everything(), function(x) sum(x) * 0.5)) %>% 
  ungroup()



stats %>% 
  left_join(stats %>% filter(Coverage == 0) %>% select(Key, IncR0 = IncR)) %>% 
  ungroup() %>% 
  mutate(Avt = 1 - IncR / IncR0) %>% 
  ggplot() + 
  geom_line(aes(x = Coverage, y = Avt, colour = p_comorb, group = p_comorb))



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
  
  
  
  incs <- post %>% 
    group_by(Key, Scenario) %>% 
    left_join(p_nontb) %>% 
    mutate(ratio = pop0 / Pop) %>% 
    summarise(
      Range = diff(range(Time)),
      Pop_RiskHi = sum(Pop_RiskHi * ratio) * 0.5,
      AccInc = sum(Pop * IncR * ratio) * 0.5,
      TP_Reached = sum(N_ACF_Reached * ratio) * 0.5,
      FP_Reached = sum(Pop * P_ACF_NonTB_Reached * ratio) * 0.5,
      ACF_Reached = TP_Reached + FP_Reached,
      PrTB_ACF = TP_Reached / max(ACF_Reached, 1e-10),
      Pop = pop0
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
      PrAvt = (1 - AccInc1 / AccInc0),
      dAvt = AccInc1 - AccInc0
    ) %>% 
    extract(Scenario, c("ACF", "Coverage", "r_acf"), "(Mod|High), (Focus|Universal), R_ACF=(\\S+)") %>% 
    mutate(
      PrPop = ACF_Reached / Pop / Range,
      PrTar = ifelse(Coverage == "Focus", ACF_Reached / Pop_RiskHi / Range, PrPop),
      N_Reached = ACF_Reached / Range
    ) %>% 
    group_by(r_acf, ACF, Coverage) %>% 
    summarise(
      across(c(PrPop, PrTar, PrAvt, N_Reached, dAvt), list(
        M = median,
        L = function(x) quantile(x, 0.25),
        U = function(x) quantile(x, 0.75)
      ))
    )
  
  impacts %>% mutate(File = d)
})) %>% 
  left_join(scs)



g_avt <- sims %>% 
  mutate(
    p_hi = sprintf("Pr(High risk)=%s", scales::percent(P_Comorb, accuracy=.1)),
    OR_ComorbTB = factor(OR_ComorbTB)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = N_Reached_M, ymin = PrAvt_L, ymax = PrAvt_U), alpha = 0.1) +
  geom_line(aes(x = N_Reached_M, y = PrAvt_M)) +
  scale_y_continuous("Averted cases, %", labels = scales::percent) + 
  scale_x_continuous("Number of people screened per year", 
                     labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  facet_grid(.~p_hi, scales = 'free_x') +
  expand_limits(x = 0, y = 0) +
  labs(caption = 'Population size: 1 million')


g_avt



g_yield <- sims %>% 
  mutate(
    p_hi = sprintf("Pr(High risk)=%s", scales::percent(P_Comorb, accuracy=.1)),
    OR_ComorbTB = factor(OR_ComorbTB),
    yield = - (N_Reached_M * 17.53) / dAvt_M 
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = yield, ymin = N_Reached_L, ymax = N_Reached_U), alpha = 0.1) +
  geom_line(aes(x = yield, y = N_Reached_M)) +
  scale_x_continuous("Cost per case averted, $") + 
  scale_y_continuous("Number of people screened per year", 
                     labels = scales::unit_format(unit = "K", scale = 1e-3)) + 
  facet_grid(.~p_hi, scales = 'free_x') +
  expand_limits(y = 0) +
  labs(caption = 'Population size: 1 million')


g_yield


ggsave(g_avt, filename = here::here("docs", "g_avt.png"), width = 6, height = 3)
ggsave(g_yield, filename = here::here("docs", "g_yield.png"), width = 6, height = 3)

