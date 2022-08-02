library(tidyverse)



dirs <- dir("out")
dirs <- dirs[startsWith(dirs, "dy_")]


dt <- 0.5

summ <- bind_rows(lapply(dirs, function(d) {
  sims <- read.csv(here::here("out", d, "Runs_Intv.csv")) %>% 
    mutate(Year = Time + 0.5) %>% 
    select(Year, Pop, IncR, Scenario, Key) %>% 
    mutate(
      Inc = Pop * IncR
    ) %>% 
    group_by(Scenario, Key) %>% 
    summarise(Inc = sum(Inc) * dt)
  
  sims  %>% 
    filter(Scenario != "Baseline") %>% 
    left_join(sims %>% filter(Scenario == "Baseline") %>% 
                ungroup() %>% 
                select(Inc0 = Inc, Key)) %>% 
    mutate(
      Avt = 1 - Inc / Inc0
    ) %>% 
    group_by(Scenario) %>% 
    summarise(
      M = median(Avt),
      L = quantile(Avt, 0.25),
      U = quantile(Avt, 0.75)
    ) %>% 
    ungroup() %>% 
    mutate(Baseline = d)
}))



g <- summ %>% 
  extract(Scenario, c("ACF", "Coverage"), "(High|Mod), Cov=(\\d+\\%)") %>%  
  extract(Baseline, c("PrComorb", "OR_prev_TB"), "dy_(\\d+\\%)_(\\S+)") %>% 
  mutate(
    PrComorb = factor(PrComorb, c("5%", "10%", "20%")),
    Coverage = factor(Coverage, c("20%", "50%", "100%"))
  ) %>% 
  ggplot() +
  geom_pointrange(aes(x = PrComorb, y=M, ymin=L, ymax=U, colour=OR_prev_TB)) +
  scale_x_discrete("People in the high risk group, %") +
  scale_y_continuous("Averted incident cases during 2020-2030, %", labels = scales::percent) +
  scale_color_discrete("Odds Ratio of TB prevalence, high risk // low risk") +
  expand_limits(y = 0) +
  facet_grid(ACF~Coverage, labeller = labeller(.cols=label_both, .rows=label_both)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  labs(caption = "High risk group: rates of primary progression/reactivation/relapse")

g

ggsave(g, filename = "docs/avt_inc.pdf", width=9, height=6)




summ_epi <- bind_rows(lapply(c("dy_5%_2", "dy_20%_5"), function(d) {
  sims <- read.csv(here::here("out", d, "Runs_Intv.csv")) %>% 
    mutate(Year = Time + 0.5) %>% 
    select(Year, IncR, IncR_RiskLo, IncR_RiskHi, 
           Prev, Prev_RiskLo, Prev_RiskHi, Scenario, Key) %>% 
    pivot_longer(-c(Year, Scenario, Key)) %>% 
    extract(name, c("Index", "Group"), "(IncR|Prev)(\\w*)") %>% 
    mutate(Group = case_when(
      Group == "_RiskLo" ~ "Low Risk",
      Group == "_RiskHi" ~ "High Risk",
      T ~ "Overall"
    )) %>% 
    group_by(Year, Scenario, Index, Group) %>% 
    summarise(
      M = median(value),
      L = quantile(value, 0.25),
      U = quantile(value, 0.75)
    ) %>% 
    ungroup() %>% 
    mutate(Baseline = d)
  sims
}))


summ_epi %>% 
  extract(Scenario, c("ACF", "Coverage"), "(High|Mod), Cov=(\\d+\\%)", remove = F) %>%  
  extract(Baseline, c("PrComorb", "OR_prev_TB"), "dy_(\\d+\\%)_(\\S+)", remove = F) %>% 
  mutate(
    ACF = ifelse(is.na(ACF), "None", ACF),
    Coverage = ifelse(is.na(Coverage), "0%", Coverage),
    PrComorb = factor(PrComorb, c("5%", "10%", "20%")),
    Coverage = factor(Coverage, c("0%", "20%", "50%", "100%")),
  ) %>% 
  filter(
    (Coverage %in% c("0%", "50%")) & 
      ((PrComorb == "5%" & OR_prev_TB == "2") | (PrComorb == "20%" & OR_prev_TB == "5"))
  ) %>% 
  filter(Group == "Overall") %>%
  ggplot(aes(x = Year)) +
  geom_ribbon(aes(ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = M, colour = Scenario)) +
  scale_y_continuous("Incidence, per 100 000", labels = scales::number_format(scale = 1e5)) +
  facet_grid(Index~PrComorb + OR_prev_TB, labeller = labeller(.cols=label_both)) +
  theme_bw() +
  theme(legend.position = "bottom")
  # ggplot() +
  # geom_line(aes(x = Year, y = M, colour = Scenario)) +
  # scale_y_continuous("Per 100 000", labels = scales::number_format(scale = 1e5)) +
  # expand_limits(y = 0) +
  # facet_grid(OR_prev_TB~Index+Group, scales = "free_y")



