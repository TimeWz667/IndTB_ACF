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
  ggplot() +
  geom_pointrange(aes(x = PrComorb, y=M, ymin=L, ymax=U, colour=ACF)) +
  scale_x_discrete("People in the high risk group, %") +
  scale_y_continuous("Averted incident cases during 2020-2030, %", labels = scales::percent) +
  scale_color_discrete("RR of being screened, high risk // low risk") +
  expand_limits(y = 0) +
  facet_grid(OR_prev_TB~Coverage, labeller = labeller(.cols=label_both)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  labs(caption = "High risk group: higher susceptibility and rates of primary progression/reactivation/relapse")

g

ggsave(g, filename = "docs/avt_inc.pdf", width=9, height=6)






