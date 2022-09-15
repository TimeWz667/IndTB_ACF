library(tidyverse)
library(jsonlite)


theme_set(theme_bw())



targets <- jsonlite::read_json(here::here("data", "targets.json"))
targets <- bind_rows(targets$All) %>% filter(Index %in% c("Prev", "ARTI", "PrSym"))


pars <- read_csv(here::here("out", "dy_free", "Post.csv")) %>% 
  rename(Key = `...1`)

sims <- read_csv(here::here("out", "dy_free", "Runs_IntvE.csv")) %>% 
  mutate(Year = Time)


sims_acc <- sims %>%
  group_by(Key, Scenario) %>%
  summarise(
    AccInc = sum(Pop * IncR) * 0.5,
    AccPrev = sum(Pop * Prev) * 0.5,
    AccMor = sum(Pop * MorR) * 0.5,
    N_ACF_Reached = sum(N_ACF_Reached) / Pop[1] * 0.5
  )


save(sims_acc, file=here::here("out", "dy_free", "Sims_Acc.rdata") )




sims_acc %>% 
  filter(Scenario != "Baseline") %>% 
  left_join(sims_acc %>% filter(Scenario == "Baseline") %>% 
              select(Key, AccInc0 = AccInc, AccPrev0 = AccPrev, AccMor0 = AccMor)) %>% 
  left_join(pars %>% select(Key, PrHi = p_comorb, RR = rr_risk_comorb)) %>% 
  mutate(
    PrHi = cut(PrHi, c(0, 0.1, 0.2, 0.4, 0.5)),
    AvtInc = 1 - AccInc / AccInc0, 
    AvtPrev = 1 - AccPrev / AccPrev0, 
    AvtMor = 1 - AccMor / AccMor0,
    Cost = N_ACF_Reached * 19.52
  ) %>% 
  extract(Scenario, c("Acc", "Tested", "Type"), "(\\w+), Budget=(\\S+), (Risk|Universal)", convert = T) %>% 
  select(Key, Acc, Tested, Type, Eff = AvtInc, Cost, PrHi, RR) %>%
  filter(Acc == 'High') %>% 
  mutate(
    Gp = case_when(
      (PrHi == '(0,0.1]') & (Type == 'Risk') ~ 'In vulnerable only, 0.7%',
      (PrHi == '(0.1,0.2]') & (Type == 'Risk') ~ 'In vulnerable only, 17.8%',
      Type == 'Universal' ~ 'Untargeted',
      T ~ 'Na'
    )
  ) %>% 
  filter(Gp != 'Na') %>% 
  group_by(Gp, Tested) %>% 
  summarise(Eff = mean(Eff), Cost = mean(Cost) * 2) %>% 
  # filter(Cost > 20) %>% 
  ggplot() +
  geom_line(aes(x = Cost, y = Eff, colour = Gp)) +
  # geom_point(aes(x = Cost, y = Eff, colour = Gp)) +
  scale_x_continuous("Incremental cost, 2023-2030, USD million") +
  # geom_line(aes(x = Tested, y = Eff, colour = Gp)) +
  # scale_x_continuous("Proportional population reached by ACF, per year", labels = scales::percent) +
  scale_y_continuous("Percentage cases averted, 2023-2030", labels = scales::percent) + 
  scale_color_discrete("ACF strategy")
                            

