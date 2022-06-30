library(tidyverse)
library(jsonlite)


theme_set(theme_bw())



targets <- jsonlite::read_json(here::here("data", "targets.json"))
targets <- bind_rows(targets$All) %>% filter(Index %in% c("Prev", "ARTI", "PrSym"))


pars <- read_csv(here::here("out", "dy", "Post.csv")) %>% 
  rename(Key = `...1`)

sims <- read_csv(here::here("out", "dy", "Runs_Intv.csv")) %>% 
  mutate(Year = Time)



sims_epi <- sims %>% 
  select(Year, Scenario, IncR, MorR) %>% 
  pivot_longer(-c(Year, Scenario), names_to = "Index") %>% 
  group_by(Year, Scenario, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)  
  ) %>% 
  extract(Scenario, c("Acc", "Tested", "Type"), "(\\w+), Budget=(\\S+), (Risk|Universal)", 
          remove = F, convert = F)


sims_epi %>% 
  ggplot(aes(x = Year)) +
#  geom_ribbon(aes(ymin = L, ymax = U, fill = Scenario), alpha = 0.3) +
  geom_line(aes(y = M, colour = Tested)) +
  scale_y_continuous("per 100k", labels = scales::number_format(scale = 1e5)) +
  facet_grid(Index~Type+Acc, scales = "free_y")
  


sims %>% 
  select(Year, Scenario, Key, IncR, MorR) %>% 
  pivot_longer(-c(Year, Scenario, Key), names_to = "Index") %>% 
  left_join(pars %>% select(Key, PrHi = p_comorb, RR = rr_risk_comorb)) %>% 
  mutate(PrHi = cut(PrHi, c(0, 0.1, 0.25, 0.5))) %>% 
  group_by(Year, Scenario, Index) %>% 
  summarise(
    M = mean(value)
  ) %>% 
  extract(Scenario, c("Acc", "Tested", "Type"), "(\\w+), Budget=(\\S+), (Risk|Universal)", 
          remove = F, convert = T) %>% 
  filter(Scenario != "Baseline" & Acc == "High") %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = M, colour = Tested, group = Tested)) +
  scale_y_continuous("per 100k", labels = scales::number_format(scale = 1e5)) +
  scale_colour_viridis_c() +
  facet_grid(Index~Type+Acc, scales = "free_y")




temp <- sims %>%
  mutate(Inc = IncR * Pop, Mor = MorR * Pop) %>% 
  select(Year, Key, Scenario, Inc, Mor)

g_avt <- temp %>%
  filter(Scenario != "Baseline") %>% 
  left_join(temp %>% filter(Scenario == "Baseline") %>% select(Key, Year, Inc0 = Inc, Mor0 = Mor)) %>% 
  left_join(pars %>% select(Key, PrHi = p_comorb)) %>% 
  mutate(
    PrHi = cut(PrHi, c(0, 0.1, 0.2, 0.5)),
    AvtInc = 1 - Inc / Inc0,
    AvtMor = 1 - Mor / Mor0
  ) %>% 
  group_by(Scenario, Year, PrHi) %>% 
  summarise(across(starts_with("Avt"), mean)) %>% 
  extract(Scenario, c("Acc", "Tested", "Type"), "(\\w+), Budget=(\\S+), (Risk|Universal)", 
          remove = F, convert = T) %>% 
  filter(Acc == "High") %>% 
  ggplot() +
  geom_line(aes(x = Year, y = AvtInc, colour = Tested, group = Tested)) +
  scale_y_continuous("Incident cases averted, %", labels = scales::percent) +
  scale_x_continuous("Year", breaks = c(2020, 2025, 2030)) +
  scale_colour_viridis_c("People reached by the ACF, % population", breaks = c(0.01, 0.03, 0.05), labels = scales::percent) +
  facet_grid(PrHi~Type + Acc, labeller = labeller(PrHi=c("(0,0.1]"="0-10%", "(0.1,0.2]"="10-20%", "(0.2,0.5]"="20-50%"),
                                                  Type=c("Risk"="Risk-based ACF", "Universal"="Community-wise ACF"))) +
  theme(legend.position = "bottom")



# sims_acc <- sims %>% 
#   group_by(Key, Scenario) %>% 
#   summarise(
#     AccInc = sum(Pop * IncR) * 0.5,
#     AccPrev = sum(Pop * Prev) * 0.5,
#     AccMor = sum(Pop * MorR) * 0.5
#   )
# 
#   
# 
# 
# sims_acc %>% write_csv(here::here("out", "dy", "Sims_Acc.csv"))


sims_acc <- read_csv(here::here("out", "dy", "Sims_Acc.csv"))



g_eff <- sims_acc %>% 
  filter(Scenario != "Baseline") %>% 
  left_join(sims_acc %>% filter(Scenario == "Baseline") %>% 
              select(Key, AccInc0 = AccInc, AccPrev0 = AccPrev, AccMor0 = AccMor)) %>% 
  left_join(pars %>% select(Key, PrHi = p_comorb, RR = rr_risk_comorb)) %>% 
  mutate(
    PrHi = cut(PrHi, c(0, 0.1, 0.2, 0.5)),
    AvtInc = 1 - AccInc / AccInc0, 
    AvtPrev = 1 - AccPrev / AccPrev0, 
    AvtMor = 1 - AccMor / AccMor0
  ) %>% 
  extract(Scenario, c("Acc", "Tested", "Type"), "(\\w+), Budget=(\\S+), (Risk|Universal)", convert = T) %>% 
  select(Key, Acc, Tested, Type, Eff = AvtInc, PrHi, RR) %>% 
  pivot_wider(names_from = Type, values_from = Eff) %>% 
  group_by(Acc, Tested, PrHi) %>% 
  summarise(Risk = mean(Risk > Universal)) %>% 
  ggplot() +
  geom_line(aes(x = Tested, y = Risk, colour = PrHi)) +
  scale_y_continuous("Prob. Risk-based ACF prevent more incident cases, %", labels = scales::percent) +
  scale_x_continuous("People reached by the ACF, percentage population", labels = scales::percent) + 
  scale_colour_discrete("High risk group share", 
                        labels = c("(0,0.1]"="0-10%", "(0.1,0.2]"="10-20%", "(0.2,0.5]"="20-50%"),
                        guide = guide_legend(reverse = TRUE)) + 
  facet_grid(Acc~.)



ggsave(g_eff, filename = here::here("docs", "g_eff.png"), width = 7, height = 8)
ggsave(g_avt, filename = here::here("docs", "g_avt.png"), width = 7, height = 8)




