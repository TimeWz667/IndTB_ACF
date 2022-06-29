library(tidyverse)
library(jsonlite)


theme_set(theme_bw())



targets <- jsonlite::read_json(here::here("data", "targets.json"))
targets <- bind_rows(targets$All) %>% filter(Index %in% c("Prev", "ARTI", "PrSym"))


pars <- read_csv(here::here("out", "dy", "Post.csv")) %>% 
  rename(Key = `...1`)

sims <- read_csv(here::here("out", "dy", "Runs_Post.csv")) %>% 
  mutate(Year = Time + 0.5)


sims %>% 
  filter(Year > 2000) %>% 
  select(Year, MorR, IncR, Prev, ARTI, PrSym, PrDR_Inc) %>% 
  pivot_longer(-Year, names_to = "Index") %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = median(value), 
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  ggplot() + 
  geom_ribbon(aes(x = Year, ymin = L, ymax = U, fill = "Simulation"), alpha = 0.4) +
  geom_line(aes(x = Year, y = M, colour = "Simulation")) + 
  geom_pointrange(data = targets, aes(x = Year, y = M, ymin = L, ymax = U)) + 
  facet_wrap(.~Index, scales = "free_y") +
  expand_limits(y = 0) +
  guides(fill = guide_none())



linked <- sims %>% 
  filter(Year == 2020) %>% 
  select(Key, Prev) %>% 
  left_join(pars %>% select(Key, p_comorb, rr_risk_comorb))


linked %>% 
  pivot_longer(-c(Key, Prev), names_to = "Parameter") %>% 
  ggplot() +
  geom_point(aes(x = value, y = Prev)) +
  scale_y_continuous("Prevalence, per 100k", labels = scales::number_format(scale = 1e5)) +
  scale_x_continuous("") +
  facet_grid(.~Parameter, scales = "free_x")


linked %>% 
  ggplot() +
  geom_point(aes(x = p_comorb, y = rr_risk_comorb))



