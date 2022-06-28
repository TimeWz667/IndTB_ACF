library(tidyverse)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))




ds <- dir("out")
ds <- ds[startsWith(ds, "dy_sc")]
ds <- setNames(ds, ds)


post <- bind_rows(lapply(ds, function(d) {
  ys <- read_csv(here::here("out", d, "Runs_Post.csv"))
  ys <- ys %>% mutate(Scenario = d)
  ys
})) %>% 
  select(Time, MorR, IncR, Prev, Scenario) %>% 
  group_by(Time, Scenario) %>% 
  summarise(
    MorR_M = mean(MorR),
    MorR_L = quantile(MorR, 0.25),
    MorR_U = quantile(MorR, 0.75),
    IncR_M = mean(IncR),
    IncR_L = quantile(IncR, 0.25),
    IncR_U = quantile(IncR, 0.75),    
    Prev_M = mean(Prev),
    Prev_L = quantile(Prev, 0.25),
    Prev_U = quantile(Prev, 0.75)
  )



post %>% 
  filter(Time > 2000) %>% 
  ggplot() + 
  geom_ribbon(aes(x = Time, ymin = IncR_L, ymax = IncR_U), alpha = 0.3) +
  geom_line(aes(x = Time, y = IncR_M)) + 
  scale_y_continuous("Incidence, per 100k", labels = scales::number_format(scale = 1e5)) + 
  scale_x_continuous("Year", breaks = seq(2000, 2020, 10)) + 
  expand_limits(y = 0) + 
  facet_grid(.~Scenario)


post %>% 
  filter(Time > 2000) %>% 
  ggplot() + 
  geom_ribbon(aes(x = Time, ymin = Prev_L, ymax = Prev_U), alpha = 0.3) +
  geom_line(aes(x = Time, y = Prev_M)) + 
  scale_y_continuous("Prevalence, per 100k", labels = scales::number_format(scale = 1e5)) + 
  scale_x_continuous("Year", breaks = seq(2000, 2020, 10)) + 
  expand_limits(y = 0) + 
  facet_grid(.~Scenario)


post %>% 
  filter(Time > 2000) %>% 
  ggplot() + 
  geom_ribbon(aes(x = Time, ymin = MorR_L, ymax = MorR_U), alpha = 0.3) +
  geom_line(aes(x = Time, y = MorR_M)) + 
  scale_y_continuous("Mortality, per 100k", labels = scales::number_format(scale = 1e5)) + 
  scale_x_continuous("Year", breaks = seq(2000, 2020, 10)) + 
  expand_limits(y = 0) + 
  facet_grid(.~Scenario)



