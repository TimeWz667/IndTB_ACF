library(tidyverse)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))




ds <- dir("out")
ds <- ds[startsWith(ds, "dy_sc")]
ds <- setNames(ds, ds)
ds <- ds[ds != "dy_sc2-1"]



post <- read_csv(here::here("out", ds[1], "Runs_ACF_t.csv"))
colnames(post)



incs <- post %>% 
  select(Year = Time, Pop, IncR, Key, Scenario, Key) %>% 
  group_by(Key, Scenario) %>% 
  arrange(Year) %>% 
  mutate(
    Inc = Pop * IncR,
    AccInc = cumsum(Inc)
  ) %>% 
  ungroup()




incs0 <- incs %>% filter(Scenario == "Baseline") %>% select(Year, Key, AccInc0 = AccInc)
incs1 <- incs %>% filter(Scenario != "Baseline") %>% select(Year, Key, AccInc1 = AccInc, Scenario)



g <- incs1 %>% 
  left_join(incs0) %>% 
  mutate(
    AvtInc = 1 - AccInc1 / AccInc0
  ) %>% 
  group_by(Scenario, Year) %>% 
  summarise(
    M = median(AvtInc),
    L = quantile(AvtInc, 0.25),
    U = quantile(AvtInc, 0.75)
  ) %>% 
  extract(Scenario, c("ACF", "Target", "Rate"), "(\\w+), (Universal|Focus), R_ACF=(\\S+)") %>% 
  mutate(
    Target = ifelse(Target == "Focus", "High risk population", "Universal")
  ) %>% 
  ggplot() + 
  geom_line(aes(x = Year, y = M, colour=Rate)) +
  scale_y_continuous("Averted incident cases, %", labels = scales::percent) +
  scale_x_continuous("Year", breaks = c(2020, 2025, 2030)) + 
  facet_grid(.~Target+ACF) +
  labs(caption = "OR. TB = 2.38, Pr(High risk) = 17.9%") 
#  labs(caption = "OR. TB = 34.7, Pr(High risk) = 0.7%")


ggsave(g, filename = here::here("out", "Sc-1-1.png"), width = 7, height = 5)
