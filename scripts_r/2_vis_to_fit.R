library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


targets <- jsonlite::read_json(here::here("data", "targets.json"))
targets <- bind_rows(targets)



folders <- c("main")


for(folder in folders) {
  sims <- read_csv(here::here("out", folder, "Sim_ToFit.csv"))
  
  
  gs <- list(
    sims %>% 
      ggplot() +
      geom_line(aes(x = Time, y = Prev, group = Key, colour = "Simulation"), alpha = 0.2) +
      expand_limits(y = 0) +
      geom_pointrange(data = targets %>% filter(Index == "Prev"), 
                      aes(x = 2020, y = M, ymin = L, ymax = U)) + 
      scale_y_continuous("Per 100 000", labels=scales::number_format(scale=1e5)) +
      labs(subtitle = "Untreated prevalent TB") + theme(legend.position = "none"),
    sims %>% 
      ggplot() +
      geom_line(aes(x = Time, y = PrAsym, group = Key, colour = "Simulation"), alpha = 0.2) +
      expand_limits(y = 0) +
      geom_pointrange(data = targets %>% filter(Index == "PrAsym"), 
                      aes(x = 2020, y = M, ymin = L, ymax = U)) + 
      scale_y_continuous("Percentage", labels=scales::percent) +
      labs(subtitle = "Untreated TB, asymptomatic") + theme(legend.position = "none"),
    sims %>% 
      ggplot() +
      geom_line(aes(x = Time, y = PrExCS, group = Key, colour = "Simulation"), alpha = 0.2) +
      expand_limits(y = 0) +
      geom_pointrange(data = targets %>% filter(Index == "PrExCS"), 
                      aes(x = 2020, y = M, ymin = L, ymax = U)) + 
      scale_y_continuous("Percentage", labels=scales::percent) +
      labs(subtitle = "Untreated TB, sought care") + theme(legend.position = "none"),
    sims %>% 
      ggplot() +
      geom_line(aes(x = Time, y = ARTI_Sp, group = Key, colour = "Simulation"), alpha = 0.2) +
      expand_limits(y = 0) +
      geom_pointrange(data = targets %>% filter(Index == "ARTI"), 
                      aes(x = 2020, y = M, ymin = L, ymax = U)) + 
      scale_y_continuous("Percentage", labels=scales::percent) +
      labs(subtitle = "Annual risk of TB infection") + theme(legend.position = "none"),
    sims %>% 
      ggplot() +
      geom_line(aes(x = Time, y = PrDR_CNR, group = Key, colour = "Simulation"), alpha = 0.2) +
      expand_limits(y = 0) +
      geom_pointrange(data = targets %>% filter(Index == "PrDR_CNR"), 
                      aes(x = 2020, y = M, ymin = L, ymax = U)) + 
      scale_y_continuous("Percentage", labels=scales::percent) +
      labs(subtitle = "DR-TB among case notification") + theme(legend.position = "none")
  )
  
  
  g_all <- ggarrange(gs[[1]], gs[[4]], gs[[2]], gs[[5]], gs[[3]], NULL, ncol=2, nrow=3)
  
  ggsave(g_all, filename = here::here("out", folder, "g_tofit.png"), width = 6, height = 7)
}


ggarrange(gs[[1]], gs[[4]], gs[[2]], gs[[5]], gs[[3]], NULL, ncol=2, nrow=3)

