library(tidyverse)
library(jsonlite)


theme_set(theme_bw())


ds <- dir("out")
ds <- ds[startsWith(ds, "dy_")]

print(ds)


for (d in ds){
  print(d)
  
  pars <- read_csv(here::here("out", d, "Post.csv")) %>% 
    rename(Key = `...1`)
  
  sims <- read_csv(here::here("out", d, "Runs_IntvE.csv")) %>% 
    mutate(Year = Time)
  
  sims_acc <- sims %>% 
    select(Time, Pop, IncR, MorR, matches("N_ACF_(\\w+)_Reached"), starts_with("PPV_"), Scenario, Key) %>% 
    # filter(Key < 20) %>% 
    group_by(Key, Scenario) %>% 
    mutate(
      Pop = Pop / Pop[1] * 1e5,
      N_ACF_TB_Reached = N_ACF_TB_Reached / Pop[1] * 1e5,
      N_ACF_NonTB_Reached = N_ACF_NonTB_Reached / Pop[1] * 1e5
    ) %>% 
    summarise(
      Inc = sum(IncR * Pop) * 0.5,
      Mor = sum(MorR * Pop) * 0.5,
      Reached_TB = sum(N_ACF_TB_Reached) * 0.5,
      Reached_NonTB = sum(N_ACF_NonTB_Reached) * 0.5,
      PPV_Acf = PPV_Acf[length(PPV_Acf)],
      PPV_Pcf = PPV_Pcf[length(PPV_Pcf)]
    )
  
  
  save(sims_acc, pars, file = here::here("out", d, "Acc_Intv.rdata"))
}
