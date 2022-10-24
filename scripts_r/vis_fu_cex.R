library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


pop.size <- 3e6


stats <- local({
  cost <- read_csv(here::here("data", "cost.csv"))
  cost <- as.list(setNames(cost$Vul, cost$Item))
  
  stats <- read_csv(here::here("out", "dy_hi", "Sim_VulACF_cex_stats.csv"))[-1] %>% 
    mutate(
      N_Screened = ACF_Vul_Screened,
      N_Confirmed = ACF_Vul_Confirmed,
      N_Fl_DS = ACF_Vul_DS_Fl,
      N_Fl_DR = ACF_Vul_DR_Fl,
      N_Sl_DR = ACF_Vul_DR_Sl,
      N_TPT = ACF_Vul_Yield - N_Fl_DS - N_Fl_DR - N_Sl_DR,
      C_Screened = ACF_Vul_Vul * cost$Vul + ACF_Vul_Sym * cost$Sym + ACF_Vul_CXR * cost$CXR,
      C_Confirmed = ACF_Vul_Xpert * cost$Xpert,
      C_Tx = (N_Fl_DS + N_Fl_DR + N_TPT) * cost$Tx_Fl + (N_Fl_DR + N_Sl_DR) * cost$Tx_Sl,
      C_Total = C_Screened + C_Confirmed + C_Tx,
      N_Screened = ACF_VulFu_Screened,
      N_Confirmed = ACF_VulFu_Confirmed,
      N_Fl_DS = ACF_VulFu_DS_Fl,
      N_Fl_DR = ACF_VulFu_DR_Fl,
      N_Sl_DR = ACF_VulFu_DR_Sl,
      N_TPT = ACF_VulFu_Yield - N_Fl_DS - N_Fl_DR - N_Sl_DR,
      C_ScreenedFu = ACF_VulFu_Sym * cost$Sym + ACF_VulFu_CXR * cost$CXR,
      C_ConfirmedFu = ACF_VulFu_Xpert * cost$Xpert,
      C_TxFu = (N_Fl_DS + N_Fl_DR + N_TPT) * cost$Tx_Fl + (N_Fl_DR + N_Sl_DR) * cost$Tx_Sl,
      C_TotalFu = C_ScreenedFu + C_ConfirmedFu + C_TxFu,
      PPV_Vul = ifelse(ACF_Vul_Yield > 0, ACF_Vul_TP / ACF_Vul_Yield, 0),
      PPV_Vulfu = ifelse(ACF_VulFu_Yield > 0, ACF_VulFu_TP / ACF_VulFu_Yield, 0)
    ) %>% 
    select(Key, Scenario, FollowUp, Duration, Coverage, Pop0, Inc1 = IncR, Mor1 = MorR,
           N_Screened, N_Confirmed, starts_with("C_"), starts_with('PPV'))
  
  
  stats <- stats %>% 
    filter(Scenario != "Baseline") %>% 
    left_join(stats %>% 
                filter(Scenario == "Baseline") %>%
                mutate(
                  C0_Total = C_Total + C_TotalFu,
                  C0_Test = C_Screened + C_Confirmed + C_ScreenedFu + C_ConfirmedFu,
                ) %>% 
                select(Key, Inc0 = Inc1, Mor0  = Mor1, C0_Total, C0_Test)) %>% 
    mutate(
      C1_Total = C_Total + C_TotalFu,
      C1_Test = C_Screened + C_Confirmed + C_ScreenedFu + C_ConfirmedFu,
      dE = Inc0 - Inc1,
      AvtInc = dE / Inc0,
      dC_Total = C1_Total - C0_Total,
      dC_Test = C1_Test - C0_Test,
      CER_Total = dC_Total / dE,
      CER_Test = dC_Test / dE,
      dC_Total = dC_Total / Pop0 * pop.size,
      dC_Test = dC_Test / Pop0 * pop.size 
    )
  
  stats %>% 
    group_by(Duration, Coverage) %>% 
    summarise(across(c(AvtInc, CER_Total, CER_Test, dC_Total, dC_Test, starts_with("PPV"), starts_with("C_")), list(
      M = function(x) mean(x, na.rm=T),
      L = function(x) quantile(x, 0.025, na.rm=T),
      U = function(x) quantile(x, 0.975, na.rm=T)
    ))) %>% 
    ungroup() %>% 
    mutate(
      #n_fu = ifelse(FollowUp == 0, 0, 1 / FollowUp)
    )
}) %>% 
  pivot_longer(-c(Duration, Coverage)) %>% 
  extract(name, c("Index", "Stats"), "(\\S+)_(M|L|U)") %>% 
  pivot_wider(names_from = Stats) %>% 
  arrange(Index, Coverage, Duration)



stats %>% 
  filter(Index == "AvtInc") %>% 
  ggplot() +
  geom_contour(aes(x = Duration, y= Coverage, z = M))


stats %>% 
  filter(Index %in% c("AvtInc", "dC_Total")) %>% 
  select(Duration = Duration, Coverage, name = Index, value = M) %>% 
  pivot_wider() %>% 
  ggplot() +
  geom_line(aes(x = dC_Total, y = AvtInc, colour = as.factor(Duration))) + 
  geom_point(aes(x = dC_Total, y = AvtInc, colour = as.factor(Duration)))








