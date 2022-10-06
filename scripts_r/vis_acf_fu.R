library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


pop.size <- 3e6


stats <- local({
  cost <- read_csv(here::here("data", "cost.csv"))
  cost <- as.list(setNames(cost$M, cost$Item))
  
  stats <- read_csv(here::here("out", "dy_hi", "Sim_VulACF_fudur_stats.csv"))[-1] %>% 
    mutate(
      N_Screened = ACF_Vul_Screened,
      N_Confirmed = ACF_Vul_Confirmed,
      N_Fl_DS = ACF_Vul_DS_Fl,
      N_Fl_DR = ACF_Vul_DR_Fl,
      N_Sl_DR = ACF_Vul_DR_Sl,
      N_TPT = ACF_Vul_Yield - N_Fl_DS - N_Fl_DR - N_Sl_DR,
      C_Screened = N_Screened * cost$Vul + cost$CXR,
      C_Confirmed = N_Confirmed * cost$Xpert,
      C_Tx = (N_Fl_DS + N_Fl_DR + N_TPT) * cost$Tx_Fl + (N_Fl_DR + N_Sl_DR) * cost$Tx_Sl,
      C_Total = C_Screened + C_Confirmed + C_Tx,
      N_Screened = ACF_Vulfu_Screened,
      N_Confirmed = ACF_Vulfu_Confirmed,
      N_Fl_DS = ACF_Vulfu_DS_Fl,
      N_Fl_DR = ACF_Vulfu_DR_Fl,
      N_Sl_DR = ACF_Vulfu_DR_Sl,
      N_TPT = ACF_Vulfu_Yield - N_Fl_DS - N_Fl_DR - N_Sl_DR,
      C_ScreenedFu = N_Screened * cost$Vul + cost$CXR,
      C_ConfirmedFu = N_Confirmed * cost$Xpert,
      C_TxFu = (N_Fl_DS + N_Fl_DR + N_TPT) * cost$Tx_Fl + (N_Fl_DR + N_Sl_DR) * cost$Tx_Sl,
      C_TotalFu = C_ScreenedFu + C_ConfirmedFu + C_TxFu,
      PPV_Vul = ifelse(ACF_Vul_Yield > 0, ACF_Vul_TP / ACF_Vul_Yield, 0),
      PPV_Vulfu = ifelse(ACF_Vulfu_Yield > 0, ACF_Vulfu_TP / ACF_Vulfu_Yield, 0),
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
    group_by(FollowUp, Duration) %>% 
    summarise(across(c(AvtInc, CER_Total, CER_Test, dC_Total, starts_with("PPV"), starts_with("C_")), list(
      M = mean,
      L = function(x) quantile(x, 0.025),
      U = function(x) quantile(x, 0.975)
    ))) %>% 
    ungroup() %>% 
    mutate(
      n_fu = ifelse(FollowUp == 0, 0, 1 / FollowUp)
    )
}) %>% 
  pivot_longer(-c(FollowUp, Duration, n_fu)) %>% 
  extract(name, c("Index", "Stats"), "(\\S+)_(M|L|U)") %>% 
  pivot_wider(names_from = Stats) %>% 
  arrange(Index, Duration, n_fu)




s0 <- stats %>% 
  filter(Duration == 0) %>% 
  mutate(Duration = as.factor(Duration))

s1 <- stats %>% 
  filter(Duration %in% 1:6) %>% 
  mutate(Duration = as.factor(Duration))





stats %>% 
  filter(startsWith(Index, "C_")) %>% 
  filter(!(Index %in% c("C_Total", "C_TotalFu"))) %>% 
  filter(Duration %in% 3:4) %>% 
  ggplot() +
  geom_histogram(aes(x = Duration, y = M, fill = Index), stat = "identity") +
  facet_wrap(n_fu~.)




g_fudur <- s1 %>% 
  filter(Index == "AvtInc") %>% 
  ggplot() +
  geom_line(aes(x = n_fu, y = M, colour = Duration, group = Duration)) +
  geom_point(aes(x = n_fu, y = M, colour = Duration, group = Duration)) +
  geom_hline(data = s0 %>% filter(Index == "AvtInc"), aes(yintercept = M[1]), linetype = 2) + 
  geom_hline(aes(yintercept = L[1]), linetype = 3) + 
  geom_hline(aes(yintercept = U[1]), linetype = 3) + 
  geom_text(aes(x = 6, y = M[1]), vjust = 1.5, hjust = 1, label = "No follow-up") + 
  scale_y_continuous("Incident case averted, %, 2023-2030", labels = scales::percent) +
  scale_x_continuous("Tests per year", breaks = seq(0, 6, 1)) +
  scale_colour_discrete("Follow-up period, year", guide = guide_legend(reverse = T)) + 
  expand_limits(y = 0, x = 1) + 
  theme(legend.position = c(1, 0), legend.justification = c(1.05, -0.05))


g_fudur



g_fudur_ce <- s1 %>% 
  filter(Index == "CER_Total") %>% 
  ggplot() +
  geom_line(aes(x = n_fu, y = M, colour = Duration, group = Duration)) +
  geom_hline(data = s0 %>% filter(Index == "CER_Total"), aes(yintercept = M[1]), linetype = 2) + 
  geom_hline(aes(yintercept = L[1]), linetype = 3) + 
  geom_hline(aes(yintercept = U[1]), linetype = 3) + 
  geom_text(aes(x = 6, y = M[1]), vjust = 1.5, hjust = 1, label = "No follow-up") + 
  scale_y_continuous("Cost per averted case, US dollar") + 
  scale_x_continuous("Tests per year", breaks = seq(0, 6, 1)) +
  scale_colour_discrete("Follow-up period, year") + 
  expand_limits(y = 0, x = 1) +
  theme(legend.position = "None")


g_fudur_ce


g_fudur_cost <- s1 %>% 
  filter(Index %in% c("AvtInc", "dC_Total")) %>% 
  select(n_fu, Duration, Index, M) %>% 
  pivot_wider(names_from = Index, values_from = M) %>% 
  ggplot() +
  geom_line(aes(y = AvtInc, x = dC_Total, colour = as.factor(Duration))) +
  geom_point(aes(y = AvtInc, x = dC_Total, colour = as.factor(Duration))) +
  geom_hline(data = s0 %>% filter(Index == "AvtInc"), aes(yintercept = M[1]), linetype = 2) + 
  geom_hline(data = s0 %>% filter(Index == "AvtInc"), aes(yintercept = L[1]), linetype = 3) + 
  geom_hline(data = s0 %>% filter(Index == "AvtInc"), aes(yintercept = U[1]), linetype = 3) + 
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", breaks=seq(5, 25, 5) * 1e6, labels = scales::number_format(scale = 1e-6)) + 
  # scale_y_continuous("Incident case averted, %, 2023-2030", labels = scales::percent) +
  scale_y_continuous("", labels = scales::percent) +
  scale_colour_discrete("Follow-up period, year", guide = guide_legend(reverse = T)) + 
  expand_limits(y = 0) +
  theme(legend.position = "None")


g_fudur_cost  


g_bind1 <- ggpubr::ggarrange(g_fudur, g_fudur_ce, nrow = 1)

g_bind2 <- ggpubr::ggarrange(g_fudur, g_fudur_cost, nrow = 1)


ggsave(g_fudur, filename = here::here("docs", "g_fu_durinc.png"), width = 7, height = 5)
ggsave(g_fudur_ce, filename = here::here("docs", "g_fu_durce.png"), width = 7, height = 5)

ggsave(g_bind1, filename = here::here("docs", "g_fu_dur1.png"), width = 9, height = 5)
ggsave(g_bind2, filename = here::here("docs", "g_fu_dur2.png"), width = 9, height = 5)

