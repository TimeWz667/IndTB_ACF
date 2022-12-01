library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


pop.size <- 3e6


stats <- local({
  cost <- read_csv(here::here("data", "cost.csv"))
  cost <- as.list(setNames(cost$Vul, cost$Item))
  
  stats <- read_csv(here::here("out", "main", "Sim_VulACF_fudur_stats.csv"))[-1] %>% 
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
      PPV_VulFu = ifelse(ACF_VulFu_Yield > 0, ACF_VulFu_TP / ACF_VulFu_Yield, 0),
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
      dC_Test = dC_Test / Pop0 * pop.size,
      across(starts_with("C_"), function(x) x / Pop0 * pop.size)
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



txt <- s1 %>% 
  filter(Index %in% c("AvtInc", "dC_Total")) %>% 
  filter(Duration == 4 & n_fu > 0) %>% 
  select(n_fu, Duration, Index, M) %>% 
  pivot_wider(names_from = Index, values_from = M)


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
  geom_text(data = txt, aes(y = AvtInc * 1.02, x = dC_Total * 0.99, label = n_fu)) +
  scale_x_continuous("Total ACF cost, in millions of 2019 USD", breaks=seq(8, 20, 2) * 1e6, labels = scales::number_format(scale = 1e-6, accuracy = 1)) + 
  # scale_y_continuous("Incident case averted, %, 2023-2030", labels = scales::percent) +
  scale_y_continuous("Incident case averted, %, 2023-2030", labels = scales::percent) +
  scale_colour_discrete("Follow-up period", guide = guide_legend(reverse = T), 
                        labels=function(x){paste0(x, ifelse(x=="1", " year", " years"))}) + 
  expand_limits(y = 0.05, x = 7e6) +
  theme(legend.position = "None") +
  labs(caption = "*Coverage=20% of total population per year\n*Numbers annotate the number of follow-up screening per year")


g_fudur_cost


g_fudur_cpart <- stats %>% 
  filter((n_fu == 1 & Duration ==4 ) | (n_fu == 4 & Duration == 1)) %>% 
  filter(startsWith(Index, "C_")) %>% 
  filter(!(Index %in% c("C_Total", "C_TotalFu"))) %>% 
  mutate(
    Type = ifelse(endsWith(Index, "Fu"), "Follow-up", "Initial ACF"),
    Scenario = ifelse(n_fu == 2, "4 years X 1", "1 years X 4"),
    Source = case_when(
      startsWith(Index, "C_Confirm") ~ "Confirmation",
      startsWith(Index, "C_Screened") ~ "Screening",
      startsWith(Index, "C_Tx") ~ "Treatment",
      T ~ "Total"
    ),
    Source = factor(Source, c("Screening", "Confirmation", "Treatment"))
  ) %>%
  ggplot() +
  geom_histogram(aes(x = Scenario, y = M, fill = Type), stat = "identity") +
  scale_y_continuous("ACF cost by type, in millions of 2019 USD", labels = scales::number_format(scale = 1e-6)) +
  scale_x_discrete("Follow-up duration X tests per year") +
  facet_wrap(.~Source, scales="free_y",) +
  labs(caption = "*Slum population size of 3 million assumed")


g_fudur_ppv <- stats %>% 
  filter((n_fu == 1 & Duration ==4) | (n_fu == 4 & Duration == 1)) %>% 
  filter(startsWith(Index, "PPV_"))  %>% 
  mutate(
    Type = ifelse(endsWith(Index, "Fu"), "Follow-up", "Initial ACF"),
    Type = factor(Type, c("Initial ACF", "Follow-up")),
    Scenario = ifelse(n_fu == 2, "4 years X 1", "1 years X 4")
  ) %>% 
  ggplot() +
  geom_histogram(aes(x = Scenario, y = M), stat = "identity", alpha = 0.7) +
  scale_y_continuous("PPV of case-yields", labels = scales::percent) +
  scale_x_discrete("Follow-up duration X tests per year") +
  facet_wrap(.~Type) +
  expand_limits(y = c(0, 0.5)) +
  labs(caption = "*Slum population size of 3 million assumed")


g_fudur_cost + 
  theme(legend.position = c(1, 0), legend.justification = c(1.05, -0.05))


g_bind1 <- ggpubr::ggarrange(g_fudur, g_fudur_ce, nrow = 1)

g_bind2 <- ggpubr::ggarrange(g_fudur, g_fudur_cost, nrow = 1)


ggsave(g_fudur, filename = here::here("docs", "figs", "g_fudur_inc.png"), width = 7, height = 5)
ggsave(g_fudur_cost + 
         theme(legend.position = c(1, 0), legend.justification = c(1.05, -0.05)), 
       filename = here::here("docs", "g_fudur_cost.png"), width = 7, height = 5)

ggsave(g_fudur_cpart, filename = here::here("docs", "figs", "g_fudur_cpart.png"), width = 9, height = 5)
ggsave(g_fudur_ppv, filename = here::here("docs", "figs", "g_fudur_ppv.png"), width = 7, height = 5)


ggsave(g_bind1, filename = here::here("docs", "figs", "g_fudur1.png"), width = 9, height = 5)
ggsave(g_bind2, filename = here::here("docs", "figs", "g_fudur2.png"), width = 9, height = 5)

