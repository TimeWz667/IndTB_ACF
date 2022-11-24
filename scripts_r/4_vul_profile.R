library(tidyverse)



n_vul <- 5277
n_mapped <- 34214



profile <- read_csv(here::here("data", "profile_vulnerable.csv"))
rr <- read_csv(here::here("data", "rr_vulnerable.csv"))

targets <- read_csv(here::here("data", "targets.csv"))
targets <- as.list(setNames(targets$M, targets$Index))


pop <- local({
  x <- profile %>% 
    mutate(
      Dust = str_detect(RF, "Dust"),
      Diabetes = str_detect(RF, "Diabetes"),
      Tobacco = str_detect(RF, "Tobacco"),
      HCW = str_detect(RF, "HCW"),
      Alc = str_detect(RF, "Alcohol"),
      CLD = str_detect(RF, "Chronic lung disease"),
      PastTB = str_detect(RF, "Past TB"),
      Contact = str_detect(RF, "HH Contact"),
      CLD_CKD = str_detect(RF, "CLD/CKD"),
      Cancer = str_detect(RF, "Cancer"),
      Total = Total * n_vul / sum(Total),
      Key = 1:n()
    ) %>% 
    pivot_longer(Dust:Cancer, names_to = "Name", values_to = "Exist") %>% 
    left_join(rr %>% select(Name, RR)) %>% 
    mutate(RR = ifelse(Exist, RR, 0)) %>% 
    filter(Exist) %>% 
    group_by(Key, Total) %>% 
    summarise(
      add = sum(RR),
      mul = prod(RR),
      ladd = 0.1 * sum(RR) + 0.9 * max(RR)
    )
  
  x <- bind_rows(x, tibble(Key = max(x$Key) + 1, Total = n_mapped - n_vul, add = 1, mul = 1, ladd = 1))
  
  x %>% ungroup()
})


populate <- function(df) {
  df %>% 
    mutate(
      Prev = targets$Prev * Total * RR / sum(Total * RR),
      TB_A = Prev * targets$PrAsym,
      TB_S = Prev * (1 - targets$PrAsym),
      NonTB = Total - Prev,
      NonTB_S = pmax(Total * 0.036 - TB_S, 0),
      NonTB_A = pmax(NonTB - NonTB_S, 0)
    ) %>% 
    select(Key, Total, RR, TB_A, TB_S, NonTB_A, NonTB_S)
}


calc_roc <- function(df, lab) {
  cs <- sort(df$RR)
  cs <- c(0, cs, max(cs) + 0.1)
  
  bind_rows(lapply(cs, function(cutoff) {
    df %>% 
      summarise(
        CutOff = cutoff,
        Sens = sum((TB_S + TB_A) * (RR > cutoff)) / sum((TB_S + TB_A)),
        Spec = sum((NonTB_S + NonTB_A) * (RR <= cutoff)) / sum((NonTB_S + NonTB_A))
      )
  })) %>% mutate(Assumption = lab)
}


calc_one <- function(df, lab) {
  df %>% 
    summarise(
      Sens = sum((TB_S + TB_A) * (RR > 1)) / sum((TB_S + TB_A)),
      Spec = sum((NonTB_S + NonTB_A) * (RR <= 1)) / sum((NonTB_S + NonTB_A))
    ) %>% mutate(Assumption = lab)
}


pop_ladd <- pop %>% 
  mutate(Total = Total / sum(Total)) %>% 
  select(Key, Total, RR = ladd) %>% 
  populate() %>% 
  arrange(RR)


pop_add <- pop %>% 
  mutate(Total = Total / sum(Total)) %>% 
  select(Key, Total, RR = add) %>% 
  populate() %>% 
  arrange(RR)


pop_mul <- pop %>% 
  mutate(Total = Total / sum(Total)) %>% 
  select(Key, Total, RR = mul) %>% 
  populate() %>% 
  arrange(RR)


rocs <- bind_rows(
  pop_ladd %>% calc_roc("ladd"),
  pop_add %>% calc_roc("add"),
  pop_mul %>% calc_roc("mul")
)


rocs0 <- bind_rows(
  pop_ladd %>% calc_one("ladd"),
  pop_add %>% calc_one("add"),
  pop_mul %>% calc_one("mul")
)


g_roc <- ggplot(rocs) +
  geom_line(aes(x = 1 - Spec, y = Sens, colour = Assumption)) +
  geom_point(data = rocs0, aes(x = 1 - Spec, y = Sens, colour = Assumption)) + 
  geom_abline(slope = 1, linetype = 2) +
  scale_x_continuous("1 - Specificity", labels = scales::percent) +
  scale_y_continuous("Sensitivity", labels = scales::percent) + 
  scale_color_discrete("Assumption for \nCombining RR.", 
                       label = c(add = "Additive", mul = "Multiplicative", ladd = "Less than Additive")) +
  theme(legend.position = c(1, 0), legend.justification = c(1.1, -0.1))

g_roc


ggsave(g_roc, filename = here::here("docs", "figs", "g_rocs.png"), width = 5, height = 5)  

