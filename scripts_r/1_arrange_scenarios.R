
dir.create(here::here("out", "main"), showWarnings = F)


for (file in dir(here::here("out", "dy_ladd"))) {
  print(file)
  file.copy(here::here("out", "dy_ladd", file), here::here("out", "main"), overwrite=TRUE)
}


dir.create(here::here("out", "sens"), showWarnings = F)


file.copy(here::here("out", "dy_ladd", "Sens_Pars_stats.csv"), here::here("out", "sens"), overwrite=TRUE)

file.copy(here::here("out", "dy_ladd", "Sim_VulACF_budget_stats.csv"), here::here("out", "sens"), overwrite=TRUE)
file.rename(here::here("out", "sens", "Sim_VulACF_budget_stats.csv"), 
            here::here("out", "sens", "Sim_VulACF_budget_stats_ladd.csv"))

file.copy(here::here("out", "dy_add", "Sim_VulACF_budget_stats.csv"), here::here("out", "sens"), overwrite=TRUE)
file.rename(here::here("out", "sens", "Sim_VulACF_budget_stats.csv"), 
            here::here("out", "sens", "Sim_VulACF_budget_stats_add.csv"))


file.copy(here::here("out", "dy_mul", "Sim_VulACF_budget_stats.csv"), here::here("out", "sens"), overwrite=TRUE)
file.rename(here::here("out", "sens", "Sim_VulACF_budget_stats.csv"), 
            here::here("out", "sens", "Sim_VulACF_budget_stats_mul.csv"))

