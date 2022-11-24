


dir.create(here::here("out", "main"), showWarnings = F)


for (file in dir(here::here("out", "dy_ladd"))) {
  print(file)
  file.copy(here::here("out", "dy_ladd", file), here::here("out", "main"))
}
