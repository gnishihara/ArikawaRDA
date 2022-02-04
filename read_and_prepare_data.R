################################################################################
# Pew Marine Fellow
# Develop environmental health index
# Greg Nishihara
# February 2022
################################################################################
# Read packages #########################################################
library(tidyverse) # tidyverse!
library(readxl) # エクセルファイルの読み込み
library(lubridate) # 時刻データ処理
library(emmeans) # ペアごとの比較に使う
library(statmod) # qresiduals()
library(broom)   # tidy() と glance()
library(ggpubr)     # 論文用 ggplot テーマ
library(patchwork)  # 複数図のアレンジ
library(lemon)      # ggplot の詳細設定
library(magick)     # imagemagick の R API
library(showtext) 　# ggplot のフォントを指定する
library(furrr)
library(gnnlab)

number_of_cpu_cores = parallel::detectCores()
plan(multisession, workers = number_of_cpu_cores %/% 2) 

if(!any(str_detect(font_families(), "notosans"))) {
  # 英語フォント
  font_add_google("Noto Sans","notosans")
}
if(!any(str_detect(font_families(), "notosans-jp"))) {
  # 日本語フォント
  font_add_google("Noto Sans JP","notosans^jp")
}
theme_pubr(base_size = 10, base_family = "notosans") |> theme_set() # ggplot のデフォルト設定
showtext_auto() # 自動的にggplot にフォントを入れるようにする

# 自分のRサーバは設定は日本語になっているので、時間 locale を設定する
if (str_detect(Sys.getlocale(category = "LC_TIME"), "ja_JP.UTF-8")) {
  Sys.setlocale("LC_TIME", "en_US.UTF-8") # This is to set the server time locate to en_US.UTF-8
}

################################################################################
# 関数の定義 ###################################################################

boltzman = function(temperature, mean_temperature) {
  boltzman_constant = 8.617e-5 # ev / K
  activation_e = 0.65 #eV
  iK = 1/(temperature + 273.15)
  iKm = 1/(mean_temperature + 273.15)
  # boltzman_iK = exp(-1 / boltzman_constant * (iK - iKm))
  exp(-1 * (iK - iKm) / boltzman_constant)
}


masstransfer = function(windspeed, temperature, salinity, oxygen, height = 1) {
  # 大気と海面における酸素の輸送の計算
  # marelac パッケージが必要
  # height は m, salinity は PSU
  calc_k600 = function(windspeed, height) {
    # Crusius and Wanninkhof 2003 L&O 48
    cf = 1 + sqrt(1.3e-3)/0.4 * (log(10/height))
    U10 = cf * windspeed # m/s
    0.228*U10^2.2 + 0.168 # cm/h
  }
  k600 = calc_k600(windspeed, height)
  SCoxygen = marelac::gas_schmidt(temperature, species = "O2")
  a = ifelse(windspeed < 3, 2/3, 1/2)
  kx = k600 * (600/SCoxygen)^a # cm / h
  o2sat = marelac::gas_O2sat(salinity, temperature, method="Weiss")
  kx/100 * (o2sat - oxygen) # g / m2 / h
}

calculate_rate=function(data, k = 10, byhour = TRUE){
  # dO/dt の計算
  # ここで溶存酸素濃度のデータから酸素変動速度の計算をする
  # mgcv パッケージが必要
  out = mgcv::gam(mgla ~ s(H, k = k, bs = "cs"),  data = data, family = mgcv::scat(link = "identity"))
  
  if(byhour) {
    data2 = tibble(H = seq(0, 23, by = 1))
  } else {
    data2 = tibble(H = seq(0, 24-(1/6), by = 1/6))
  }
  # ここGAMの接線を求めている
  y1 = predict(out, newdata = data2)
  eps = 1e-6
  y2 = predict(out, newdata = tibble(H=data2$H+eps))
  data %>% mutate(rate=(y2-y1)/eps, mgl_predict = y1) 
}

calculate_dawn_dusk = function(datetime, gpscoord) {
  # 光のデータが十分じゃない時、日中の長さを求められないので、
  # 薄暮と薄明は crepuscule で求める
  # maptools のパッケージが必要
  # solarDep = 6 is civil twilight （市民薄明）
  # solarDep = 18 is astronomical twilight （天文薄明）
  
  if(!is.matrix(gpscoord)) {stop("gpscoord は行列として渡す [x, y]")}
  
  tz(datetime) = "Japan"
  dawn = maptools::crepuscule(gpscoord, datetime, solarDep = 6, direction = "dawn", POSIXct = T)[,2]
  dusk = maptools::crepuscule(gpscoord, datetime, solarDep = 6, direction = "dusk", POSIXct = T)[,2]
  tz(dawn) = "UCT"
  tz(dusk) = "UCT"
  tibble(daylight = interval(dawn, dusk))
}


se = function(x, na.rm=FALSE) {
  # 標準誤差
  N = sum(!is.na(x))
  sd(x, na.rm) / sqrt(N - 1)
} 

# 関数の処理中にエラーが出た場合、結果をNAにする
try_calc_mt = possibly(masstransfer, otherwise = NA)
try_calc_rate = possibly(calculate_rate, otherwise = NA)

# Shin-kamigoto
KAWATE = "~/Lab_Data/kawatea/"

#データの読み込み---------------------
#調査期間
weather = read_rds("~/Lab_Data/weather/20150101_20220201_Nagasaki_JMA_Data.rds")

period_file = dir(KAWATE, "period_info", full = TRUE)
kikan = read_csv(period_file)
kikan = kikan |> mutate(interval = interval(start_date, end_date)) |> select(-c(comment, remarks))
kikan = kikan |> mutate(location = str_replace(location, "ariakwa", "arikawa"))
kikan = kikan |> distinct() |> group_nest(location)


process_cem = function(x) {
  dset = read_alec(x)
  dset %>%
    mutate(datetime = floor_date(datetime, "10 minutes")) %>%
    group_by(datetime) %>%
    summarise(across(c(speed, dir, ew, ns, temperature), mean))
}

cem = tibble(fnames = dir(str_glue("{KAWATE}/CEM"), full = TRUE, pattern = "arik")) %>%
  mutate(data = future_map(fnames, process_cem)) %>% unnest(data)

cem = cem %>% group_by(datetime) %>% 
  summarise(across(c(speed, dir, ew, ns, temperature),
                   mean))

cku = tibble(fnames = dir(str_glue("{KAWATE}/CKU"), full = TRUE, pattern = "arik")) %>%
  mutate(data = future_map(fnames, read_alec)) %>% unnest(data)

cku = cku %>% group_by(datetime) %>% 
  summarise(across(c(chla, turbidity, temperature),
                   mean))


microstation = tibble(fnames = dir(str_glue("{KAWATE}/Microstation"), full = TRUE, pattern = "arik")) %>%
  mutate(data = future_map(fnames, read_onset)) %>% unnest(data) %>% 
  mutate(datetime = floor_date(datetime, "minutes")) %>% 
  select(-fnames) %>% distinct()

oxygen = tibble(fnames = dir(str_glue("{KAWATE}/Oxygen/"), full = TRUE), 
                type = str_extract(fnames, "xlxs|csv")) %>%
  mutate(data = future_map(fnames, read_onset))

oxygen = oxygen %>% mutate(fnames = basename(fnames)) %>%
  separate(fnames, into = c("type", "id", "location", "position", "date")) %>%
  filter(!str_detect(position, "calibra"))

oxygen = oxygen %>%
  unnest(data) %>%
  select(location, position, datetime, oxygen = mgl, temperature) %>%
  mutate(Date = as_date(datetime))

oxygen |> pull(location) |> unique()

oxygen = oxygen |> mutate(location = recode(location, amamo = "arikawaamamo", arikawa = "arikawaamamo")) |> 
  group_nest(location)

oxygen = full_join(oxygen, kikan, by = "location") |> rename(data = data.x, kikan = data.y)

oxy_all = 
  oxygen |> 
  mutate(data = future_map2(data, kikan, function(X, K) {
    intervals = K |> pull(interval)
    X |> ungroup() |> 
      filter(map_lgl(datetime, ~any(.x %within% intervals))) |> 
      drop_na() |> 
      select(datetime, position, oxygen,　temperature, Date)
  }))


light = 
  tibble(fnames =dir(str_glue("{KAWATE}/Light/"), full = TRUE)) |> 
  mutate(data = future_map(fnames, read_odyssey))

light = light %>% mutate(fnames = basename(fnames)) %>% 
  separate(fnames, into = c("type", "id", "location", "position", "date"))

light = light %>% unnest(data) %>% 
  select(location,id, position, datetime, ppfd) %>% 
  mutate(Date = as.Date(datetime))

light = light |> group_nest(location)

light = fuzzyjoin::fuzzy_full_join(light, kikan, by = "location", str_detect)

light = light |> rename(location = location.x, data = data.x, kikan = data.y)

light = light |> mutate(data = future_map2(data, kikan, function(X, K) {
  intervals = K |> pull(interval)
  X |> ungroup() |> 
    filter(map_lgl(datetime, ~any(.x %within% intervals))) |> 
    drop_na(ppfd) |> 
    select(id, position, Date, datetime, ppfd)
})) |> 
  select(-location.y, -kikan) |> 
  unnest(data)

# Get light calibration coefficients
fnames = dir(str_glue("{KAWATE}/Odyssey_Calibration"), pattern = "[Cc]alibration*.*[Cc][Ss][Vv]", full.names = TRUE)
calib = tibble(fnames) |> 
  mutate(has_coefficients = map_lgl(fnames, \(x) {
    read_lines(x, n_max = 1) |> str_detect("coefficient")
    })) |> 
  filter(has_coefficients) |> 
  mutate(data = map(fnames, read_csv)) |> 
  unnest(data) |> mutate(id = as.character(id)) |> 
  group_by(id) |> 
  summarise(cf = mean(calibration_coefficient))

light_all = left_join(light, calib, by = "id") %>% 
  mutate(ppfd.water = ppfd*cf) %>% select(-cf)

# Join all of the data

tmp = full_join(
  cem %>% select(datetime, speed, dir, ew, ns, temperature.cem = temperature),
  cku %>% select(datetime, chla, turbidity, temperature.cku = temperature),
  by = "datetime")
# tmp = full_join(tmp, depth, by = "datetime")

tmp = full_join(
  tmp,
  microstation %>% select(datetime, wind, gust, ppfd.microstation = ppfd),
  by = "datetime"
)

envdata = full_join(
  tmp, 
  weather %>% select(datetime, rain, temperature.jma = temperature_air, 
                     wind.jma = wind, gust.jma = gust),
  by = "datetime"
)


folder = "~/Lab_Data/nishihara/"
write_rds(oxy_all,   str_glue("{folder}all_oxygen_data.rds"))
write_rds(envdata,   str_glue("{folder}all_other_data.rds"))
write_rds(light_all, str_glue("{folder}all_light_data.rds"))




















