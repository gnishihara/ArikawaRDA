# Calculate the NEP data for Arikawa
# 2022 Feb
# greg
# 
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


################################################################################
# Read data files  #############################################################
labdatafolder = "~/Lab_Data"
oxy_all   = read_rds(str_glue("{labdatafolder}/nishihara/all_oxygen_data.rds"))
envdata   = read_rds(str_glue("{labdatafolder}/nishihara/all_other_data.rds"))
light_all = read_rds(str_glue("{labdatafolder}/nishihara/all_light_data.rds"))
quadrat = read_csv(str_glue("{labdatafolder}/bellezad/write_out/quadrat_data.csv"))
urchin  = read_csv(str_glue("{labdatafolder}/bellezad/write_out/urchin_data.csv"))


# Arikawa Bay GPS coordinates from Google Maps
yx = c(32.998570964410106, 129.10076644664358) 
gpscoord = matrix(rev(yx), ncol = 2)

# Check the times and fix them
oxy_all = oxy_all |> select(location, data) |> unnest(data) 
oxy_all   = oxy_all |> mutate(datetime = floor_date(datetime, "10 minutes"))
light_all = light_all |> mutate(datetime = floor_date(datetime, "10 minutes"))

# Remove mushima data
oxy_all = oxy_all |> filter(!str_detect(location, "mushima"))

# Determine the temperature and salinity corrected oxygen concentrations
oxy_all = oxy_all |> rename(mgl = oxygen)  |> 
  mutate(H = hour(datetime) + minute(datetime) / 60) |> 
  mutate(mglper = mgl /  marelac::gas_O2sat(S = 0,  t = temperature, method="Weiss")) |> 
  mutate(mgla = mglper * marelac::gas_O2sat(S = 35, t = temperature, method="Weiss"))

# Remove data with less than 50 data points
# Remove data that does not include the 00:00:00 and 23:50:00
oxy_all = oxy_all |> 
  group_by(location, position, Date) |> 
  filter(n() > 50) |> 
  group_nest() |> 
  mutate(duration = future_map_dbl(data, \(x) {
    difftime(max(x$datetime), min(x$datetime), units = "hours") |> as.numeric()
  })) |> 
  filter(duration > 23+(4/6)) |> select(-duration)

# This procedure can take a few hours. #########################################
# Calculate the NEP rates
cfilename = str_glue("{labdatafolder}/nishihara/calculated_rates.rds")
if(!file.exists(cfilename)) {
  plan(multisession, workers = 16+8) 
  oxy_all = oxy_all |>  
    mutate(rate = future_map(data, try_calc_rate, byhour = FALSE, k = 48))
  oxy_all |> write_rds(cfilename)
} else {
  oxy_all = read_rds(cfilename)
}

oxy_all = oxy_all |>  select(-data) |>  unnest(rate)
alldata = inner_join(oxy_all, envdata, by = "datetime")
alldata = alldata |> select(-c(dir, rain))
alldata = alldata |> mutate(ppfd.microstation = ppfd.microstation - 1.2)

afiledata = str_glue("{labdatafolder}/nishihara/alldata3.rds")
if(!file.exists(afiledata)) {
  X = alldata |> 
    mutate(across(matches("mgl|temperature|speed|chla|turb|wind|gust|ppfd"), ~log(.+0.1))) |> 
    dplyr::select(where(is.numeric))
  
  # USE PCA to impute missing values.
  # Better than MICE?
  # This can take a few minutes
  xfilename = str_glue("{labdatafolder}/nishihara/imputed_pca_data.rds")
  if(!file.exists(xfilename)) {
    Xout = missMDA::imputePCA(X)
    Xout |>  write_rds(xfilename)
  } else {
    Xout = read_rds(xfilename)
  }
  
  alldata2 = alldata |> dplyr::select(location, position, datetime, date = Date) |> bind_cols(as_tibble(Xout$completeObs))
  alldata2 = alldata2 |> 
    mutate(across(matches("mgl|temperature|speed|chla|turb|wind|gust|ppfd"), ~ round(exp(.) - 0.1, 3) ))
  alldata2 = alldata2 |> mutate(wind = ifelse(wind < 0, 0, wind))
  
  alldata2 = 
    alldata2 |> 
    group_nest(location, position) |> 
    mutate(data = future_map(data, \(df) {
      df |> mutate(mt = try_calc_mt(wind, temperature, 35, mgla)) 
    })) |> unnest(data)
  
  alldata2 = alldata2 |> 
    group_by(location, datetime) |> 
    summarise(across(c(temperature, rate), mean),
              across(c(speed, ew, ns,
                       temperature.cem, temperature.cku, wind, gust, ppfd.microstation, 
                       temperature.jma, wind.jma, gust.jma, mt), first))
  
  alldata2 = alldata2 |> ungroup() |> 
    mutate(date = as_date(datetime)) |> 
    group_nest(location, date) |> 
    mutate(daylight = calculate_dawn_dusk(date, gpscoord = gpscoord)) |> 
    unnest(data)
  
  alldata3 = alldata2  |> 
    mutate(isday = datetime %within% daylight, .after = "datetime") |> 
    group_by(location, date) |> 
    summarise(across(c(temperature, speed, wind, gust, wind.jma, gust.jma), mean, na.rm=T),
              across(c(ppfd.microstation), sum),
              nep = sum(rate - mt),
              er  = sum((rate - mt) * (!isday)) / sum(!isday) * 144)
  
  alldata3 = alldata3 |> ungroup() |> mutate(gep = nep - er)
  alldata3 = alldata3 |> mutate(month = floor_date(date, "month")) 
  alldata3 |> write_rds(afiledata)
} else {
  alldata3 =  read_rds(afiledata)
}



