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
alldata3 =  read_rds(str_glue("{labdatafolder}/nishihara/alldata3.rds"))


alldata3 = alldata3 |> mutate(month = month(date), year = year(date)) |> mutate(x = year + (month - 1)/12) |> 
  mutate(location = factor(location))

alldata3 = alldata3 |> 
  mutate(x = as.numeric(date) - min(as.numeric(date))) 

library(mgcv)

gout = gam(gep ~ s(x, k = 50,  bs = "gp") + s(month, k = 10, bs = "gp") + s(month, k = 10, bs = "gp", by = interaction(location)) + year + location, data = alldata3)
summary(gout)
x = seq(ymd("2017-01-01"), ymd("2021-12-31"), by = "day")

pdata = alldata3 |> 
  expand(location) |> 
  mutate(x = list(x)) |> 
  unnest(x) |> 
  mutate(year = year(x),
         month = month(x),
         x = as.numeric(x) - min(as.numeric(x))) 

pdata = pdata |> mutate(fit = predict(gout, newdata = pdata))

ggplot(alldata3) + 
  geom_point(aes(x = x,  y = gep, color = location)) +
  geom_line(aes(x = x, y = fit, color = location), 
            data = pdata, 
            size = 3)

################################################################################
alldata3 = alldata3 |> rename(ppfd = ppfd.microstation)
Xa = alldata3 |> dplyr::select(-c(location, date, month, gep, nep, er))
Xb = alldata3 |> dplyr::select( c(location, date, month, gep, nep, er))

# Linear Discriminant Analysis #################################################
lout = MASS::lda(location ~ temperature + nep + speed + wind + gust + ppfd, 
                 data = alldata3 |> mutate(across(where(is.numeric), ~scale(.)[,1])))
scores  = predict(lout, alldata3 |> mutate(across(where(is.numeric), ~scale(.)[,1])))

# LDA output
ldscores = scores$x |> as_tibble() |> 
  mutate(location = alldata3$location,
         year = alldata3$year,
         month = alldata3$month)

ldscores |> 
  ggplot() +
  geom_hex(aes(x = LD1, y = LD2), bins = 51) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_viridis_c(end = 0.9) + 
  facet_grid(cols = vars(location))
  
  


################################################################################
# Build data for RDA analysis

alldata4 = bind_cols(alldata3, as_tibble(scores$x)) |> mutate(year = year(date))

# Summarize the data into monthly values

adata = alldata4 |> filter(str_detect(location, "garamo"))
udata = urchin   |> filter(str_detect(location, "Garamo")) |> mutate(date = as_date(datetime))
qdata = quadrat  |> filter(str_detect(location, "Garamo")) |> mutate(date = as_date(datetime))

udata = udata |> group_by(date, species) |> summarise(counts = n(), .groups = "drop")
qdata = qdata |>
  mutate(element = str_to_sentence(element)) |> 
  group_by(date) |> 
  mutate(n = sum(count)) |> 
  group_by(date, element) |> 
  summarise(count = sum(count), n = first(n)) |> 
  mutate(rate = count / n)
qdata = qdata |> dplyr::select(-c(n, rate)) |> pivot_wider(names_from = element, values_from = count) |> ungroup()
udata = udata |> pivot_wider(names_from = species, values_from = counts) |> mutate(across(where(is.numeric), ~replace_na(., 0)))
adata = adata |> 
  mutate(date = as_date(floor_date(date, "month"))) |> 
  group_by(date) |> 
  summarise(across(where(is.numeric), mean))
aqudata = inner_join(qdata, udata, by = "date") |> inner_join(adata, by = "date")

# Convert to percent coverage or percent abundance
# Create the community matrix
Y = aqudata |> dplyr::select(matches("Cor|Macro|Turf|Ds|Hc")) |> 
  mutate(across(matches("Cor|algae|Turf"), ~. / 1000)) |>
  rowwise() |> 
  mutate(N = sum(c_across(!matches("Cor|algae|Turf")))) |> 
  mutate(across(!matches("Cor|algae|Turf|^N$"), ~./N)) |> dplyr::select(-N) |> ungroup() |> 
  as.matrix()

colnames(Y) = c("Coralline algae", "Macroalgae", "Turf algae", 
                "D. savignyi", "D. setosum", "H. crassipina")

# Hellinger transformation of the community matrix
library(vegan)

Yhat = decostand(Y, "hellinger") 
aqudata = aqudata |> 
  rename(Temperature = temperature, 
         Speed = speed, 
         Wind = wind,
         PPFD = ppfd.microstation,
         GEP = gep)

# RDA
routh_null = rda(Yhat ~ 1, data = aqudata)
routh      = rda(Yhat ~ Temperature + Speed + Wind + PPFD + GEP, data = aqudata)

# Create data for the annotations
z = anova(routh) |> as_tibble()
fval = sprintf("F['(%d, %d)'] == %0.4f", z$Df[1], z$Df[2], z$F[1])
r2 = sprintf("R[adj]^{2}==%0.4f", RsquareAdj(routh)$adj.r.squared)

sites =   fortify(routh, axes = 1:2, scaling = 2) |> filter(str_detect(Score, "sites")) |> mutate(date = aqudata$date) |> mutate(date = as.factor(date))
species = fortify(routh, axes = 1:2, scaling = 2) |> filter(str_detect(Score, "species"))
biplot  = fortify(routh, axes = 1:2, scaling = 2) |> filter(str_detect(Score, "biplot"))

# Environmental variables
envivar = biplot |> 
  mutate(sRDA1 = 1.5*RDA1 / sqrt((RDA1^2+RDA2^2))) |> 
  mutate(sRDA2 = 1.5*RDA2 / sqrt((RDA1^2+RDA2^2))) |> 
  mutate(theta = atan(sRDA2/ sRDA1)) |>  print() |> 
  mutate(theta = 180 * theta / pi + ifelse(RDA1 < 0, theta , theta),
         hjust = ifelse(RDA1 > 0, 1, 0),
         vjust = ifelse(str_detect(Label, "GEP"), 1, 0))

# Species community matrix
specvar = species |> 
  mutate(sRDA1 = 1.5*RDA1 / sqrt((RDA1^2+RDA2^2))) |> 
  mutate(sRDA2 = 1.5*RDA2 / sqrt((RDA1^2+RDA2^2))) |> 
  mutate(theta = atan(sRDA2/ sRDA1)) |>  print() |> 
  mutate(theta = 180 * theta / pi + ifelse(RDA1 < 0, theta , theta),
         hjust = ifelse(RDA1 > 0, 1, 0),
         vjust = ifelse(str_detect(Label, "crass"), 1, 0)) |> 
  mutate(Label = ifelse(str_detect(Label, "\\."), 
                        str_glue("italic('{Label}')"), 
                        str_glue("plain('{Label}')")))

# The ggplot
ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(x = RDA1, y = RDA2), data = sites2, size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = sRDA1, yend = sRDA2, color = Label), alpha = 0.25, size = 2.5, data = envivar) +
  geom_segment(aes(x = 0, y = 0, xend = sRDA1, yend = sRDA2, color = Label), alpha = 0.25, size = 2.5, data = specvar) +
  geom_text(aes(x = sRDA1,  y = sRDA2, label = Label, color = Label, angle = theta, hjust = hjust, vjust = vjust), data = envivar, family = "notosans")  +
  geom_text(aes(x = sRDA1,  y = sRDA2, label = Label, color = Label, angle = theta, hjust = hjust, vjust = vjust), data = specvar, family = "notosans", parse = TRUE)  +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, color = Label), data = biplot, arrow = arrow(15, unit(1.5, "mm"))) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, color = Label), data = species, arrow = arrow(15, unit(1.5, "mm"))) +
  annotate("text",x = -1, y = 1, label = r2, vjust = 0, hjust = 0, parse = T) +
  annotate("text", x = -1, y = 1, label = fval, vjust = 1,hjust = 0, parse = T) +
  scale_color_viridis_d(end = 0.9) + guides(color = "none") + 
  scale_x_continuous(parse(text = "RDA[1]")) +
  scale_y_continuous(parse(text = "RDA[2]")) +
  coord_equal() +
  labs(caption = "Correlation biplot") + theme(legend.position = "left") +
  theme(panel.background =  element_rect(color = "black", size = 1),
        axis.line = element_blank())

pdfname = str_glue("RDA-plot-{today()}.pdf")
pngname = str_replace(pdfname, "pdf", "png")
ggsave(pdfname, width = 160, height = 160, units = "mm")
magick::image_read_pdf(pdfname, density = 600) |> magick::image_write(pngname, density = 300)



