# MARSS
# Multivariate Autoregressive State-Space Modeling with R
# https://nwfsc-timeseries.github.io/MARSS/#documentation
# 
# Playing around with the code. Nothing serious.
# 2022 Feb
# greg
# 

library(tidyverse)
library(MARSS)
library(data.table)
library(tsibble)

dset = oxy_all |> slice(1) |> unnest(data) |> select(-kikan)
dset = dset |> mutate(Date = floor_date(Date, "month"))

dset = as.data.table(dset)

oxygen = dset[, by = c("position", "Date"),
       lapply(.SD, mean, na.rm = TRUE),
       .SDcols = "oxygen"] |> as_tibble() |> 
  arrange(Date)
oxygen = oxygen |> pivot_wider(names_from = position, values_from = oxygen)

dat = tsibble(oxygen) |> mutate(Date = yearmonth(Date))

temperature = dset[, by = c("Date"),
                lapply(.SD, mean, na.rm = TRUE),
                .SDcols = "temperature"] |> as_tibble()

covariate = tsibble(temperature) |> mutate(Date = yearmonth(Date)) 
covariate = as.ts(covariate) |> as.vector() |> print()
covariate = (covariate - mean(covariate, na.rm = T)) / sd(covariate, na.rm = T)

covariate = nafill(covariate, "nocb")
t(covariate)
dat = as.ts(dat)

modlist = list(B = matrix(1),
               U = matrix("u"),
               Q = matrix("q"),
               Z = matrix(1,6,1),
               A = "scaling", 
               R = "diagonal and equal",
               d = t(covariate),
               D = "unconstrained",
               x0 = matrix("mu"), 
               tinitx = 0)

fit = MARSS(dat, model = modlist)
resids = MARSSresiduals(fit, type = "tt1")
summary(fit)

zz = oxygen |> pull(Date) |> range() 
plot(resids$model.residuals[3, ])


dout = tibble(fit = fit$states[1, ]) |> mutate(n = 1:n(),
                                               Date = seq(zz[1], zz[2], by = "month"))

o2 = oxygen |> mutate(n = 1:n()) |> pivot_longer(-c(n, Date))

tmp = o2 |> select(Date, n) |> distinct()

ggplot() + 
  geom_point(aes(x = Date, y = value, color = name), data = o2) +
  geom_line(aes(x = Date, y = fit), data = dout)  +
  geom_point(aes(x = Date, y = fit), data = dout)  




