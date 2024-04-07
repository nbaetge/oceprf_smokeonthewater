Shipmet for Inline Analysis
================
Nick Baetge
compiled most recently on 04 April, 2024

``` r
library(tidyverse)
library(lubridate)
library(ncdf4)
library(patchwork)
```

# import and tidy data

``` r
metdata_path <-
  "/Users/nicholasbaetge/github/oceprf_ash/raw/r0_shipmet.csv"
r0_tsg_path <-
  "/Users/nicholasbaetge/github/oceprf_ash/InLineData/raw/TSG/tsg_081923.txt"
p0_shipmet <-
  "/Users/nicholasbaetge/github/oceprf_ash/prod/p0_shipmet.csv"
```

## shipmet for inline analysis

``` r
shipmet <- read_csv(metdata_path) %>%
  select(LA, LO, ate, Time, TT, TC, SA, SD, SP) %>%
  rename(
    lat = LA,
    lon = LO,
    ymd = ate,
    hms = Time,
    t1 = TT,
    c1 = TC,
    s = SA,
    sigma = SD,
    knots = SP
  ) %>%
  mutate(date = ymd(ymd)) %>%
  mutate(dmy = format(date, "%d%m%y")) %>%
  select(lat, lon, dmy, hms:knots)
```

# Save data

``` r
write_tsv(shipmet, r0_tsg_path)
write_csv(shipmet, p0_shipmet)
```
