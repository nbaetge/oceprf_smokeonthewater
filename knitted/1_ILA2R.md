1: Tidy Matlab product from Inline Analysis
================
Nick Baetge
compiled most recently on 20 November, 2024

``` r
library(tidyverse)
library(hms)
library(lubridate)
library(purrr)
library(patchwork)
library(janitor)
library(kohonen)
library(factoextra)
```

Here we tidy the bio-optical data from the inline system and create a
composite z-score index for high and low biomass stations. We’ll also
add in FRR and influx data

# data paths

``` r
acs_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r1_acs/"
bb_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r1_bb3/"
frr_data_path <- 
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r1_frr/frr_data.txt"
frr_blank_path <- 
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r1_frr/frr_blanks.txt"
influx_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r6_phyto.xlsx"
# data contains cells in units of per L 
shipmet_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p0_shipmet.csv"
exploc_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r1_exploc.csv"
prod_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p1_insitu_biooptics.csv"
index_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p1_insitu_index.csv"
bcodmo_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/bco_dmo/BIO-OPTICS.csv"
```

``` r
# wl <- 470, 532, 660 (bb wavelengths)
acs <-
  list.files(path = acs_path,
             pattern = "*.csv",
             full.names = T) %>%
  map_df(~ read_csv(.))
bb <-
  list.files(path = bb_path,
             pattern = "*.csv",
             full.names = T) %>%
  map_df(~ read_csv(.)) 
```

# identify locations and times near casts

``` r
cast_times <- read_csv(exploc_path) %>%
  select(stn, date, lat, lon) %>%
  rename(cast_lat = lat,
         cast_lon = lon) %>%
  mutate(cast_date = mdy(date)) %>%
  select(-date) %>%
  pivot_wider(names_from = stn,
              values_from = c(cast_date, cast_lat, cast_lon))


stn_times <- read_csv(shipmet_path) %>%
  mutate(
    date = dmy(dmy),
    datetime = paste(dmy, "T", hms, sep = ""),
    datetime = dmy_hms(datetime),
    local_datetime = with_tz(datetime, "America/Los_Angeles"),
    local_hms = as_hms(local_datetime)
  ) %>%
  bind_cols(., cast_times) %>%
  mutate(
    stn = if_else(
      # near(lat, cast_lat_4, 0.1) &
        # near(lon, cast_lon_4, 0.1) & 
        near(date, cast_date_4, 0.0625),
      4,
      NA
    ),
    stn = if_else(
      # near(lat, cast_lat_5, 0.1) &
        # near(lon, cast_lon_5, 0.1) & 
        near(date, cast_date_5, 0.0625),
      5,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_6, 0.1) &
        # near(lon, cast_lon_6, 0.1) & 
        near(date, cast_date_6, 0.0625),
      6,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_7, 0.1) &
        # near(lon, cast_lon_7, 0.1) & 
        near(date, cast_date_7, 0.0625),
      7,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_8, 0.1) &
        # near(lon, cast_lon_8, 0.1) & 
        near(date, cast_date_8, 0.0625),
      8,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_9, 0.1) &
        # near(lon, cast_lon_9, 0.1) & 
        near(date, cast_date_9, 0.0625),
      9,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_10, 0.1) &
        # near(lon, cast_lon_10, 0.1) &
        near(date, cast_date_10, 0.0625),
      10,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_11, 0.1) &
        # near(lon, cast_lon_11, 0.1) &
        near(date, cast_date_11, 0.0625),
      11,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_12, 0.1) &
        # near(lon, cast_lon_12, 0.1) &
        near(date, cast_date_12, 0.0625),
      12,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_13, 0.1) &
        # near(lon, cast_lon_13, 0.1) &
        near(date, cast_date_13, 0.0625),
      13,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_14, 0.1) &
        # near(lon, cast_lon_14, 0.1) &
        near(date, cast_date_14, 0.0625),
      14,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_15, 0.1) &
        # near(lon, cast_lon_15, 0.1) &
        near(date, cast_date_15, 0.0625) ,
      15,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_16, 0.1) &
        # near(lon, cast_lon_16, 0.1) &
        near(date, cast_date_16, 0.0625) ,
      16,
      stn
    ),
    stn = if_else(
      # near(lat, cast_lat_17, 0.1) &
        # near(lon, cast_lon_17, 0.1) &
        near(date, cast_date_17, 0.0625),
      17,
      stn
    )
  ) %>%
  drop_na(stn) %>%
  select(stn, datetime) %>%
  group_by(stn) %>%
  mutate(bin_datetime = round_date(datetime, "30 min")) %>%
  select(-datetime) %>%
  distinct() %>%
  filter(bin_datetime %in% c(first(bin_datetime), last(bin_datetime))) %>%
  mutate(start = bin_datetime[[1]],
         end = bin_datetime[[length(bin_datetime)]]) %>%
  select(stn, start, end) %>%
  distinct()
# mutate(date = date(start)) %>%
# mutate(stn_start = as_hms(start),
#        stn_end = as_hms(end)) %>%
# ungroup() %>%
# select(stn, date, stn_start, stn_end) 
```

# function to convert matlab time to R time

``` r
matlab2POS = function(x, timez = "UTC") {
    days = x - 719529   # 719529 = days from 1-1-0000 to 1-1-1970
    secs = days * 86400 # 86400 seconds in a day
    # Converting the secs value to a POSIXct value in the UTC time zone, then converts that to a time/date string that 
    # should lose the time zone, and then it performs a second as.POSIXct() conversion on the time/date string to get a POSIXct value in the user's specified timezone.
    return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
            tz = 'UTC'), format = '%Y-%m-%d %H:%M', 
            tz = 'UTC', usetz = FALSE), tz = timez))
}
```

# tidy frr data

``` r
frr <- read_table(frr_data_path) %>%
  select(DATE, TIME, Fo, Fm) %>%
  mutate_at(vars(TIME), as.character) %>%
  mutate(dt = mdy_hms(paste(DATE, TIME))) %>%
  select(dt, Fo, Fm) %>%
  full_join(
    .,
    read_table(frr_blank_path) %>%
      select(DATE, TIME, Fo) %>%
      mutate_at(vars(TIME), as.character) %>%
      mutate(dt = mdy_hms(paste(DATE, TIME))) %>%
      select(dt, Fo) %>%
      rename(blank_Fo = Fo)
  ) %>%
  arrange(dt) %>%
  mutate(
    blank_Fo = zoo::na.approx(blank_Fo, dt, na.rm = F),
    corr_Fo = Fo - blank_Fo,
    corr_Fm = Fm - blank_Fo,
    Fv = corr_Fm - corr_Fo,
    Fv_Fm = Fv / corr_Fm
  ) %>%
  drop_na(Fo) %>%
  mutate(dt = round_date(dt, "60 s")) %>%
  select(dt, Fv_Fm) 
```

# tidy influx data

``` r
cells <- readxl::read_xlsx(influx_path, sheet = 2) %>%
  select(influx_id:nano_tf) %>% 
  drop_na() %>% 
  # filter(trt == "CTL") %>% 
  select(influx_id, contains("t0")) %>% 
  rename(stn = influx_id) %>% 
  rename(syn = 2, 
         pico = 3, 
         nano = 4) %>% 
  mutate(nano_syn = nano/syn,
         pico_syn = pico/syn) %>% 
  select(stn, nano_syn, pico_syn) %>% 
  group_by(stn) %>% 
  summarize(across(
    .cols = c(nano_syn, pico_syn),
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  )) %>% 
  ungroup()
```

    ## Warning: There was 1 warning in `summarize()`.
    ## ℹ In argument: `across(...)`.
    ## ℹ In group 1: `stn = 4`.
    ## Caused by warning:
    ## ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
    ## Supply arguments directly to `.fns` through an anonymous function instead.
    ## 
    ##   # Previously
    ##   across(a:b, mean, na.rm = TRUE)
    ## 
    ##   # Now
    ##   across(a:b, \(x) mean(x, na.rm = TRUE))

``` r
cells_n <- readxl::read_xlsx(influx_path, sheet = 2) %>%
  select(influx_id:nano_tf) %>% 
  drop_na() %>% 
  # filter(trt == "CTL") %>% 
  select(influx_id, contains("t0")) %>% 
  rename(stn = influx_id) %>% 
  rename(syn = 2, 
         pico = 3, 
         nano = 4) %>% 
  mutate(nano_syn = nano/syn,
         pico_syn = pico/syn) %>% 
  select(stn, nano_syn, pico_syn) %>%
  group_by(stn) %>% 
  summarise(cells_n = n()) %>% 
  ungroup()
```

# tidy bio-optical data

restrict data to within 3 hours of cast time

``` r
tidy <- acs %>%
  select(dt,
         ap_16,
         ap_29,
         ap_56,
         cp_16,
         cp_29,
         cp_56,
         poc,
         chl_ap676lh,
         gamma) %>%
  mutate(
    bp_16 = cp_16 - ap_16,
    bp_29 = cp_29 - ap_29,
    bp_56 = cp_56 - ap_56,
    .after = dt
  ) %>%
  select(-c(ap_16:cp_56)) %>%
  rename(
    bp_470 = bp_16,
    bp_532 = bp_29,
    bp_660 = bp_56,
    poc_cp_660 = poc,
    gamma_cp = gamma
  ) %>%
  left_join(
    bb %>%
      select(dt, bbp_1:bbp_3, gamma_bbp, poc_1:cphyto_3) %>%
      rename(
        bbp_470 = bbp_1,
        bbp_532 = bbp_2,
        bbp_660 = bbp_3,
        poc_bbp_470 = poc_1,
        poc_bbp_532 = poc_2,
        poc_bbp_660 = poc_3,
        cphyto_470 = cphyto_1,
        cphyto_532 = cphyto_2,
        cphyto_660 = cphyto_3
      ),
    .
  ) %>%
  mutate(
    bbb_470 = bbp_470 / bp_470,
    bbb_532 = bbp_532 / bp_532,
    bbb_660 = bbp_660 / bp_660,
    poc_chl = poc_cp_660 / chl_ap676lh, 
  ) %>%
  select(dt, contains(c("gamma", "bbb", "poc", "cphyto", "ap676", "poc_chl"))) %>%
  rename(m_dt = dt) %>%
  mutate(dt = matlab2POS(m_dt), .before = m_dt) %>%
  select(-m_dt) %>%
  mutate(date = date(dt),
         time = as_hms(dt),
         .after = dt) %>%
  left_join(., frr) %>%
  mutate(stn = ifelse(between(dt, stn_times$start[stn_times$stn == 4], stn_times$end[stn_times$stn == 4]), 4, NA),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 5], stn_times$end[stn_times$stn == 5]), 5, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 6], stn_times$end[stn_times$stn == 6]), 6, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 7], stn_times$end[stn_times$stn == 7]), 7, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 8], stn_times$end[stn_times$stn == 8]), 8, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 9], stn_times$end[stn_times$stn == 9]), 9, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 10], stn_times$end[stn_times$stn == 10]), 10, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 11], stn_times$end[stn_times$stn == 11]), 11, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 12], stn_times$end[stn_times$stn == 12]), 12, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 13], stn_times$end[stn_times$stn == 13]), 13, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 14], stn_times$end[stn_times$stn == 14]), 14, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 15], stn_times$end[stn_times$stn == 15]), 15, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 16], stn_times$end[stn_times$stn == 16]), 16, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 17], stn_times$end[stn_times$stn == 17]), 17, stn)
    ) %>% 
  drop_na(stn) 
```

# summarize data

``` r
acs_n <- tidy %>%
  group_by(stn) %>%
  drop_na(chl_ap676lh) %>%
  summarise(acs_n = n())

summary <- tidy %>%
  select(stn, gamma_bbp:Fv_Fm) %>%
  group_by(stn) %>%
  summarize(across(
    .cols = gamma_bbp:Fv_Fm,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  )) %>%
  left_join(read_csv(exploc_path), .) %>% 
  left_join(., cells) %>% 
  select(
    stn,
    mean_chl_ap676lh,
    sd_chl_ap676lh,
    mean_gamma_cp,
    sd_gamma_cp,
    mean_bbb_532,
    sd_bbb_532,
    mean_poc_cp_660,
    sd_poc_cp_660,
    mean_poc_chl,
    sd_poc_chl,
    mean_nano_syn,
    sd_nano_syn,
    mean_pico_syn,
    sd_pico_syn,
  )  
```

## self-organizing maps

``` r
som <- summary %>%
  select(
    stn,
    mean_chl_ap676lh,
    mean_gamma_cp,
    mean_bbb_532,
    mean_poc_cp_660,
    mean_nano_syn,
    mean_pico_syn,
  ) %>%
  rename_all( ~ stringr::str_replace(., "mean_", "")) %>%
  # mutate(across(c(2:7), ~ (. - mean(., na.rm = T)) / sd(., na.rm = T),  .names = "z_{col}")) %>%
  # mutate_at(vars(z_gamma_cp, z_bbb_532), ~ . * -1) %>% # reverse these z-scores as they are inversely related to particle size
  # mutate(composite_z = rowSums(across(c(7:11))),
  #        biomass = ifelse(composite_z < 0, "Lower", "Higher")) %>% 
  #  left_join(read_csv(exploc_path), .) %>% 
  # left_join(., acs_n) %>% 
  # left_join(., cells_n) %>% 
  # arrange(composite_z) 
  mutate_at("stn", as_factor)
```

``` r
# make a train data sets that scaled and convert them to be a matrix beccause kohonen function accepts numeric matrix
som.train <- as.matrix(scale(som[,-1]))
```

``` r
# make stable sampling
RNGkind(sample.kind = "Rounding")
```

    ## Warning in RNGkind(sample.kind = "Rounding"): non-uniform 'Rounding' sampler
    ## used

``` r
# make a SOM grid
set.seed(100)
som.grid <- somgrid(xdim =2, ydim = 2, topo = "hexagonal") #user input of 2 clusters

# make a SOM model
set.seed(100)
som.model <- som(som.train, som.grid, rlen = 100, radius = 2.5, keep.data = TRUE,
                  dist.fcts = "euclidean")
```

``` r
set.seed(100)
fviz_nbclust(som.train, kmeans, method = "silhouette", k.max = 6 )
```

![](/Users/nicholasbaetge/github/oceprf_smokeonthewater/knitted/1_ILA2R_files/figure-gfm/find%20optimum%20number%20of%20clusters-1.png)<!-- -->

``` r
clust <- kmeans(som.train, 2)
```

``` r
svg(file = "~/github/oceprf_smokeonthewater/prod/figs/SOM_clusters.svg", width = 4, height = 4) 
plot(som.model, type = "codes", bgcol = rainbow(9)[clust$cluster], main = "Cluster Map")
add.cluster.boundaries(som.model, clust$cluster)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
clusters <- som.model$unit.classif
som_results <- som %>% 
  mutate(cluster = clusters) %>% 
  mutate(stn = as.character(stn),
         stn = as.numeric(stn)) %>% 
  left_join(read_csv(exploc_path), .) %>%
  left_join(., acs_n) %>%
  left_join(., cells_n) %>%
  arrange(cluster)
```

    ## Rows: 14 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): date, exp
    ## dbl (3): stn, lat, lon
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Joining with `by = join_by(stn)`
    ## Joining with `by = join_by(stn)`
    ## Joining with `by = join_by(stn)`

``` r
index <-  summary %>%
  select(
    stn,
    mean_chl_ap676lh,
    mean_gamma_cp,
    mean_bbb_532,
    mean_poc_cp_660,
    mean_nano_syn,
    mean_pico_syn,
  ) %>%
  rename_all( ~ stringr::str_replace(., "mean_", "")) %>%
  mutate(across(c(2:7), ~ (. - mean(., na.rm = T)) / sd(., na.rm = T),  .names = "z_{col}")) %>%
  mutate_at(vars(z_gamma_cp, z_bbb_532), ~ . * -1) %>% # reverse these z-scores as they are inversely related to particle size
  mutate(composite_z = rowSums(across(c(7:11)))) %>%
  select(stn, composite_z, contains("z")) %>% 
  left_join(som_results, .) %>% 
    arrange(composite_z)  %>% 
  mutate(cluster = ifelse(cluster == 1, "High biomass", "Low biomass"))
```

    ## Joining with `by = join_by(stn)`

## data for bco-dmo

``` r
matlab2POS_2 = function(x, timez = "UTC") {
    days = x - 719529   # 719529 = days from 1-1-0000 to 1-1-1970
    secs = days * 86400 # 86400 seconds in a day
    # Converting the secs value to a POSIXct value in the UTC time zone, then converts that to a time/date string that 
    # should lose the time zone, and then it performs a second as.POSIXct() conversion on the time/date string to get a POSIXct value in the user's specified timezone.
    return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
            tz = 'UTC'), format = '%Y-%m-%d %H:%M:%S', 
            tz = 'UTC', usetz = FALSE), tz = timez))
}

bco_dmo <-  acs %>%
  select(dt,
         ap_16,
         ap_29,
         ap_56,
         cp_16,
         cp_29,
         cp_56,
         poc,
         chl_ap676lh,
         gamma) %>%
  rename(
    ap_470 = ap_16,
    ap_532 = ap_29,
    ap_660 = ap_56,
    cp_470 = cp_16,
    cp_532 = cp_29,
    cp_660 = cp_56,
    poc_cp_660 = poc,
    gamma_cp = gamma
  ) %>%
  left_join(
    bb %>%
      select(dt, bbp_1:bbp_3) %>%
      rename(
        bbp_470 = bbp_1,
        bbp_532 = bbp_2,
        bbp_660 = bbp_3
      ),
    .
  ) %>%
 rename(m_dt = dt) %>%
  mutate(dt = matlab2POS_2(m_dt), .before = m_dt) %>%
  select(-m_dt) %>%
  mutate(stn = ifelse(between(dt, stn_times$start[stn_times$stn == 4], stn_times$end[stn_times$stn == 4]), 4, NA),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 5], stn_times$end[stn_times$stn == 5]), 5, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 6], stn_times$end[stn_times$stn == 6]), 6, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 7], stn_times$end[stn_times$stn == 7]), 7, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 8], stn_times$end[stn_times$stn == 8]), 8, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 9], stn_times$end[stn_times$stn == 9]), 9, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 10], stn_times$end[stn_times$stn == 10]), 10, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 11], stn_times$end[stn_times$stn == 11]), 11, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 12], stn_times$end[stn_times$stn == 12]), 12, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 13], stn_times$end[stn_times$stn == 13]), 13, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 14], stn_times$end[stn_times$stn == 14]), 14, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 15], stn_times$end[stn_times$stn == 15]), 15, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 16], stn_times$end[stn_times$stn == 16]), 16, stn),
         stn = ifelse(between(dt, stn_times$start[stn_times$stn == 17], stn_times$end[stn_times$stn == 17]), 17, stn)
    ) %>% 
  drop_na(stn) %>% 
  left_join(.,  read_csv(exploc_path)) %>% 
  select(stn, lat, lon, dt, everything(), -date, -exp) 
```

    ## Joining with `by = join_by(dt)`
    ## Rows: 14 Columns: 5
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (2): date, exp dbl (3): stn, lat, lon
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Joining with `by = join_by(stn)`

# save

``` r
write_csv(summary, prod_path)
write_csv(index, index_path)
write_csv(bco_dmo, bcodmo_path)
```
