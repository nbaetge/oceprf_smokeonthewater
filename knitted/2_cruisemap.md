Cruise map with satellite chl and ship sst
================
Nick Baetge
compiled most recently on 20 November, 2024

``` r
library(tidyverse)
library(ggOceanMaps)
library(lubridate)
library(ncdf4)
library(patchwork)
library(gt)
```

# import and tidy data

Level-2 VIIRS_SNPP files were downloaded from [Ocean Color
Web](https://oceancolor.gsfc.nasa.gov/cgi/browse.pl?sen=amod) and were
processed to Level-3 with a spatial binning of 1.1 km and temporal
binning coinciding with the cruise legs (7/28/23 - 8/9/23 (11 d),
8/12/23 - 8/19/23 (7 d))

``` r
index_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p1_insitu_index.csv"
meta_path <- "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p1_insitu_biooptics.csv"
metdata_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r0_shipmet.csv"
v1_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r2_viirs_scene1_20230728_20230731_chl.nc"
v2_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r2_viirs_scene2_20230801_20230809_chl.nc"
v3_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r2_viirs_scene3_20230812_20230819_chl.nc"
table_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/knitted/2_cruisemap_files/insitu_table.html"
suptable_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/knitted/2_cruisemap_files/zscore_table.html"
```

## viirs chl data (for ggplot raster)

``` r
# function for extracting data from nc files
nc_extract <- function(sat_path) {
  sat <- nc_open(sat_path)
  ret <- list()
  ret$lat <- ncvar_get(sat, "lat")
  ret$lon <- ncvar_get(sat, "lon")
  ret$chla <- ncvar_get(sat, "chlor_a")
  nc_close(sat)
  
  ret
}

# extract data from nc files
viirs1 <- nc_extract(v1_path)
viirs2 <- nc_extract(v2_path)
viirs3 <- nc_extract(v3_path)

# tidy data
melt_chl <- function(L) {
  dimnames(L$chla) <- list(lon = L$lon, lat = L$lat)
  ret <- reshape2::melt(L$chla, value.name = "chla")
}

viirs1_chla <-  melt_chl(viirs1)
viirs2_chla <-  melt_chl(viirs2)
viirs3_chla <-  melt_chl(viirs3)

#combine the first two scences (viirs 1 contains data from 7/28 - 7/30 and viirs 2 contains data from 8/1 - 8/9, together they make up the first leg of the cruise)
viirs1_2_chla <- viirs1_chla %>%
  rename(v1 = chla) %>%
  left_join(., viirs2_chla %>% rename(v2 = chla)) %>%
  pivot_longer(c(v1, v2), names_to = "scene", values_to = "val") %>%
  group_by(lon, lat) %>%
  mutate(chla = mean(val, na.rm = T)) %>%
  ungroup() %>%
  select(lon, lat, chla) %>%
  distinct()
```

## shipmet data (for ggplot ship tracks)

``` r
met <- read_csv(metdata_path) %>%
  select(Time, ZD, ST, LA, LO) %>%
  rename(sst_c = ST,
         lat = LA,
         lon = LO) %>%
  mutate(across(where(is.numeric), ~ na_if(., -99))) %>%
  mutate(dt = as_datetime(ZD)) %>%
  select(dt, lat, lon, sst_c) %>%
  mutate(date = as_date(dt), .after = dt) %>%
  mutate(viirs =
           case_when(between(
             date, ymd("2023-07-27"), ymd("2023-08-12")
           ) ~ 1,
           between(
             date, ymd("2023-08-13"), ymd("2023-08-21")
           ) ~ 2))
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

## experiment locations

``` r
index <- read_csv(index_path) %>%
  mutate(date = mdy(date)) %>%
  mutate(viirs = case_when(between(
    date, ymd("2023-07-27"), ymd("2023-08-12")
  ) ~ 1,
  between(
    date, ymd("2023-08-13"), ymd("2023-08-21")
  ) ~ 2)) 
```

# Plot data

``` r
region <-
  data.frame(lon = c(-116.7, -122.5, -122.5, -116.7),
             lat = c(32, 35.5, 35.5, 32))
chl_breaks <- c(0, 0.01, 0.1, 1, 10, 100)
```

``` r
custom.theme <- theme(
  plot.title = element_text(size = 35, face = "bold"),
  legend.position = "right", legend.key.height = unit(1.5, "cm"),
  legend.title = element_text(size = 22),
  legend.key.width = unit(0.7, "cm"),
  legend.key.size = unit(1, "cm"),
  legend.text = element_text(size = 18),
  legend.box = "vertical", 
  legend.box.just = "left",
  legend.justification = "left",
  axis.title = element_text(size = 35, face = "bold"),
  axis.text = element_text(size = 28),
  panel.spacing.x = unit(2, "cm"),
  strip.text.x = element_text(
    size = 35,
    color = "white",
    face = "bold"
  ),
  strip.text.y = element_text(
    size = 35,
    color = "white",
    face = "bold"
  ),
  strip.background = element_rect(
    color = "black",
    fill = alpha('black', 0.7),
    linewidth = 1.5,
    linetype = "solid"
  )
)
```

``` r
leg1 <- basemap(data = region, bathymetry = F) +
  geom_raster(data = viirs1_2_chla,
              aes(x = lon, y = lat, fill = chla),
              interpolate = TRUE, alpha = 0.85) +
  viridis::scale_fill_viridis(
    breaks = chl_breaks,
    labels = chl_breaks,
    trans = scales::pseudo_log_trans(sigma = 0.001),
    na.value = NA
  ) +
  # geom_path(
  #   data = met %>% filter(viirs == 1),
  #   aes(x = lon, y = lat),
  #   color = "white",
  #   linewidth = 4
  # )  +
  geom_path(
    data = met %>% filter(viirs == 1),
    aes(x = lon, y = lat, color = sst_c),
    linewidth = 3
  )  +
  viridis::scale_color_viridis(option = "H", na.value = NA) +
  labs(
    y = "Latitude",
    x = "Longitude",
    fill = expression(Chl ~ italic(a) ~ (µg ~ L ^ -1)),
    color = expression(SST ~ ("˚C"))
  ) +
  ggrepel::geom_label_repel(
    data = index %>% filter(viirs == 1) %>% filter(cluster == "High biomass"),
    aes(x = lon, y = lat, label = exp),
    alpha = 0.75,
    fontface = 'bold',
    fill = "#92EA33",
    color = "black",
    size = 10,
    box.padding = 0.4,
    min.segment.length = 1,
    segment.size = 1,
    xlim = c(-120, -116),
    ylim = c(35, 36)
  ) +
  ggrepel::geom_label_repel(
    data = index %>% filter(viirs == 1) %>% filter(cluster == "Low biomass"),
    aes(x = lon, y = lat, label = exp),
    alpha = 0.75,
    fontface = 'bold',
    fill = "#0947EA",
    color = "black",
    size = 10,
    box.padding = 0.4,
    min.segment.length = 1,
    segment.size = 1,
    xlim = c(-119, -116),
    ylim = c(34, 34.5)
  ) +
  geom_point(
    data = index %>% filter(viirs == 1) %>% filter(cluster == "High biomass") %>% drop_na(exp),
    aes(x = lon, y = lat),
    color = "#92EA33",
    size = 8,
    shape = 21,
    stroke = 3
  ) +
  geom_point(
    data = index %>% filter(viirs == 1) %>% filter(cluster == "Low biomass") %>% drop_na(exp),
    aes(x = lon, y = lat),
    color =  "#0947EA",
    size = 8,
    shape = 21,
    stroke = 3
  ) +
  theme_test() +
  custom.theme +
  guides(fill = guide_colorbar(order = 1), col = guide_colorbar(order = 2)) +
  ggtitle("07/28/23 - 08/09/23")
```

``` r
leg2 <- basemap(data = region, bathymetry = F) +
  geom_raster(data = viirs3_chla,
              aes(x = lon, y = lat, fill = chla),
              interpolate = TRUE, alpha = 0.85) +
  viridis::scale_fill_viridis(
    breaks = chl_breaks,
    labels = chl_breaks,
    trans = scales::pseudo_log_trans(sigma = 0.001),
    na.value = NA
  ) +
  # geom_path(
  #   data = met %>% filter(viirs == 2),
  #   aes(x = lon, y = lat),
  #   color = "white",
  #   linewidth = 6
  # )  +
  geom_path(
    data = met %>% filter(viirs == 2),
    aes(x = lon, y = lat, color = sst_c),
    linewidth = 3
  )  +
  viridis::scale_color_viridis(option = "H",  na.value = NA) +
  labs(
    y = "Latitude",
    x = "Longitude",
    fill = expression(Chl ~ italic(a) ~ (µg ~ L ^ -1)),
    color = expression(SST ~ ("˚C"))
  ) +
    
   ggrepel::geom_label_repel(
    data = index %>% filter(viirs == 2) %>% filter(cluster == "High biomass"),
    aes(x = lon, y = lat, label = exp),
    alpha = 0.75,
    fontface = 'bold',
    fill = "#92EA33",
    color = "black",
    size = 10,
    box.padding = 0.8,
    min.segment.length = 1,
    segment.size = 1,
     xlim = c(-120.5, -114),
    ylim = c(5, 36)
  ) +
  ggrepel::geom_label_repel(
    data = index %>% filter(viirs == 2) %>% filter(cluster == "Low biomass") ,
    aes(x = lon, y = lat, label = exp),
    alpha = 0.75,
    fontface = 'bold',
    fill =  "#0947EA",
    color = "black",
    size = 10,
    box.padding = 0.8,
    min.segment.length = 1,
    segment.size = 1,
    xlim = c(-121, -114),
    ylim = c(33.5, 35),
    seed = 1
  ) +
  
   geom_point(
    data = index %>% filter(viirs == 2) %>% filter(cluster == "High biomass") %>% drop_na(exp),
    aes(x = lon, y = lat),
    color =  "#92EA33",
    size = 8,
    shape = 21,
    stroke = 3
  ) +
  geom_point(
    data = index %>% filter(viirs == 2) %>% filter(cluster == "Low biomass")  %>% drop_na(exp),
    aes(x = lon, y = lat),
    color =  "#0947EA",
    size = 8,
    shape = 21,
    stroke = 3
  ) +
  theme_test() +
  custom.theme +
  guides(fill = guide_colorbar(order = 1), col = guide_colorbar(order = 2)) +
  ggtitle("08/12/23 - 08/19/23")
```

``` r
maps <- (leg1) + (leg2 + guides(color = "none", fill = "none"))  + plot_layout(guides = "collect")

ggsave(maps, file = "~/github/oceprf_smokeonthewater/prod/figs/maps.svg", width = 30, height = 14)
```

# Tables

``` r
table_data <- left_join(read_csv(index_path), read_csv(meta_path) %>% select(stn, contains("sd"))) %>% 
  select(exp, stn, date, lat, lon, everything(), -stn) %>% 
  mutate_at(vars(contains(c("lat", "lon", "z", "chl", "gamma", "poc", "nano_syn", "pico_syn"))), round, 2) %>% 
  mutate_at(vars(contains("bbb")), round, 4) %>% 
  arrange(composite_z) %>% 
  select(exp, date, lat, lon, cluster, acs_n, chl_ap676lh, poc_cp_660, gamma_cp, bbb_532, cells_n, everything(), -contains("poc_chl"), -contains("z")) 
```

    ## Rows: 14 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (3): date, exp, cluster
    ## dbl (18): stn, lat, lon, chl_ap676lh, gamma_cp, bbb_532, poc_cp_660, nano_sy...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 14 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (15): stn, mean_chl_ap676lh, sd_chl_ap676lh, mean_gamma_cp, sd_gamma_cp,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Joining with `by = join_by(stn)`

``` r
gt_tbl <- gt(table_data ) 

table <- 
  gt_tbl |>
  # tab_header(
  #   title = md("**Pre-existing in situ biomass conditions**"),
  # ) |>
  cols_label(
    exp = html("Experimental site"),
    date = html("Date"),
    lat = html("Latitude, <br>&deg;N"),
    lon = html("Longitude, <br>&deg;W"),
    cluster = md("**SOM cluster**"),
    chl_ap676lh = html("Chl<sub>a<sub>p</sub>(676)lh<br>(mg m<sup>-3</sup>)"),
    gamma_cp = html("&gamma;"),
    bbb_532 = html("b<sub>bp</sub>/b<sub>p</sub>(532)"),
    poc_cp_660 = html("POC<sub>c<sub>p</sub>(660)<br>(mg m<sup>-3</sup>)"),
    acs_n = html("n"),
    nano_syn = html("nanoeukaryotes:<br><i>Synechococcus</i>"),
    pico_syn = html("picoeukaryotes:<br><i>Synechococcus</i>"),
    cells_n = html("n")
    
  ) |>
  cols_merge_uncert(
    col_val = chl_ap676lh,
    col_uncert = sd_chl_ap676lh
  ) |>
  cols_merge_uncert(
    col_val = gamma_cp,
    col_uncert = sd_gamma_cp
  ) |>
  cols_merge_uncert(
    col_val = bbb_532,
    col_uncert = sd_bbb_532
  )|>
  cols_merge_uncert(
    col_val = poc_cp_660,
    col_uncert = sd_poc_cp_660
  )|>
   cols_merge_uncert(
    col_val = nano_syn,
    col_uncert = sd_nano_syn
  )|>
  cols_merge_uncert(
    col_val = pico_syn,
    col_uncert = sd_pico_syn
  )|>
  tab_spanner(
    label = html("Inline bio-optics"),
    columns = c(chl_ap676lh, poc_cp_660, gamma_cp, bbb_532, acs_n)
  ) |>
  tab_spanner(
    label = html("Flow cytometry"),
    columns = c(nano_syn, pico_syn, cells_n )
  ) 

gtsave(table, table_path)
```

``` r
suptable_data <- left_join(read_csv(index_path), read_csv(meta_path) %>% select(stn, contains("sd"))) %>% 
  select(exp, date, lat, lon, composite_z, contains("z"))  %>% 
  mutate_at(vars(contains(c("lat", "lon", "z"))), round, 2) %>% 
  arrange(composite_z)
```

    ## Rows: 14 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (3): date, exp, cluster
    ## dbl (18): stn, lat, lon, chl_ap676lh, gamma_cp, bbb_532, poc_cp_660, nano_sy...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 14 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (15): stn, mean_chl_ap676lh, sd_chl_ap676lh, mean_gamma_cp, sd_gamma_cp,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Joining with `by = join_by(stn)`

``` r
gt_sup <- gt(suptable_data ) 

suptable <- 
  gt_sup |>
  tab_header(
    title = md("**Table S1.** In situ biomass conditions from composite z-scores"),
  ) |>
  cols_label(
    exp = html("Experimental site"),
    date = html("Date"),
    lat = html("Latitude, <br>&deg;N"),
    lon = html("Longitude, <br>&deg;W"),
    composite_z = md("**Composite z-score**"),
    z_chl_ap676lh = html("Chl<sub>a<sub>p</sub>(676)lh"),
    z_gamma_cp = html("&gamma;"),
    z_bbb_532 = html("b<sub>bp</sub>/b<sub>p</sub>(532)"),
    z_poc_cp_660 = html("POC<sub>c<sub>p</sub>(660)"),
    z_nano_syn = html("nanoeukaryotes:<br><i>Synechococcus</i>"),
    z_pico_syn = html("picoeukaryotes:<br><i>Synechococcus</i>")
    
  ) |>
  tab_spanner(
    label = html("Inline bio-optics z-scores"),
    columns = c(z_chl_ap676lh, z_poc_cp_660, z_gamma_cp, z_bbb_532)
  ) |>
  tab_spanner(
    label = html("Flow cytometry z-scores"),
    columns = c(z_nano_syn, z_pico_syn)
  ) |>
  tab_source_note(source_note = md(
    "Z-scores for each variable were calculated as the experimental site estimate minus the mean of estimates across all experimental sites, then divided by the standard deviation of the estimates across all experimental sites. Z-scores for γ and bbp:bp were reversed as both metrics are inversely related to particle size. A composite z-score, or biomass index, is calculated as the sum of the z-scores for each site."
  )) 

gtsave(suptable, suptable_path)
```
