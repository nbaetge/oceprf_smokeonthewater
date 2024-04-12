Bacteria and DOM experiments
================
Nick Baetge
compiled most recently on 12 April, 2024

``` r
library(tidyverse)
library(hms)
library(lubridate)
library(purrr)
library(patchwork)
library(janitor)
```

``` r
prod_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p5_bacteria_dom_exps.csv"
```

# Bacterial cell data

## Import and tidy data

``` r
fcm_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p4_fcm_plate_data"
filter_vols <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r5_filter_vols.xlsx"
time_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r5_bact_dom_sampletimes.csv"
```

``` r
dt <- read_csv(time_path) %>%
  mutate_at(vars(datetime_utc), ymd_hms) %>%
  rename(dt = 3) %>%
  group_by(exp) %>%
  mutate(
    interv = interval(first(dt), dt),
    dur = as.duration(interv),
    days = round(time_length(dur, unit = "day"), 2)
  ) %>%
  select(exp, tp, dt, days)
```

import data and clean based on run notes:

- B2 run 1 (and dil): timepoint 1 and 2 are switched
- B3_14 run 1 dil: B14-2-4-0 didn’t add sample so data is in well H12
- B16 run 2 dil: B16-V2-4-4 didn’t add sample so data is in well G1
- B3_14 run 2 B14_V1_3 was a 2 ml sample instead of a 1 ml sample, also
  extra 1.2 vial (could be 1.2 GF75)
- Filtrates run 2: missing vial for B3 1.2 GF75

``` r
fcm <-
  list.files(path = fcm_path,
             pattern = "*.csv",
             full.names = T) %>%
  map_df( ~ read_csv(., col_types = cols(.default = "c"))) %>%
  mutate_at(vars(tp, dil, cells_ml), as.numeric) %>%
  mutate(tp = ifelse(
    name %in% c(
      "B2_dil-1.fcs",
      "B2-7.fcs",
      "B2_dil-2.fcs",
      "B2-8.fcs",
      "B2_dil-3.fcs",
      "B2-9.fcs",
      "B2_dil-4.fcs",
      "B2-10.fcs"
    ),
    2,
    tp
  ),
  tp = ifelse(
    name %in% c("B2-11.fcs", "B2-12.fcs", "B2-13.fcs", "B2-14.fcs"),
    1,
    tp
  )) %>%   # sample switch
  mutate(btl = ifelse(name == "B3_14-1.fcs", "1.2GF75", btl)) %>%    # sample switch
  filter(!name %in% c(
    "B3_14_dil-33.fcs",
    "B15-2.fcs",
    "B16_r2-25.fcs",
    "B16_r2-25.fcs"
  )) %>%    # samples had 0 events
  mutate(
    trt = case_when(
      btl %in% c("1", "V1", "1.2T", "1.2DNA", "1.2GF75", "1DNA", "1GF75") ~ 1,
      btl %in% c("2", "V2", "2DNA", "2GF75") ~ 2,
      btl %in% c("3", "V3", "3DNA", "3GF75") ~ 3,
      btl %in% c("4", "V4", "4DNA", "4GF75") ~ 4
    ),
    type = ifelse(btl %in% c("1", "2", "3", "4"), "5L Bottle", "40 mL Vial"),
    desc = case_when(
      trt == 1 ~ "Control",
      trt == 2 ~ "Thomas Fire Ash",
      trt == 3 ~ "Low Temp. Ash",
      trt == 4 ~ "High Temp. Ash"
    ),
    .after = btl
  ) %>%
  arrange(exp, trt, tp) 
```

## Subtract sybr stain

``` r
sybr <- fcm %>%
  filter(exp == "TE+SYBR") %>%
  mutate(med_sybr = round(median(cells_ml))) %>%
  select(med_sybr) %>%
  distinct() %>%
  as.numeric()
```

``` r
fcm_corr <- fcm %>%
  select(name, exp, btl, trt, desc, type, tp, dil, cells_ml) %>%
  arrange(exp, btl, tp) %>%
  filter(!exp %in% c("TE", "TE+SYBR")) %>%
  mutate(cells_ml = cells_ml - sybr,
         cells_ml = ifelse(cells_ml < 0, 0, cells_ml)) 
```

## Assess vial v bottle cell dynamics

``` r
reg_data <- fcm %>%
  filter(btl %in% c("1", "2", "3", "4", "V1", "V2", "V3", "V4")) %>%
  group_by(exp, trt, btl, type, tp) %>%
  mutate(mean_cells_ml = round(mean(cells_ml))) %>%
  ungroup() %>%
  select(exp, trt, desc, type, tp, dil, mean_cells_ml) %>%
  distinct() %>%
  pivot_wider(
    id_cols = c(exp, trt, desc, tp, dil),
    names_from = type,
    values_from = mean_cells_ml
  ) %>%
  drop_na(7)

reg <-
  lmodel2::lmodel2(reg_data$`40 mL Vial` ~ reg_data$`5L Bottle`, nperm = 99)

reg
```

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2::lmodel2(formula = reg_data$`40 mL Vial` ~ reg_data$`5L
    ## Bottle`, nperm = 99)
    ## 
    ## n = 113   r = 0.9833924   r-square = 0.9670606 
    ## Parametric P-values:   2-tailed = 4.155702e-84    1-tailed = 2.077851e-84 
    ## Angle between the two OLS regression lines = 0.9520367 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept     Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  45497.37 0.8677869        40.95102              0.01
    ## 2     MA  35837.72 0.8805902        41.36683              0.01
    ## 3    SMA  34440.51 0.8824421        41.42654                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS       21756.78        69237.96  0.8376644   0.8979094
    ## 2     MA       12417.52        58558.30  0.8504753   0.9116323
    ## 3    SMA       11326.38        56779.10  0.8528336   0.9130786
    ## 
    ## Eigenvalues: 92447880184 761915031 
    ## 
    ## H statistic used for computing C.I. of MA: 0.0002964096

The cell dynamics occurring in the bottle and the vials track well with
each other. we can use averages between the two for subsequent
measurements. so cells per ml will reflect the average of two subsamples
from the vials and two subsamples from the bottle (n = 4)

## Extract experimental data

``` r
# apply cell counts from bottles at T0 as vial values for T0
fcm_40mlt0 <- fcm %>%
  filter(btl %in% c("1", "2", "3", "4")) %>%
  group_by(exp, trt, btl, type, tp) %>%
  mutate(mean_cells_ml = round(mean(cells_ml)),
         sd_cells_ml = round(sd(cells_ml))) %>%
  ungroup() %>%
  filter(tp == 0) %>%
  mutate(type = "40 mL Vial",
         btl = paste("V", btl, sep = ""))

#extract experimental data and merge with datetime data
fcm_exp <- fcm %>%
  filter(btl %in% c("1", "2", "3", "4", "V1", "V2", "V3", "V4")) %>%
  group_by(exp, trt, btl, type, tp) %>%
  mutate(mean_cells_ml = round(mean(cells_ml)),
         sd_cells_ml = round(sd(cells_ml))) %>%
  ungroup() %>%
  bind_rows(., fcm_40mlt0) %>%
  select(-name) %>%
  arrange(exp, btl, tp) %>%
  left_join(., dt) %>%
  select(exp:tp, dt, days, everything()) %>%
  mutate(hours = days * 24, .after = days)

#calculate mean data
fcm_exp_means <- fcm_exp %>%
  select(-c(type, dil, btl, mean_cells_ml, sd_cells_ml)) %>%
  arrange(exp, trt, tp) %>%
  group_by(exp, trt, tp) %>%
  mutate(mean_cells_ml = round(mean(cells_ml)),
         sd_cells_ml = round(sd(cells_ml))) %>%
  select(-cells_ml) %>%
  distinct() %>% 
  ungroup() 
```

## Estimate growth rates (d-1) in exponential phase

``` r
fcm_gr <- fcm_exp %>%
  filter(
    exp == "B1" & tp <= 3 |
      exp == "B2" & tp <= 2 |
      exp == "B3" & tp <= 2 |
      exp == "B14" & tp <= 2 |
      exp == "B15" &  tp %in% c(1, 2, 3) |
      exp == "B16" & tp %in% c(1, 2, 3)
  ) %>%
  group_by(exp, btl, trt) %>%
  mutate(gr_d = (log(last(mean_cells_ml)) - log(first(mean_cells_ml))) / (last(days) - first(days))) %>%
  select(exp:desc, gr_d) %>%
  distinct() %>%
  group_by(exp, trt, desc) %>%
  summarise_at(vars(gr_d), list(mean = mean, sd = sd)) %>%
  ungroup() %>%
  rename(mean_gr_d = mean,
         sd_gr_d = sd)

p5_fcm <- fcm_gr  %>%
  left_join(fcm_exp_means, .) %>%
  select(exp,
         desc,
         tp,
         dt,
         days,
         hours,
         mean_cells_ml,
         sd_cells_ml,
         mean_gr_d,
         sd_gr_d) %>%
  distinct()
```

## Estimate filter retention for bacterial organic carbon and dna samples

``` r
fcm_gf75_dna <- fcm %>%
  filter(tp %in% c(0, 2, 5)) %>%
  mutate(cells_ml = ifelse(
    name %in% c(
      "Filtrates-21.fcs",
      "Filtrates_r2-21.fcs",
      "Filtrates_r2-19.fcs",
      "Filtrates-19.fcs"
    ),
    NA,
    cells_ml
  )) %>%  # suspect samples. these gf75 filtrates had more cells than any other sample in the dataset.
  arrange(exp, tp, btl) %>%
  group_by(exp, btl, tp) %>%
  mutate(mean_cells_ml = round(mean(cells_ml))) %>%
  ungroup() %>%
  select(-c(name, dil, cells_ml)) %>%
  distinct() %>%
  filter(!btl %in% c(1, 2, 3, 4) | !tp == 0) %>%
  group_by(exp, trt, tp) %>%
  mutate(
    total = ifelse(
      str_detect(btl, "DNA") |
        str_detect(btl, "GF75"),
      NA,
      mean_cells_ml
    ),
    gf75_filtrate = ifelse(str_detect(btl, "GF75"), mean_cells_ml, NA),
    dna_filtrate = ifelse(str_detect(btl, "DNA"), mean_cells_ml, NA)
  ) %>%
  fill(c(total:dna_filtrate), .direction = "updown") %>%
  ungroup() %>%
  select(exp, trt, desc, tp, total:dna_filtrate) %>%
  distinct() %>%
  group_by(exp, trt, tp) %>%
  mutate(mean_total = round(mean(total, na.rm = T)), .after = tp) %>%
  select(-total) %>%
  distinct() %>%
  rename(total = mean_total) %>%
  left_join(readxl::read_xlsx(filter_vols, sheet = 1), .) %>%
  mutate(
    gf75_total = total * gf75_vol,
    dna_total = total * dna_vol,
    .before = gf75_filtrate
  ) %>% #total cells in 1L
  mutate_at(vars(dna_filtrate), ~ . * dna_vol) %>%
  mutate_at(vars(gf75_filtrate), ~ . * gf75_vol) %>%
  select(-total) %>%
  mutate(
    gf75_filter = gf75_total - gf75_filtrate,
    dna_filter = dna_total - dna_filtrate,
    gf75_percent_retention = round(gf75_filter / gf75_total * 100),
    dna_percent_retention = round(dna_filter / dna_total * 100)
  )
```

``` r
summary(fcm_gf75_dna)
```

    ##      exp                 trt            desc                 tp       
    ##  Length:38          Min.   :1.000   Length:38          Min.   :0.000  
    ##  Class :character   1st Qu.:1.000   Class :character   1st Qu.:2.000  
    ##  Mode  :character   Median :2.000   Mode  :character   Median :2.000  
    ##                     Mean   :1.842                      Mean   :2.947  
    ##                     3rd Qu.:2.000                      3rd Qu.:5.000  
    ##                     Max.   :4.000                      Max.   :5.000  
    ##                                                                       
    ##     gf75_vol       dna_vol         gf75_total          dna_total        
    ##  Min.   :1000   Min.   : 500.0   Min.   :2.910e+08   Min.   :1.455e+08  
    ##  1st Qu.:1000   1st Qu.: 500.0   1st Qu.:5.542e+08   1st Qu.:3.183e+08  
    ##  Median :1000   Median : 500.0   Median :6.937e+08   Median :3.794e+08  
    ##  Mean   :1019   Mean   : 632.9   Mean   :6.736e+08   Mean   :4.007e+08  
    ##  3rd Qu.:1000   3rd Qu.: 500.0   3rd Qu.:8.084e+08   3rd Qu.:4.285e+08  
    ##  Max.   :1500   Max.   :1500.0   Max.   :1.019e+09   Max.   :1.034e+09  
    ##                                                                         
    ##  gf75_filtrate        dna_filtrate        gf75_filter       
    ##  Min.   : 46126000   Min.   :  6418000   Min.   :197302000  
    ##  1st Qu.:104019250   1st Qu.: 10254250   1st Qu.:453859500  
    ##  Median :121111000   Median : 15845500   Median :561488000  
    ##  Mean   :135120104   Mean   : 32787957   Mean   :556954385  
    ##  3rd Qu.:155165000   3rd Qu.: 34574588   3rd Qu.:661688500  
    ##  Max.   :293678000   Max.   :256902000   Max.   :910566950  
    ##  NA's   :2                               NA's   :2          
    ##    dna_filter        gf75_percent_retention dna_percent_retention
    ##  Min.   :108963000   Min.   :61.00          Min.   :73.00        
    ##  1st Qu.:278489250   1st Qu.:76.50          1st Qu.:91.50        
    ##  Median :358163000   Median :82.00          Median :95.00        
    ##  Mean   :367894593   Mean   :80.22          Mean   :92.55        
    ##  3rd Qu.:403498500   3rd Qu.:84.25          3rd Qu.:97.00        
    ##  Max.   :831114000   Max.   :90.00          Max.   :98.00        
    ##                      NA's   :2

``` r
sd(fcm_gf75_dna$gf75_percent_retention, na.rm = T)
```

    ## [1] 7.593648

``` r
sd(fcm_gf75_dna$dna_percent_retention)
```

    ## [1] 7.127023

# Bacterial organic carbon data

## Import and tidy data

``` r
chn_path <- "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r5_chn.csv"
```

``` r
p5_boc <- read_csv(chn_path) %>% 
  left_join(., readxl::read_xlsx(filter_vols, sheet = 2)) %>% 
  select(1, 2, 4, 5, 7) %>% 
  rename(rundate = 1,
         exp_sample_tp_rep = 2, 
         c_ug = 3, 
         n_ug = 4) %>% 
  separate(., 2, into = c("exp", "sample", "tp", "rep"), sep = "-") %>% 
  separate(., 3, into = c("sample", "filter"), sep = -1) %>% 
  group_by(exp, sample, tp, rep) %>% 
   mutate(raw_c = round(sum(c_ug),1),
         raw_n = round(sum(n_ug),1)) %>% 
  select(-c(filter, c_ug, n_ug)) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(c_blank = round(mean(raw_c[exp == "BLK"], na.rm = T), 1),
         n_blank = round(mean(raw_n[exp == "BLK"], na.rm = T), 1),
         c_blank_sd = round(sd(raw_c[exp == "BLK"], na.rm = T), 1),
         n_blank_sd = round(sd(raw_n[exp == "BLK"], na.rm = T), 1)) %>% 
  filter(!exp == "BLK") %>% 
  mutate(c_ug = raw_c - c_blank,
         n_ug = raw_n - n_blank) %>% 
  filter(!c_ug < 0) %>% 
  mutate(c_mg = c_ug / 10^3,
         n_mg = n_ug / 10^3,
         c_umol = c_mg * 10^3 / 12,
         n_umol = n_mg * 10^3 / 14,
         vol_m3 = vol / 10^6,
         c_mg_m3 = round(c_mg / vol_m3, 1),
         n_mg_m3 = round(n_mg / vol_m3, 1),
         c_umol_l = c_mg_m3/12000,
         n_umol_l = n_mg_m3/14000,
         cn =  round(c_umol_l/n_umol_l, 1)) %>% 
  select(exp, sample, tp, c_mg:cn) %>% 
  rename(trt = sample) %>% 
  group_by(exp, trt, tp) %>% 
  summarize(across(
    .cols = c_mg:cn,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  )) %>% 
  ungroup() %>% 
  mutate_at(vars(tp), as.numeric) %>% 
  mutate(trt = ifelse(trt == 1.2, 1, trt)) %>% 
  group_by(exp, trt) %>% 
  group_modify(~ add_row(., tp = 0 )) %>% 
  arrange(exp,  tp) %>% 
  group_by(exp, tp) %>% 
  fill(c(3:23), .direction = "down") %>% 
  ungroup() %>% 
  distinct() %>% 
  arrange(exp, trt, tp) %>% 
  mutate_at(vars(trt), as.numeric) %>% 
  left_join(., fcm_gf75_dna %>% select(exp, trt, desc, tp, gf75_filter)) %>% 
  group_by(exp, tp) %>% 
  fill(gf75_filter, .direction = "downup") %>% 
  ungroup() %>% 
  rename(gf75_cells = gf75_filter) %>% 
  select(exp, trt, desc, tp, gf75_cells, everything()) %>% 
  mutate(fgC_cell = (mean_c_umol * 12 * 10^9)/gf75_cells,
         fgN_cell = (mean_n_umol * 14 * 10^9)/gf75_cells ) %>% 
  mutate(desc = case_when(trt == 1 ~ "Control",
                          trt == 2 ~ "Thomas Fire Ash",
                          trt == 3 ~ "Low Temp. Ash",
                          trt == 4 ~ "High Temp. Ash")) %>% 
  select(exp, trt, desc, tp, fgC_cell, fgN_cell) %>% 
  left_join(p5_fcm,.) %>%
  mutate(boc = (mean_cells_ml * fgC_cell) / (12 * 10^6),
         bon = (mean_cells_ml * fgN_cell) / (14 * 10^6)) %>% 
  select(-trt)  %>% 
  rename(trt = desc)
```

c_blank = 15.3 +/- 5.1, n = 9 n_blank = 7.1 +/- 1, n = 9

# DOC data

## Import and tidy data

``` r
doc_path <- "//Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r5_doc.xlsx"
```

``` r
doc <- readxl::read_xlsx(doc_path, skip = 27) %>%
  drop_na() %>%
  select(1, 3) %>%
  rename(id = 1,
         toc = 2) %>%
  separate(id, into = c("exp", "btl", "tp"), sep = "-") %>%
  group_by(exp, btl, tp) %>%
  mutate(mean_toc = mean(toc),
         sd_toc = sd(toc),
         tp = as.numeric(tp)) %>%
  ungroup() %>%
  select(-toc) %>%
  distinct() %>%
  mutate(
    trt =  case_when(
      btl == 1 ~ "Control",
      btl == 2 ~ "Thomas Fire Ash",
      btl == 3 ~ "Low Temp. Ash",
      btl == 4 ~ "High Temp. Ash"
    ),
    .after = btl
  ) %>%
  select(-btl) %>%
  left_join(p5_boc, .) %>%
  mutate(doc = mean_toc - boc)

#calculate boc rate of change as slope
boc_rate <- doc %>% 
  select(exp, trt, tp, days, boc) %>%
  drop_na(boc) %>%
  group_by(exp, trt) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  group_split(group_id) %>%
  map( ~ lm(boc ~ days, data = .x)) %>%
  map_df(broom::tidy, .id = "group_id") %>%
  filter(term == "days") %>%
  select(group_id, estimate) %>%
  rename(boc_rate = estimate) %>%
  mutate_at(vars(group_id), as.numeric) %>%
  ungroup() %>% 
  left_join(
    .,
    doc %>%
      select(exp, trt) %>% distinct() %>%
      group_by(exp, trt) %>%
      mutate(group_id = cur_group_id()) %>%
      ungroup()
  ) %>% 
  select(-group_id)

#calculate doc rate of change as slope
doc_rate <- doc %>%
  select(exp, trt, tp, days, doc) %>%
  drop_na(doc) %>%
  group_by(exp, trt) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  group_split(group_id) %>%
  map( ~ lm(doc ~ days, data = .x)) %>%
  map_df(broom::tidy, .id = "group_id") %>%
  filter(term == "days") %>%
  select(group_id, estimate) %>%
  rename(doc_rate = estimate) %>%
  mutate_at(vars(group_id), as.numeric) %>%
  mutate_at(vars(doc_rate), abs) %>%
  ungroup() %>% 
  left_join(
    .,
    doc %>%
      select(exp, trt) %>% distinct() %>%
      group_by(exp, trt) %>%
      mutate(group_id = cur_group_id()) %>%
      ungroup()
  ) %>% 
  select(-group_id)

#calculate bges and doc bioavailability
p5_bge_bioav <- doc %>% 
  left_join(., boc_rate) %>% 
  left_join(., doc_rate) %>% 
  group_by(exp, trt) %>% 
  mutate(doc_rem = ifelse(is.na(doc[tp == 5]), first(doc) - doc[tp == 2], first(doc) - last(doc) ),
         doc_bdl = ifelse(doc_rem <= 1.4, T, F),
         bge = ifelse(doc_bdl == F, boc_rate/doc_rate, NA)) %>% 
  ungroup() %>% 
  group_by(exp, tp) %>% 
  mutate(control_doc = doc[trt == "Control"]) %>% 
  ungroup() %>% 
  mutate(amend_doc = doc - control_doc) %>% 
  group_by(exp, trt) %>% 
  mutate(amend_doc_bioav = ifelse(trt != "Control", (amend_doc[tp == 0] - amend_doc[tp == 5]) / amend_doc[tp == 0]* 100, NA )) %>% 
  ungroup() %>% 
  mutate(exp = paste("P", exp, sep = ""))
```

# Merge data with contextual in situ data

``` r
context_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p1_insitu_index.csv"

p5_prod <- read_csv(context_path) %>% 
  select(stn, exp, lat, lon, biomass, composite_z) %>% 
  left_join(., p5_bge_bioav) %>% 
  arrange(biomass, composite_z) %>% 
  drop_na(trt) %>% 
  mutate_at(vars(composite_z), round, 2) %>% 
  arrange(desc(composite_z), .by_group = TRUE) 
levels(p5_prod$composite_z)
```

    ## NULL

# Save data

``` r
write_csv(p5_prod, prod_path)
```

# Plot data

``` r
custom.theme <- theme(
  legend.position = "top",
  legend.title = element_text(size = 35),
  legend.key.size = unit(3, "cm"),
  legend.key.spacing.x = unit(0.75, "cm"),
  legend.text = element_text(size = 38),
  legend.box = "vertical",
  legend.margin = margin(),
  axis.title = element_text(size = 35, face = "bold"),
  axis.text = element_text(size = 28),
  panel.spacing.x = unit(0.5, "cm"),
  panel.spacing.y = unit(1, "cm"),
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

pal4 = c("#2C1B21", "#6aa42c", "#5990e7", "#b95f83")
pal3 = c("#6aa42c", "#5990e7", "#b95f83")
plot_levels = c(
  "Control",
  "Thomas Fire Ash",
  "Low Temp. Ash",
  "High Temp. Ash",
  "DOC",
  "bold(Leachate)",
  "bold(Amendment)",
  "bold(Lower~biomass~index)",
  "bold(Higher~biomass~index)",
  "bold(TOC~(µmol~C~L^-1))",
  "bold(Bacterial~C~(µmol~C~L^-1))",
  "bold(Cells~(L^-1))",
  "bold(Growth~rate~(d^-1))",
  "bold(BGE~(dimensionless))",
  "bold(+DOC~bioav.~('%'))",
  "PB16",
  "PB15",
  "PB1",
  "PB14",
  "PB2",
  "PB3"
) 
```

``` r
fig4_data <- p5_prod %>%
  select(exp, trt, biomass, composite_z, days, mean_toc, boc) %>%
  rename(toc = mean_toc) %>%
  pivot_longer(toc:boc, names_to = "var", values_to = "val") %>%
  mutate(
    var = case_when(
      var == "toc" ~ "bold(TOC~(µmol~C~L^-1))",
      var == "boc" ~ "bold(Bacterial~C~(µmol~C~L^-1))"
    )
  ) %>%
  left_join(
    .,
    p5_prod %>%
      select(exp, trt, biomass, days,  composite_z, sd_toc) %>%
      mutate(var = "bold(TOC~(µmol~C~L^-1))") %>%
      rename(sd = sd_toc)
  ) %>%
  bind_rows(
    .,
    p5_prod %>%
      select(exp, trt, biomass, days,  composite_z, mean_cells_ml, sd_cells_ml) %>%
      mutate_at(vars(mean_cells_ml, sd_cells_ml), ~ . * 10 ^ 3) %>%
      mutate(var = "bold(Cells~(L^-1))") %>%
      rename(val = mean_cells_ml,
             sd = sd_cells_ml)
  ) %>%
  arrange(exp, trt, days, var) %>%
  mutate(
    biomass = ifelse(
      biomass == "Lower",
      "bold(Lower~biomass~index)",
      "bold(Higher~biomass~index)"
    )
  ) 
```

``` r
curves <- ggplot(fig4_data[!is.na(fig4_data$val), ],
                 aes(
                   x = days,
                   y = val,
                   group = interaction(exp, trt, var),
                   color = factor(trt, levels = plot_levels),
                   fill = factor(trt, levels = plot_levels)
                 )) +
  geom_errorbar(
    aes(ymin = val - sd, ymax = val + sd),
    width = 0.2,
    linewidth = 0.6,
    alpha = 0.9
  ) +
  geom_line(linewidth = 1.1, alpha = 0.8) +
  geom_point(
    color = "black",
    size = 8,
    shape = 21,
    alpha = 0.8
  ) +
  scale_color_manual(values = pal4) +
  scale_fill_manual(values = pal4) +
  ggh4x::facet_nested(
    factor(var, levels = plot_levels) ~  factor(biomass, levels = plot_levels) + factor(
      composite_z,
      levels = c("-5.63", "-3.29", "-0.63", "0.39", "1.27", "1.53")
    ),
    scales = "free_y",
    labeller = label_parsed
  ) +
  labs(y = "", 
    x = expression(Day),
    fill = "",
    color = "",
    linetype = ""
  ) +
  guides(color = "none") +
  theme_linedraw() +
  custom.theme 
```

``` r
metric_data <- p5_prod %>%
  select(exp,
         trt,
         biomass,
         composite_z,
         mean_gr_d,
         bge,
         amend_doc_bioav) %>%
  distinct() %>%
  pivot_longer(mean_gr_d:amend_doc_bioav,
               names_to = "var",
               values_to = "val") %>%
  mutate(
    var = case_when(
      var == "mean_gr_d" ~ "bold(Growth~rate~(d^-1))",
      var == "bge" ~ "bold(BGE~(dimensionless))",
      var == "amend_doc_bioav" ~ "bold(+DOC~bioav.~('%'))",
    )
  ) %>%
  left_join(
    .,
    p5_prod %>%
      select(exp, trt, biomass,  composite_z, sd_gr_d) %>%
      distinct() %>%
      mutate(var = "bold(Growth~rate~(d^-1))") %>%
      rename(sd = sd_gr_d)
  ) %>%
  mutate(biomass = ifelse(
    biomass == "Lower",
    "bold(Lower~biomass~index)",
    "bold(Higher~biomass~index)"
  ))
```

``` r
metrics <- ggplot(metric_data,
                  aes(
                    x = factor(
                      composite_z,
                      levels = c("-5.63", "-3.29", "-0.63", "0.39", "1.27", "1.53")
                    ),
                    y = val,
                    fill = factor(trt, levels = plot_levels),
                    color = factor(trt, levels = plot_levels)
                  )) +
  geom_errorbar(aes(ymin = val - sd,
                    ymax = val + sd),
                width = .1,
                position = position_dodge(.2)) +
  # geom_linerange(
  #   aes(ymin = 0, ymax = val),
  #   position = position_dodge(width = 0.4),
  #   linewidth = 1,
  #   alpha = 0.8
  # ) +
  geom_point(position = position_dodge(width = 0.2),
             size = 10,
             alpha = 0.8) +
  ggh4x::facet_nested_wrap(
    ~ factor(var, levels = plot_levels),
    labeller = label_parsed,
    scales = "free",
    strip.position = "right",
    nrow = 3
  ) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  labs(y = expression(""),
       x = expression("Biomass index"),
       color = "") +
  theme_linedraw() +
  custom.theme +
  guides(fill = "none") +
  theme(panel.spacing.y = unit(0.5, "cm"))
```

``` r
(
  curves
) /(metrics  + guides(color = "none", fill = "none")) + plot_layout(heights = c(0.8, 0.7)) + plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 28))
```

![](/Users/nicholasbaetge/github/oceprf_smokeonthewater/knitted/5_bacteria_dom_exps_files/figure-gfm/Figure3-1.png)<!-- -->
