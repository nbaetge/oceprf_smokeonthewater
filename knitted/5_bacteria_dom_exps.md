Bacteria and DOM experiments
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
library(gt)
```

``` r
prod_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p5_bacteria_dom_exps.csv"
bcodmo_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/bco_dmo/BACTERIA_DOM.csv"
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
  ungroup() %>% 
  mutate(desc = case_when(trt == 1 ~ "Control",
                          trt == 2 ~ "Thomas Fire Ash",
                          trt == 3 ~ "Low Temp. Ash",
                          trt == 4 ~ "High Temp. Ash")) %>% 
  select(exp, trt, desc, tp, fgC_cell, fgN_cell) %>% 
  left_join(p5_fcm,.) %>%
  select(-trt)  %>% 
  rename(trt = desc)  %>% 
  group_by(exp, trt) %>%
  mutate(fgC_cell = ifelse(exp == "B2" & tp == 3, fgC_cell[tp == 2], fgC_cell)) %>% 
  mutate(interp_fgC_cell = zoo::na.approx(fgC_cell, dt, na.rm = F)) %>% 
  mutate(boc = (mean_cells_ml * interp_fgC_cell) / (12 * 10^6)) %>% 
  ungroup() %>% 
  mutate(interp_bocdoc = ifelse(is.na(fgC_cell) | exp == "B2" & tp == 3, T, F))
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
  mutate(doc = mean_toc - boc) %>% 
  mutate(toc_unc = sd_toc / sqrt(3),
         mean_toc_unc = round(mean(toc_unc), 1)) 


#inclusion points

incl <- doc %>% 
  select(exp, trt, tp, days) %>%
  mutate(bge_inclusion = ifelse(exp == "B16" & trt == "Control" & tp %in% c(0, 1, 2, 3, 4, 5) |
                                  exp == "B16" & trt != "Control" & tp %in% c(0, 1, 2, 3, 4, 5) |
                                  exp == "B15" & tp %in% c(0, 1, 2, 3) |
                                  exp == "B1"  & tp %in% c(0, 1, 2, 3, 4, 5) |
                                  exp == "B14" & trt != "Control" & tp %in% c(0, 1, 2, 3, 4) |
                                  exp == "B14" & trt == "Control" & tp %in% c(0, 1, 2, 3, 4, 5) |
                                  exp == "B3"  & tp %in% c(0, 1, 2) |
                                  exp == "B2" & tp %in% c(0, 1, 2, 3),
                                T, F)) %>% 
  filter(bge_inclusion == T)


#calculate boc rate of change as slope
boc_rate <- doc %>% 
  select(exp, trt, tp, days, boc) %>%
  drop_na(boc) %>%
  left_join(., incl) %>% 
  drop_na(bge_inclusion) %>% 
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
  left_join(., incl) %>% 
  drop_na(bge_inclusion) %>% 
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
         doc_bdl = ifelse(doc_rem <= 0.8, T, F),
         bge = ifelse(doc_bdl == F, boc_rate/doc_rate, NA)) %>% 
  ungroup() %>% 
  group_by(exp, tp) %>% 
  mutate(control_doc = doc[trt == "Control"]) %>% 
  ungroup() %>% 
  mutate(amend_doc = doc - control_doc) %>% 
  group_by(exp, trt) %>% 
  mutate(amend_doc_bioav = ifelse(trt != "Control", (amend_doc[tp == 0] - last(amend_doc, na_rm = T)) / amend_doc[tp == 0]* 100, NA )) %>% 
  ungroup() %>% 
  mutate(exp = paste("P", exp, sep = ""))
```

# Merge data with contextual in situ data

``` r
context_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p1_insitu_index.csv"

p5_prod <- read_csv(context_path) %>% 
  select(stn, exp, lat, lon,  cluster, composite_z) %>% 
  left_join(., p5_bge_bioav) %>% 
  arrange(cluster) %>% 
  drop_na(trt) %>% 
  left_join(., incl %>% mutate(exp = gsub("B", "PB", exp)))
```

``` r
bco_dmo <- p5_prod %>% 
  select(stn, lat, lon, dt, trt, mean_cells_ml, sd_cells_ml, boc, mean_toc, sd_toc) %>% 
  rename(bact_cells = mean_cells_ml,
         sd_bact_cells = sd_cells_ml) %>% 
  arrange(stn, trt, dt)
```

``` r
data4stats <- p5_prod %>% 
  select(exp, trt, cluster, mean_gr_d, bge, doc_rate, amend_doc, amend_doc_bioav) %>%
  distinct() %>% 
  # mutate(bge = ifelse(is.na(bge), mean(bge, na.rm = T), bge)) %>% 
  rename(bact_mu = 4) %>% 
  mutate(trt2 = ifelse(trt == "Control", "Control", "Ash leachate")) 
```

``` r
data4curves <- p5_prod %>%
  select(exp, trt, cluster, days, bge_inclusion, mean_toc, boc) %>%
  rename(toc = mean_toc) %>%
  pivot_longer(toc:boc, names_to = "var", values_to = "val") %>%
  mutate(
    var = case_when(
      var == "toc" ~ "bold(TOC~(µmol~C~L^-1))",
      var == "boc" ~ "bold(BOC~(µmol~C~L^-1))"
    )
  ) %>%
  left_join(
    .,
    p5_prod %>%
      select(exp, trt, days, bge_inclusion,  cluster, sd_toc) %>%
      mutate(var = "bold(TOC~(µmol~C~L^-1))") %>%
      rename(sd = sd_toc)
  ) %>%
   bind_rows(
    .,
    p5_prod %>%
      select(exp, trt, days, bge_inclusion,  cluster, doc, sd_toc) %>%
      mutate(var = "bold(DOC~(µmol~C~L^-1))") %>%
      rename(val = doc,
             sd = sd_toc)
  ) %>%
  bind_rows(
    .,
    p5_prod %>%
      select(exp, trt, days,  bge_inclusion, cluster, mean_cells_ml, sd_cells_ml) %>%
      mutate_at(vars(mean_cells_ml, sd_cells_ml), ~ . * 10 ^ 3) %>%
      mutate(var = "bold(Cells~(L^-1))") %>%
      rename(val = mean_cells_ml,
             sd = sd_cells_ml)
  ) %>%
  arrange(exp, trt, days, var) %>%
  mutate(bge_inclusion = ifelse(is.na(bge_inclusion), F, bge_inclusion)) %>% 
  mutate(bge_inclusion = ifelse(var %in% c("bold(TOC~(µmol~C~L^-1))", "bold(Cells~(L^-1))"), F, bge_inclusion)) %>% 
  mutate(cluster = case_when(cluster == "Low biomass" ~ "Low~biomass",
                             cluster == "High biomass" ~ "High~biomass")) 
```

# Save data

``` r
write_csv(p5_prod, prod_path)
write_csv(bco_dmo, bcodmo_path)
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
  axis.title = element_text(size = 44, face = "bold"),
  axis.text = element_text(size = 38),
  panel.spacing.x = unit(0.5, "cm"),
  panel.spacing.y = unit(1, "cm"),
  strip.text.x = element_text(
    size = 42,
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

pal4 = c("#785EF0", "#DC267F","#FFB000", "#FE6100")
pal3 = c("#DC267F","#FFB000", "#FE6100")
pal2 = c("#0947EA","#92EA33")
shape2 = c(21, 24)
plot_levels = c(
  "Control",
  "Thomas Fire Ash",
  "Low Temp. Ash",
  "High Temp. Ash",
  "Ash leachate",
  "DOC",
  "bold(Leachate)",
  "bold(Amendment)",
  "Low~biomass",
  "High~biomass",
  "Low biomass",
  "High biomass",
  "bold(Cells~(L^-1))",
  "bold(TOC~(µmol~C~L^-1))",
  "bold(BOC~(µmol~C~L^-1))",
  "bold(DOC~(µmol~C~L^-1))",
  "bold(Specific~growth~rate~(d^-1))",
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
curves <- ggplot(data4curves[!is.na(data4curves$val), ] %>% filter(var %in% c("bold(BOC~(µmol~C~L^-1))", "bold(DOC~(µmol~C~L^-1))")),
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
  geom_point( aes(shape = bge_inclusion),
    color = "black",
    size = 8,
    alpha = 0.8
  ) +
  scale_shape_manual(values = shape2, guide = guide_legend(override.aes = list(fill = "black"))) +
  scale_color_manual(values = pal4) +
  scale_fill_manual(values = pal4) +
  ggh4x::facet_nested(
    factor(var, levels = plot_levels) ~  factor(cluster, levels = plot_levels) + exp,
    scales = "free_y",
    labeller = label_parsed
  ) +
  labs(y = "", 
    x = expression(bold(Day)),
    fill = "",
    color = "",
    linetype = "",
    shape = "BGE calculation inclusion"
  ) +
  guides(color = "none", fill = guide_legend( 
    override.aes=list(shape = 21))) +
  theme_linedraw() +
  custom.theme 
```

``` r
# Fit the ANOVA model with interaction
mod1 <- aov(bact_mu ~ cluster * trt2, data = data4stats)

# Extract the ANOVA table
anova_results1 <- summary(mod1)

# Get p-values for significant terms
p_value_cluster <- anova_results1[[1]]["Pr(>F)"][1,1]
p_value_trt <- anova_results1[[1]]["Pr(>F)"][2,1]
# p_value_interaction <- anova_results2[[1]]["Pr(>F)"][3,1]

# Format p-values for annotations
formatted_p_cluster <- paste0("p = ", format.pval(p_value_cluster, digits = 2))
formatted_p_trt <- paste0("p = ", format.pval(p_value_trt, digits = 2))
# formatted_p_interaction <- paste0("p (cluster:trt) = ", format.pval(p_value_interaction, digits = 2))

cluster <- ggplot(data4stats, aes(x = factor(cluster, levels = plot_levels), y = bact_mu, fill = factor(cluster, levels = plot_levels))) +
  geom_boxplot(alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = factor(cluster, levels = plot_levels)),
    size = 7,
    position = position_jitterdodge(),
    alpha = 0.5, 
    shape = 21,
    color = "black"
  ) +
  labs(
    x = "SOM cluster",
    y = expression(bold(Bacterial~specific~growth~rate~(d^-1))),
    fill = "",
    color = ""
  ) +
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2) +
   theme_linedraw() +
  custom.theme +
  annotate(
    "text",
    x = 1,
    y = max(data4stats$bact_mu, na.rm = TRUE),
    label = formatted_p_cluster,
    size = 6,  vjust = 2
  ) +
   guides(color = "none") 


# Plot for Treatment Effect with ANOVA p-value using Box Plot
trt <-  ggplot(data4stats, aes(x = factor(trt2, levels = plot_levels), y = bact_mu, fill = factor(trt, levels = plot_levels))) +
  geom_boxplot(alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = factor(trt, levels = plot_levels)),
    size = 7,
    shape = 21, 
    position = position_jitterdodge(),
    alpha = 0.5,
    color = "black"
  ) +
  labs(
    x = "Treatment",
    y = expression(bold(Bacterial~specific~growth~rate~(d^-1))),
    fill = "",
    color = ""
  ) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
   theme_linedraw() +
  custom.theme +
  # Add p-value annotations for interaction
  annotate(
    "text",
    x = 1,
    y = max(data4stats$bact_mu, na.rm = TRUE),
    label = formatted_p_trt,
    size = 6,  vjust = 2
  ) +
   guides(color= "none") 

# Combine plots side by side
mu <- cluster + trt +  plot_layout(guides = 'collect') & theme(legend.position = "top")
```

``` r
# Fit the ANOVA model with interaction
mod2 <- aov(bge ~ cluster * trt2, data = data4stats)

# Extract the ANOVA table
anova_results2 <- summary(mod2)

# Get p-values for significant terms
p_value_cluster <- anova_results2[[1]]["Pr(>F)"][1,1]
p_value_trt <- anova_results2[[1]]["Pr(>F)"][2,1]
# p_value_interaction <- anova_results2[[1]]["Pr(>F)"][3,1]

# Format p-values for annotations
formatted_p_cluster <- paste0("p = ", format.pval(p_value_cluster, digits = 2))
formatted_p_trt <- paste0("p = ", format.pval(p_value_trt, digits = 2))
# formatted_p_interaction <- paste0("p (cluster:trt) = ", format.pval(p_value_interaction, digits = 2))

cluster2 <- ggplot(data4stats, aes(x = factor(cluster, levels = plot_levels), y = bge, fill = factor(cluster, levels = plot_levels))) +
  geom_boxplot(alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = factor(cluster, levels = plot_levels)),
    size = 7,
    position = position_jitterdodge(),
    alpha = 0.5, 
    shape = 21,
    color = "black"
  ) +
  labs(
    x = "SOM cluster",
    y = expression(bold(BGE)),
    fill = "",
    color = ""
  ) +
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2) +
   theme_linedraw() +
  custom.theme +
  annotate(
    "text",
    x = 1,
    y = max(data4stats$bge, na.rm = TRUE),
    label = formatted_p_cluster,
    size = 6,  vjust = 2
  ) +
   guides(color = "none") 


# Plot for Treatment Effect with ANOVA p-value using Box Plot
trt2 <-  ggplot(data4stats, aes(x = factor(trt2, levels = plot_levels), y = bge, fill = factor(trt, levels = plot_levels))) +
  geom_boxplot(alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = factor(trt, levels = plot_levels)),
    size = 7,
    shape = 21, 
    position = position_jitterdodge(),
    alpha = 0.5,
    color = "black"
  ) +
  labs(
    x = "Treatment",
    y = expression(bold(BGE)),
    fill = "",
    color = ""
  ) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
   theme_linedraw() +
  custom.theme +
  # Add p-value annotations for interaction
  annotate(
    "text",
    x = 1,
    y = max(data4stats$bge, na.rm = TRUE),
    label = formatted_p_trt,
    size = 6,  vjust = 2
  ) +
   guides(color= "none") 

# Combine plots side by side
bge <- cluster2 + trt2 +  plot_layout(guides = 'collect') & theme(legend.position = "top")
```

``` r
fig3 <- (
  (curves)
) /(mu)/(bge) + plot_layout(heights = c(0.8, 0.7, 0.7), guides = "collect") + plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 45), legend.position = "top", legend.box="vertical")

ggsave(fig3, file = "~/github/oceprf_smokeonthewater/prod/figs/bact.svg", width = 37, height = 48)
```

``` r
# Fit the ANOVA model with interaction
mod3 <- aov(doc_rate ~ cluster * trt2, data = data4stats)

# Extract the ANOVA table
anova_results3 <- summary(mod3)

# Get p-values for significant terms
p_value_cluster <- anova_results3[[1]]["Pr(>F)"][1,1]
p_value_trt <- anova_results3[[1]]["Pr(>F)"][2,1]
p_value_interaction <- anova_results3[[1]]["Pr(>F)"][3,1]

# Format p-values for annotations
formatted_p_cluster <- paste0("p = ", format.pval(p_value_cluster, digits = 2))
formatted_p_trt <- paste0("p = ", format.pval(p_value_trt, digits = 2))
formatted_p_interaction <- paste0("p (Cluster : Treatment) = ", format.pval(p_value_interaction, digits = 2))

cluster <- ggplot(data4stats, aes(x = factor(cluster, levels = plot_levels), y = doc_rate, fill = factor(cluster, levels = plot_levels))) +
  geom_boxplot(alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = factor(cluster, levels = plot_levels)),
    size = 7,
    position = position_jitterdodge(),
    alpha = 0.5, 
    shape = 21,
    color = "black"
  ) +
  labs(
    x = "SOM cluster",
    y = expression(bold(DOC~remineralization~rate~(µmol~C~L^-1~d^-1))),
    fill = "",
    color = ""
  ) +
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2) +
   theme_linedraw() +
  custom.theme +
   guides(color = "none") 


# Plot for Treatment Effect with ANOVA p-value using Box Plot
trt <-  ggplot(data4stats, aes(x = factor(trt2, levels = plot_levels), y = doc_rate, fill = factor(trt, levels = plot_levels))) +
  geom_boxplot(alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = factor(trt, levels = plot_levels)),
    size = 7,
    shape = 21, 
    position = position_jitterdodge(),
    alpha = 0.5,
    color = "black"
  ) +
  labs(
    x = "Treatment",
    y = expression(bold(DOC~remineralization~rate~(µmol~C~L^-1~d^-1))),
    fill = "",
    color = ""
  ) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
   theme_linedraw() +
  custom.theme +
  # Add p-value annotations for interaction
  annotate(
    "text",
    x = 1.5,  # Position for p-value of interaction
    y = max(data4stats$doc_rate, na.rm = TRUE) * 0.95,
    label = formatted_p_interaction,
    size = 6, vjust = 2
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(data4stats$doc_rate, na.rm = TRUE) * 0.85,
    label = formatted_p_cluster,
    size = 6, vjust = 2
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(data4stats$doc_rate, na.rm = TRUE) * 0.75,
    label = formatted_p_trt,
    size = 6, vjust = 2
  ) +
   guides(color= "none") 

# Combine plots side by side
doc_rate <- cluster + trt +  plot_layout(guides = 'collect') & theme(legend.position = "top")
```

``` r
# Fit the ANOVA model with interaction
mod4 <- aov(amend_doc_bioav ~ cluster * trt, data = data4stats %>% filter(!trt == "Control"))

# Extract the ANOVA table
anova_results4 <- summary(mod4)

# Get p-values for significant terms
p_value_cluster <- anova_results4[[1]]["Pr(>F)"][1,1]
p_value_trt <- anova_results4[[1]]["Pr(>F)"][2,1]
# p_value_interaction <- anova_results4[[1]]["Pr(>F)"][3,1]

# Format p-values for annotations
formatted_p_cluster <- paste0("p = ", format.pval(p_value_cluster, digits = 2))
formatted_p_trt <- paste0("p = ", format.pval(p_value_trt, digits = 2))
# formatted_p_interaction <- paste0("p (Cluster : Treatment) = ", format.pval(p_value_interaction, digits = 2))

cluster <- ggplot(data4stats, aes(x = factor(cluster, levels = plot_levels), y = amend_doc_bioav, fill = factor(cluster, levels = plot_levels))) +
  geom_boxplot(alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = factor(cluster, levels = plot_levels)),
    size = 7,
    position = position_jitterdodge(),
    alpha = 0.5, 
    shape = 21,
    color = "black"
  ) +
  labs(
    x = "SOM cluster",
    y = expression(bold(Bioavailability~of~amended~DOC~('%'))),
    fill = "",
    color = ""
  ) +
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2) +
   theme_linedraw() +
  custom.theme +
   guides(color = "none") 


# Plot for Treatment Effect with ANOVA p-value using Box Plot
trt <-  ggplot(data4stats %>% filter(trt != "Control"), aes(x = factor(trt, levels = plot_levels), y = amend_doc_bioav, fill = factor(trt, levels = plot_levels))) +
  geom_boxplot(alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = factor(trt, levels = plot_levels)),
    size = 7,
    shape = 21, 
    position = position_jitterdodge(),
    alpha = 0.5,
    color = "black"
  ) +
  labs(
    x = "Treatment",
    y = expression(bold(Bioavailability~of~amended~DOC~('%'))),
    fill = "",
    color = ""
  ) +
  scale_fill_manual(values = pal3) +
  scale_color_manual(values = pal3) +
   theme_linedraw() +
  custom.theme +
  # Add p-value annotations for interaction
  annotate(
    "text",
    x = 1.5,  # Position for p-value of interaction
    y = max(data4stats$amend_doc_bioav, na.rm = TRUE) * 0.95,
    label = formatted_p_interaction,
    size = 6, vjust = 2
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(data4stats$amend_doc_bioav, na.rm = TRUE) * 0.85,
    label = formatted_p_cluster,
    size = 6, vjust = 2
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(data4stats$amend_doc_bioav, na.rm = TRUE) * 0.75,
    label = formatted_p_trt,
    size = 6, vjust = 2
  ) +
   guides(color= "none") 

# Combine plots side by side
bioav <- cluster + trt +  plot_layout(guides = 'collect') & theme(legend.position = "top")
```

``` r
fig4 <- (
 doc_rate)/(bioav) + plot_layout(heights = c(0.7, 0.7), guides = "collect") + plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 45), legend.position = "top", legend.box="vertical")

ggsave(fig4, file = "~/github/oceprf_smokeonthewater/prod/figs/doc.svg", width = 37, height = 40)
```

# Supplementary Table

``` r
table_data <- p5_prod %>% 
  select(exp, trt, days, interp_bocdoc, bge_inclusion, interp_fgC_cell, boc_rate, doc_rate) %>% 
  mutate_at(vars(contains(c("rate"))), round, 2) %>% 
  mutate_at(vars(contains("interp_fgC_cell")), round, 1) %>% 
  mutate(bge_inclusion = ifelse(is.na(bge_inclusion), F, bge_inclusion))
  
gt_tbl <- gt(table_data) 

table <- 
  gt_tbl |>
  tab_header(
    title = md("**Table S2.** Cell-specific bacterial carbon, bacterial organic carbon rate changes, and DOC removal rates"),
  ) |>
  cols_label(
    exp = html("Experimental site"),
    trt = html("Treatment"),
    days = html("Days"), 
    interp_bocdoc = html("Interpolated cell-specific bacterial carbon"),
    bge_inclusion = html("BGE calculation inclusion"),
    interp_fgC_cell = html("Cell-specific bacterial carbon"), 
    boc_rate = html("Rate change in bacterial organic carbon"),
    doc_rate = html("Rate change in DOC"),
  ) |>
  tab_spanner(
    label = html("µmol C L<sup>-1</sup> d<sup>-1"),
    columns = c(boc_rate, doc_rate)
  ) |>
  tab_spanner(
    label = html("fg C cell<sup>-1"),
    columns = c(interp_fgC_cell)
  ) 

table_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/knitted/5_bacteria_dom_exps_files/bact_table.html"

gtsave(table, table_path)
```

# Supplementary Figure

``` r
suppfig1 <- ggplot(data4curves[!is.na(data4curves$val), ] %>% filter(!var %in% c("bold(BOC~(µmol~C~L^-1))", "bold(DOC~(µmol~C~L^-1))")),
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
  geom_point( aes(shape = bge_inclusion),
    color = "black",
    size = 8,
    alpha = 0.8
  ) +
  scale_shape_manual(values = shape2, guide = guide_legend(override.aes = list(fill = "black"))) +
  scale_color_manual(values = pal4) +
  scale_fill_manual(values = pal4) +
  ggh4x::facet_nested(
    factor(var, levels = plot_levels) ~  factor(cluster, levels = plot_levels) + exp,
    scales = "free_y",
    labeller = label_parsed
  ) +
  labs(y = "", 
    x = expression(bold(Day)),
    fill = "",
    color = "",
    linetype = "",
    shape = "BGE calculation inclusion"
  ) +
  guides(color = "none", fill = guide_legend( 
    override.aes=list(shape = 21))) +
  theme_linedraw() + 
  custom.theme 

ggsave(suppfig1, file = "~/github/oceprf_smokeonthewater/prod/figs/bact_suppl.svg", width = 37, height = 20)
```
