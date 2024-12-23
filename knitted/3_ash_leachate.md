Ash leachate
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
library(staRdom)
library(gt)
```

# IMPORT DATA & DO CALCS

``` r
data_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r3_leachate_nuts_doc_tm.xlsx"

cdom_abs_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r3_cdom_abs"

cdom_eems_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r3_cdom_eems"

prod1_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p3_ash_leachate_nuts_tm_doc.csv"

prod2_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p3_ash_leachate_cdom.csv"

prod3_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p3_ash_leachate_eems.csv"
```

## Leachate DOC from high temperature combustion (units = µmol/L )

``` r
leach_doc_data <- readxl::read_xlsx(data_path, sheet = 3) %>%
  mutate(
    id = c(1, 1, 2, 2, 3, 3, 4, 4),
    trt =  case_when(
      id == 1 ~ "Control",
      id == 4 ~ "Thomas Fire Ash",
      id == 2 ~ "Low Temp. Ash",
      id == 3 ~ "High Temp. Ash"
    ),
    .after = 1
  ) %>%
  select(3, 5) %>%
  group_by(trt) %>%
  mutate(ave_doc = mean(doc),
         sd_doc = sd(doc)) %>%
  select(-doc) %>%
  distinct() %>%
  ungroup() %>%
  mutate(blank_doc = first(ave_doc)) %>%
  mutate(corr_doc = ave_doc - blank_doc) %>% 
  select(-blank_doc)
```

## Actual DOC amendment in incubations (from 6 incubations)

``` r
amend_doc_data <- readxl::read_xlsx(data_path, sheet = 4, skip = 27) %>% 
  drop_na() %>% 
  select(1,3) %>% 
  rename(id = 1, 
         toc = 2) %>% 
  separate(id, into = c("exp", "btl", "tp"), sep = "-") %>% 
  filter(tp == 0) %>% 
  select(-tp) %>% 
  group_by(exp, btl) %>% 
  mutate(exp = paste("P", exp, sep = ""),
         mean_toc = mean(toc),
         sd_toc = sd(toc)) %>% 
  ungroup() %>% 
  select(-toc) %>% 
  distinct() %>% 
  mutate(trt =  case_when(btl == 1 ~ "Control",
                          btl == 2 ~ "Thomas Fire Ash",
                          btl == 3 ~ "Low Temp. Ash",
                          btl == 4 ~ "High Temp. Ash"),
         .after = btl) %>% 
  select(-btl) %>% 
  group_by(exp) %>% 
  mutate(control_toc = mean_toc[trt == "Control"]) %>% 
  ungroup() %>% 
  mutate(added_doc = mean_toc - control_toc) %>% 
  select(trt, added_doc) %>% 
  filter(trt != "Control") %>% 
  group_by(trt) %>% 
  summarize_at(vars(added_doc), list(mean_added_doc = mean, sd_added_doc = sd)) 
```

## Nutrients from flow injection analysis (units = µmol/L )

concentrations below the limit of detection were replaced with values
representing 1/2 the limit of detection

``` r
leach_nut_data <- readxl::read_xlsx(data_path, sheet = 1) %>%
    mutate(df = case_when(
    sample == "TFA" ~ 50,
    sample == "LTA" ~ 4,
    sample == "HTA" ~ 5,
    sample == "Blank" ~ 1
  )) %>% #account for diluting the leachate prior to analysis (20 ml vol of diluent/vol of leachate) 
  mutate(
    po4 = ifelse(po4_bdl == "T", 0.05, po4),
    sio4 = ifelse(sio4_bdl == "T", 0.1, sio4),
    no2 = ifelse(no2_bdl == "T", 0.05, no2),
    no2no3 = ifelse(no2no3_bdl == "T", 0.1, no2no3),
    nh4 = ifelse(nh4_bdl == "T", 0.1, nh4)
  ) %>%
  select(-contains("bdl")) %>%
  mutate_at(vars(po4:nh4), ~.*df) %>% 
  mutate(
    blank_po4 = first(po4),
    blank_sio4 = first(sio4),
    blank_no2 = first(no2),
    blank_no2no3 = first(no2no3),
    blank_nh4 = first(nh4)
  ) %>%
  mutate(
    corr_po4 = po4 - blank_po4,
    corr_sio4 = sio4 - blank_sio4,
    corr_no2 = no2 - blank_no2,
    corr_no2no3 = no2no3 - blank_no2no3,
    corr_nh4 = nh4 - blank_nh4
  ) %>%
  select(-contains("blank")) %>%
  mutate(across(po4:corr_nh4, ~ ifelse(.x < 0, 0, .x))) %>%
  mutate(across(po4:corr_nh4, ~ round(., 2))) %>%
  mutate(
    trt = case_when(
      sample == "Blank" ~ "Control",
      sample == "TFA" ~ "Thomas Fire Ash",
      sample == "LTA" ~ "Low Temp. Ash",
      sample == "HTA" ~ "High Temp. Ash"
    ),
    .after = sample
  )  %>%
  select(-sample, -df) 
```

## Actual nutrient amendment in incubations (from 6 incubations)

``` r
amend_nut_data <- readxl::read_xlsx(data_path, sheet = 2)  %>%
  filter(tp == 0) %>%
  mutate(
    po4 = ifelse(po4_bdl == "T", 0.05, po4),
    sio4 = ifelse(sio4_bdl == "T", 0.1, sio4),
    no2 = ifelse(no2_bdl == "T", 0.05, no2),
    no2no3 = ifelse(no2no3_bdl == "T", 0.1, no2no3),
    nh4 = ifelse(nh4_bdl == "T", 0.1, nh4)
  ) %>% 
  select(-contains("bdl")) %>%
  mutate(
    trt =  case_when(
      btl == 100 ~ "Control",
      btl == "100-TFA" ~ "Thomas Fire Ash",
      btl == "100-LTA" ~ "Low Temp. Ash",
      btl == "100-HTA" ~ "High Temp. Ash"
    ),
    .after = btl
  ) %>%
  select(-btl) %>%
  group_by(exp) %>%
  mutate(
    control_po4 = po4[trt == "Control"],
    control_sio4 = sio4[trt == "Control"],
    control_no2 = no2[trt == "Control"],
    control_no2no3 = no2no3[trt == "Control"],
    control_nh4 = nh4[trt == "Control"]
  ) %>%
  ungroup() %>%
  mutate(
    added_po4 = po4 - control_po4,
    added_sio4 = sio4 - control_sio4,
    added_no2 = no2 - control_no2,
    added_no2no3 = no2no3 - control_no2no3,
    added_nh4 = nh4 - control_nh4
  ) %>%
  select(trt, contains("added")) %>%
  filter(trt != "Control") %>%
  mutate(across(added_po4:added_nh4, ~ ifelse(.x < 0, 0, .x))) %>%
  group_by(trt) %>%
  summarize_at(vars(added_po4:added_nh4), list(mean = mean, sd = sd)) %>%
  ungroup() %>%
  mutate(across(added_po4_mean:added_nh4_sd, ~ round(., 2)))
```

## Trace metals from LC-IPMS (units = nmol/L)

``` r
leach_tm_data <- readxl::read_xlsx(data_path, sheet = 5) %>% 
  mutate(tm2 = str_sub(tm, -2),
         iso = as.numeric(str_sub(tm, 1, -3))) %>% 
  select(-tm) %>% 
  rename(tm = tm2) %>% 
  select(sample, tm:iso, ave_tm, sd_tm) %>% 
  mutate(sample =  case_when(sample =="Blank" ~ "Blank",
                          sample  == "TFA" ~ "Thomas Fire Ash",
                          sample  == "LTA" ~ "Low Temp. Ash",
                          sample  == "HTA" ~ "High Temp. Ash")) %>% 
  # mutate(mass = case_when(tm == "Mn" ~ 54.938,
  #                         tm == "Fe" ~ 55.845,
  #                         tm == "Co" ~ 58.933,
  #                         tm == "Ni" ~ 58.693,
  #                         tm == "Cu" ~ 63.546,
  #                         tm == "Zn" ~ 65.380, 
  #                         tm == "Cd" ~ 112.41,
  #                         tm == "Pb" ~ 207.20))
  mutate(nM = ave_tm / 10^6 / iso * 10^9,
         sd_nM = sd_tm / 10^6 / iso * 10^9 ) %>% 
  select(sample:tm, nM, sd_nM)
```

## Estimated trace metal addition based on DOC addition (units = nmol/L)

``` r
amend_tm_data <- leach_doc_data %>% 
  select(trt, corr_doc) %>% 
  filter(trt != "Control") %>% 
  left_join(., amend_doc_data %>% select(trt, mean_added_doc)) %>% 
  mutate(df = mean_added_doc/corr_doc) %>% 
  select(trt, df) %>% 
  left_join(leach_tm_data %>% rename(trt = sample), .) %>% 
  filter(trt != "Blank") %>% 
  mutate(added_tm = nM * df) %>% 
  select(trt, tm, added_tm)
```

## Combine leachate data

``` r
leach <- leach_tm_data %>%
  rename(
    trt = sample,
    analyte = tm,
    val = nM,
    sd = sd_nM
  ) %>%
  mutate(trt = ifelse(trt == "Blank", "Control", trt),
         facet = "bold(Trace~metals~(nmol~L^-1))") %>%
  bind_rows(
    .,
    leach_nut_data %>%
      select(trt:nh4) %>%
      pivot_longer(c(2:6), names_to = "analyte", values_to = "val") %>%
      mutate(facet = "bold(Inorganic~nutrients~(µmol~L^-1))") 
  ) %>%
  bind_rows(
    .,
    leach_doc_data %>%
      select(trt:sd_doc) %>%
      mutate(analyte = "DOC") %>%
      rename(val = ave_doc,
             sd = sd_doc) %>% 
      mutate(facet = "bold(Organic~carbon~(µmol~L^-1))")
  ) %>%
  mutate(
    analyte = case_when(
      analyte == "nh4" ~ "NH[4]",
      analyte == "no2" ~ "NO[2]",
      analyte == "no2no3" ~ "NO[2]+NO[3]",
      analyte == "po4" ~ "PO[4]",
      analyte == "sio4" ~ "SiO[4]",
      .default = as.character(analyte)
    )
  )
```

## Combine amendment data

``` r
amend <-  amend_nut_data %>%
  select(trt:added_nh4_mean) %>%
  pivot_longer(c(2:6), names_to = "analyte", values_to = "val") %>%
  mutate(analyte = gsub("added_" , "", analyte),
         analyte = gsub("_mean" , "", analyte)) %>%
  left_join(
    .,
    amend_nut_data %>%
      select(trt, added_po4_sd:added_nh4_sd) %>%
      pivot_longer(c(2:6), names_to = "analyte", values_to = "sd") %>%
      mutate(
        analyte = gsub("added_" , "", analyte),
        analyte = gsub("_sd" , "", analyte)
      )
  ) %>%
  mutate(facet = "bold(Inorganic~nutrients~(µmol~L^-1))") %>% 
  bind_rows(.,  amend_tm_data %>%
              rename(analyte = tm,
                     val = added_tm) %>% 
              mutate(facet = "bold(Trace~metals~(nmol~L^-1))")) %>%
  bind_rows(
    .,
    amend_doc_data %>%
      select(trt:sd_added_doc) %>%
      mutate(analyte = "DOC") %>%
      rename(val = mean_added_doc,
             sd = sd_added_doc) %>% 
      mutate(facet = "bold(Organic~carbon~(µmol~L^-1))")
  ) %>%
  mutate(
    analyte = case_when(
      analyte == "nh4" ~ "NH[4]",
      analyte == "no2" ~ "NO[2]",
      analyte == "no2no3" ~ "NO[2]+NO[3]",
      analyte == "po4" ~ "PO[4]",
      analyte == "sio4" ~ "SiO[4]",
      .default = as.character(analyte)
    )
  )
```

## Combine leachate and amendment data

``` r
conc <- leach %>% 
  mutate(type = "bold(Leachate)") %>% 
  bind_rows(., amend %>% 
              mutate(type = "bold(Amendment)"))
```

## CDOM Absorption

``` r
cores <- detectCores(logical = FALSE)
absorbance <- absorbance_read(cdom_abs_path, cores = cores)
```

``` r
absorbance <- abs_blcor(absorbance, wlrange = c(680, 700))
```

``` r
abs <- absorbance %>%
  mutate(HTA = HTA - Blank,
         LTA = LTA - Blank,
         TFA = TFA - Blank) %>%
  select(-Blank) %>%
  pivot_longer(c(2:4), names_to = "trt", values_to = "abs") %>%
  mutate(
    trt =  case_when(
      trt == "TFA" ~ "Thomas Fire Ash",
      trt == "LTA" ~ "Low Temp. Ash",
      trt == "HTA" ~ "High Temp. Ash"
    )
  ) %>%
  rename(wl = 1) %>%
  filter(between(wl, 250, 500)) 
```

``` r
abs_slopes <- absorbance %>%
  mutate(HTA = HTA - Blank,
         LTA = LTA - Blank,
         TFA = TFA - Blank)

slope_params <- abs_parms(abs_slopes, cuvl = 1, cores = cores) %>%
  select(sample, E2_E3, SR) %>%
  pivot_longer(c(2:3), names_to = "slope_parameter", values_to = "slope") %>%
  rename(trt = 1) %>%
  mutate(
    trt =  case_when(
      trt == "Blank" ~ "Blank",
      trt == "TFA" ~ "Thomas Fire Ash",
      trt == "LTA" ~ "Low Temp. Ash",
      trt == "HTA" ~ "High Temp. Ash"
    )
  ) %>%
  mutate(
    slope_parameter = case_when(
      slope_parameter == "E2_E3" ~ "bold(E2:E3)",
      slope_parameter == "SR" ~ "bold(S[R])"
    )
  ) %>%
  filter(trt != "Blank")
```

## CDOM Fluorescence

``` r
eem_list <- eem_read(cdom_eems_path, import_function = "fluoromax4")
```

``` r
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(15, 22, 15, 22)

eem_list <-
  eem_rem_scat(eem_list,
               remove_scatter = remove_scatter,
               remove_scatter_width = remove_scatter_width)
```

``` r
nelson <-
  tibble(
    peak = c("A", "T", "M", "C", "B", "N"),
    ex = c(260, 275, 300, 340, 275, 280),
    em = c(430, 340, 390, 440, 305, 370)
  )
```

``` r
blank <- as.data.frame(eem_list[[1]]) %>% mutate(sample = "Blank")
hta <-
  as.data.frame(eem_list[[2]]) %>% mutate(sample = "High Temp. Ash")
lta <-
  as.data.frame(eem_list[[3]]) %>% mutate(sample = "Low Temp. Ash")
tfa <-
  as.data.frame(eem_list[[4]]) %>% mutate(sample = "Thomas Fire Ash")

combined <- bind_rows(list(blank, hta, lta, tfa)) %>%
  mutate_at(vars(em, ex), as.numeric)
```

``` r
eems <-
  combined %>% mutate(blank = ifelse(sample == "Blank", value[sample == "Blank"], NA)) %>%
  group_by(em, ex) %>%
  fill(blank, .direction  = "updown") %>%
  ungroup() %>%
  mutate(blank_corr = value - blank) %>%
  filter(sample != "Blank") %>% 
  select(sample, ex, em, blank_corr) %>% 
  rename(trt = sample,
         eems = blank_corr)
```

## Combine CDOM data

``` r
cdom <- left_join(eems, abs %>% rename(ex = wl))
```

# PLOTS

``` r
custom.theme <- theme(
  legend.position = "top",
  legend.title = element_text(size = 35),
  legend.key.size = unit(3, "cm"),
  legend.key.spacing.x = unit(0.75, "cm"),
  legend.text = element_text(size = 38),
  axis.title = element_text(size = 35, face = "bold"),
  axis.text = element_text(size = 28),
  # panel.spacing.x = unit(2, "cm"),
  panel.spacing.x = unit(0.5, "cm"),
  panel.spacing.y = unit(0.25, "cm"),
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

pal4 = c("#785EF0", "#DC267F","#FFB000", "#FE6100")
pal3 = c("#DC267F","#FFB000", "#FE6100")
plot_levels = c(
  "Control",
  "Thomas Fire Ash",
  "Low Temp. Ash",
  "High Temp. Ash",
  "Mn",
  "Fe",
  "Co",
  "Ni",
  "Cu",
  "Zn",
  "Cd",
  "Pb",
  "NH[4]",
  "NO[2]",
  "NO[2]+NO[3]",
  "PO[4]",
  "SiO[4]",
  "DOC",
  "bold(Leachate)",
  "bold(Amendment)",
  "bold(E2:E3)",
  "bold(S[R])"
)
```

## Nutrients, trace metals, DOC

``` r
nuts <- ggplot(data = conc, aes(
  x = factor(analyte, levels = plot_levels),
  y = val,
  color = factor(trt, levels = plot_levels),
)) +
  geom_errorbar(aes(ymin = val - sd,
                    ymax = val + sd),
                width = .2,
                position = position_dodge(.4)) +
  # geom_linerange(
  #   aes(ymin = 0, ymax = val),
  #   position = position_dodge(width = 0.4),
  #   linewidth = 1,
  #   alpha = 0.8
  # ) +
  geom_point(position = position_dodge(width = 0.4),
             size = 10,
             alpha = 0.8) +
  scale_color_manual(values = pal4) +
  scale_fill_manual(values = pal4) +
  scale_x_discrete(labels = scales::parse_format()) +
  labs(
    y = expression(bold(Concentration)),
    x = "",
    fill = "",
    color = "",
    linetype = ""
  ) +
  # guides(color = "none") +
  ggh4x::facet_nested_wrap(factor(type, levels = plot_levels) ~ facet,
                           scales = "free",
                           labeller = label_parsed) +
  theme_linedraw() +
  custom.theme 
```

## CDOM

``` r
abs_plot <-
  ggplot(cdom %>%  filter(between(ex, 250, 500)) %>% select(trt, ex, abs) %>% distinct(),
         aes(x = ex, y = abs)) +
  geom_line(aes(color = factor(trt, levels = plot_levels)), linewidth = 4, alpha = 0.8) +
  scale_color_manual(values = pal3) +
  scale_fill_manual(values = pal3) +
  labs(
    y = expression(bold(a[CDOM] ~ (m ^ -1))),
    x = expression(bold(Wavelength ~ (nm))),
    color = "",
    linetype = ""
  ) +
  theme_linedraw() +
  custom.theme +
  theme(legend.box = "vertical", legend.margin = margin())
```

``` r
slope_plot <-
  ggplot(slope_params, aes(
    x = slope_parameter,
    y = slope,
    color = factor(trt, levels = plot_levels)
  )) +
  # geom_linerange(
  #   aes(ymin = 0, ymax = slope),
  #   position = position_dodge(width = 0.4),
  #   linewidth = 1.5,
  #   alpha = 0.8
  # ) +
  geom_point(position = position_dodge(width = 0.4),
             size = 10,
             alpha = 0.8) +
  ggh4x::facet_nested_wrap(
    ~ factor(slope_parameter, levels = plot_levels),
    scales = "free",
    labeller = label_parsed
  ) +
  scale_color_manual(values = pal3) +
  scale_x_discrete(labels = scales::parse_format()) +
  labs(y = expression(bold(a[CDOM] ~ Ratio)),
       x = "",
       color = "") +
  theme_linedraw() +
  custom.theme +
  theme(axis.text.x = element_blank())
```

``` r
eems_plot <- ggplot() +
  geom_raster(data = eems,
              aes(
                x = ex,
                y = em,
                z = eems,
                fill = eems
              ),
              interpolate = T) +
  geom_contour(colour = "black",
               size = 0.3,
               binwidth = 2) +
  labs(
    x = expression(bold(Excitation ~ (nm))),
    y = expression(bold(Emission ~ (nm))),
    fill = expression(bold(CDOM ~ Fluorescence ~ intensity ~ (ppb ~ QSE)))
  ) +
  viridis::scale_fill_viridis(option = "H", na.value = NA) +
  geom_label(
    data = nelson,
    aes(x = ex, y = em, label = peak),
    alpha = 0.7,
    size = 12
  ) +
  facet_wrap( ~ factor(trt, levels = plot_levels), scales = "free") +
  theme_test() +
  custom.theme +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width  = unit(3, "cm"),
        legend.title = element_text(vjust = 1, hjust = 1.5))
```

``` r
design <- c(
  area(1, 1, 2,3),
  area(3, 1, 3, 2),
  area(3,3),
  area(4, 1, 4, 3)
)

fig2 <- nuts / (abs_plot + guides(color = "none")) / (slope_plot + guides(color = "none")) / eems_plot + plot_layout(design = design) + plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 45))

ggsave(fig2, file = "~/github/oceprf_smokeonthewater/prod/figs/leachate.svg", width = 32, height = 35.2)
```

    ## Warning: Removed 3132 rows containing missing values or values outside the scale range
    ## (`geom_raster()`).

# SAVE

``` r
write_csv(conc, prod1_path)
write_csv(cdom, prod2_path)
write_csv(slope_params, prod3_path)
```

# Supporting pH table

``` r
ph_data <- readxl::read_xlsx("~/github/oceprf_smokeonthewater/raw/r9_pH.xlsx", sheet = 4) %>% 
  mutate_at(c(3:4), ~round(.,4))
 
gt_ph<- gt(ph_data) 

phtable <- 
  gt_ph |>
  tab_header(
    title = md("**Table S3.** Spectrophotometric pH measurements of phytoplankton dilution experiment bottles"),
  )  |>
  cols_label(
    pH = html("pH of 5 replicates"),
    sd = html("Std. deviation")
  )   |>
  tab_source_note(source_note = md(
    "Precision of 5 replicates for a single pH sample is ±0.0004 (Liu, X., M. C. Patsavas, and R. H. Byrne. 2011. Purification and Characterization of meta-Cresol Purple for Spectrophotometric Seawater pH Measurements. Environ. Sci. Technol. 45: 4862–4868. doi:10.1021/es200665d)"
  )) |>
  tab_source_note(source_note = md(
    "Accuracy of spectrophotometric pH sample is ±0.003 (Orr, J. C., J.-M. Epitalon, A. G. Dickson, and J.-P. Gattuso. 2018. Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207: 84–107. doi:10.1016/j.marchem.2018.10.006)"
  )) 

gtsave(phtable, "/Users/nicholasbaetge/github/oceprf_smokeonthewater/knitted/3_ash_leachate_files/pH_table.html")
```
