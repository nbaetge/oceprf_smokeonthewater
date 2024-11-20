DOC concentrations during ACIDD
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

## Import and tidy data

``` r
path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r7_ACIDD_biogeochem_data_.xlsx"
```

``` r
data <- readxl::read_xlsx(path) %>%
  mutate(across(.cols = everything(),
                .fns = ~ ifelse(.x == -999, NA, .x))) %>% 
  mutate(Time_Stamp = gsub("T", " ", Time_Stamp),
         local_dt = ymd_hm(Time_Stamp, tz = "America/Los_Angeles"),
         dt = with_tz(local_dt, tzone = "UTC"),
         .after = Station
         ) %>% 
  select(-Time_Stamp, -CruiseCN, -Leg) %>% 
  rename(Cast = SCN)

doc_data <- data %>% 
  select(Station, Cast, local_dt, dt, Latitude, Longitude, Z, DOC) %>% 
  drop_na(DOC) %>% 
  filter(Z < 10, !Station == 0) %>% 
  mutate(local_date = as_date(local_dt), .after = local_dt) %>% 
  mutate(group = ifelse(local_date %in% c(as_date("2017/12/17"), as_date("2017-12-18")), "12/17 - 12/18", "12/19 - 12/22")  )
```

# Stats

``` r
mod1 <- aov(DOC ~ group, data = doc_data)

# Extract the ANOVA table
anova_results1 <- summary(mod1)

# Get p-values for significant terms
p_value_cluster <- anova_results1[[1]]["Pr(>F)"][1,1]
p_value_trt <- anova_results1[[1]]["Pr(>F)"][2,1]
# p_value_interaction <- anova_results2[[1]]["Pr(>F)"][3,1]

# Format p-values for annotations
formatted_p_cluster <- paste0("p = ", format.pval(p_value_cluster, digits = 2))
formatted_p_trt <- paste0("p = ", format.pval(p_value_trt, digits = 2))
```

# Plot

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
    size = 35,
    color = "white",
    face = "bold"
  ),
  strip.text.y = element_text(
    size = 42,
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

pal2 = c( "#DC267F", "#785EF0")
```

``` r
doc_plot <- ggplot(doc_data, aes(x = group, y = DOC, group = group)) +
  geom_boxplot(
    aes(fill = group),
    alpha = 0.7,
    outlier.shape = NA,
    width = 0.4,
    position = position_dodge(0.5)) +
  geom_jitter(
    aes(fill = group),
    size = 3,
    alpha = 0.5, 
    shape = 21,
    color = "black"
  ) +
  labs(
    x = expression(bold(Date~range~(2017))),
    y = expression(bold(DOC~(µmol~C~L^-1))),
    color = "",
    fill = "",
  ) +
  scale_color_manual(values = pal2) +
  scale_fill_manual(values = pal2) +
  guides(color = "none",
         fill = "none") +
  theme_linedraw() +
  custom.theme +
   annotate(
    "text",
    x = 1,
    y = max(doc_data$DOC, na.rm = TRUE),
    label = formatted_p_cluster,
    size = 6,  vjust = 2
  ) +
   guides(color = "none")

ggsave(doc_plot, file = "~/github/oceprf_smokeonthewater/prod/figs/acidd_doc.svg", width = 12, height = 10)
```

``` r
doc_data %>% 
  filter(local_date %in% c(as_date("2017-12-17"), as_date("2017-12-18"))) %>% 
  summarize_at(vars(DOC), list(mean = mean, sd = sd))
```

    ## # A tibble: 1 × 2
    ##    mean    sd
    ##   <dbl> <dbl>
    ## 1  68.0  1.14

``` r
doc_data %>% 
  filter(!local_date %in% c(as_date("2017-12-17"), as_date("2017-12-18"))) %>% 
  summarize_at(vars(DOC), list(mean = mean, sd = sd))
```

    ## # A tibble: 1 × 2
    ##    mean    sd
    ##   <dbl> <dbl>
    ## 1  66.5  1.39

``` r
doc_data %>% 
  summarize_at(vars(DOC), list(mean = mean, sd = sd))
```

    ## # A tibble: 1 × 2
    ##    mean    sd
    ##   <dbl> <dbl>
    ## 1  66.7  1.45

diff = ~1.5 µM C between ash-impacted and not ash-impacted

“Surface water DBC concentrations observed beneath the Thomas Fire smoke
plume were similar to those observed in the Gulf of Mexico (1.2 μM;
Dittmar, 2008) and the Chukchi and Bering Seas (0.4–1.3 μM; Nakane et
al., 2017), but lower than DBC concentrations recorded for the Bohai and
Northern Yellow Seas (5.7–6.7 μM; Fang et al., 2021).”

SPE-DOC = 33.8 + 1.5µM ~49% recovery of DOC as SPE-DOC

DBC = 1 µM ratio = 0.03

~2 µM C
