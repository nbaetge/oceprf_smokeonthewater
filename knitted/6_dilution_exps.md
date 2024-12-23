Phytoplankton & zooplankton dilution experiments
================
Nick Baetge
compiled most recently on 01 June, 2024

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
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p6_dilution_exps.csv"
```

# Influx flow cytometry data

## Import and tidy data

``` r
influx_path <-
  "/Users/nicholasbaetge/github/oceprf_smokeonthewater/raw/r6_phyto.xlsx"
# data contains cells in units of per L 
```

``` r
ids <- readxl::read_xlsx(influx_path, sheet = 1) %>%
  select(exp, influx_id) %>% 
  distinct() 
```

``` r
context_path <- "/Users/nicholasbaetge/github/oceprf_smokeonthewater/prod/p1_insitu_index.csv"

context <- read_csv(context_path) %>% 
  select(stn, exp, lat, lon, biomass, composite_z)
```

``` r
apparent_gr  <- readxl::read_xlsx(influx_path, sheet = 2) %>%
  select(influx_id:nano_tf) %>%
  mutate(syn_100 = syn_t0[syn_dil == 1],
         pico_100 = pico_t0[pico_dil == 1],
         nano_100 = nano_t0[nano_dil == 1]) %>%
  group_by(influx_id, trt) %>%
  fill(c(syn_100:nano_100)) %>%
  select(-contains("dil")) %>%
  mutate(
    syn_d = 1 - (syn_t0 / syn_100),
    pico_d = 1 - (pico_t0 / pico_100),
    nano_d = 1 - (nano_t0 / nano_100),
    syn_k = log(syn_tf / syn_t0),
    pico_k = log(pico_tf / pico_t0),
    nano_k = log(nano_tf / nano_t0)
  ) %>%
  select(influx_id, trt, contains(c("_d", "_k"))) %>%
  ungroup()

rates <- apparent_gr %>%
  select(influx_id:syn_d, syn_k) %>%
  rename(d = syn_d,
         k = syn_k) %>%
  mutate(phyto = "syn", .after = "trt") %>%
  bind_rows(
    .,
    apparent_gr %>%
      select(influx_id:trt, pico_d, pico_k) %>%
      rename(d = pico_d,
             k = pico_k) %>%
      mutate(phyto = "pico", .after = "trt")
  ) %>%
  bind_rows(
    .,
    apparent_gr %>%
      select(influx_id:trt, nano_d, nano_k) %>%
      rename(d = nano_d,
             k = nano_k) %>%
      mutate(phyto = "nano", .after = "trt")
  ) %>%
  group_by(influx_id, trt, phyto) %>%
  mutate(
    g = ifelse(d != 0, (k - k[d == 0]) / d, NA),
    theory_violated = ifelse(g < 0, T, F),
    g_corr = ifelse(theory_violated == T, 0, g),
    ave_g = mean(g_corr, na.rm = T),
    sd_g = sd(g_corr, na.rm = T),
    mu = g_corr + k[d == 0],
    ave_mu = mean(mu, na.rm = T),
    sd_mu = sd(mu, na.rm = T),
    r = mu - g_corr,
    ave_r = mean(r, na.rm = T),
    sd_r = sd(r, na.rm = T)
  ) %>%
  ungroup() %>%
  left_join(ids, .) %>%
  left_join(., context) %>%
  mutate(amend = ifelse(trt == "CTL", "Control", "Ash leachate"),
         .after = trt) %>%
  mutate(
    ash = case_when(
      trt == "CTL" ~ "Control",
      trt == "TFA" ~ "Thomas Fire Ash",
      trt == "LTA" ~ "Lab Ash",
      trt == "HTA" ~ "Lab Ash"
    ),
    trt = case_when(
      trt == "CTL" ~ "Control",
      trt == "TFA" ~ "Thomas Fire Ash",
      trt == "LTA" ~ "Low Temp. Ash",
      trt == "HTA" ~ "High Temp. Ash"
    ),
    phyto = case_when(
      phyto == "syn" ~ "italic(Synechococcus)",
      phyto == "pico" ~ "Picoeukaryotes",
      phyto == "nano" ~ "Nanoeukaryotes"
    )
  ) %>%
  mutate(
    biomass = ifelse(
      biomass == "Lower",
      "bold(Lower~biomass~index)",
      "bold(Higher~biomass~index)"
    )
  ) %>%
  select(exp:trt, ash, everything())
```

# Plots

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

pal4 = c("#2C1B21", "#6aa42c", "#5990e7", "#b95f83")
pal3 = c("#6aa42c", "#5990e7", "#b95f83")
pal2 = c("#2C1B21", "white")
plot_levels = c(
  "Control",
  "Ash leachate",
  "Thomas Fire Ash",
  "Lab Ash",
  "Low Temp. Ash",
  "High Temp. Ash",
  "DOC",
  "bold(Leachate)",
  "bold(Amendment)",
  "bold(Lower~biomass~index)",
  "bold(Biomass~index~'<'~-1)",
  "bold(Higher~biomass~index)",
  "bold(-1~'<'~Biomass~index~'<'~1)",
  "bold(Biomass~index~'>'~1)",
  "bold(Division~(d^-1))",
  "bold(Grazing~(d^-1))",
  "bold(Accumulation~(d^-1))",
  "P13",
  "PB16",
  "PB15",
  "P17",
  "P1.5",
  "P12",
  "PB1",
  "PB14",
  "PB2",
  "PB3",
  "P2.5",
  "P4",
  "P3.5",
  "P5",
  "italic(Synechococcus)",
  "Picoeukaryotes",
  "Nanoeukaryotes"
)
```

``` r
p_rates <- rates %>%
  select(exp:phyto, biomass, composite_z, mu, g_corr, r) %>%
  pivot_longer(c(mu, g_corr, r), names_to = "rate", values_to = "val") %>%
  mutate(
    rate = case_when(
      rate == "mu" ~ "bold(Division~(d^-1))",
      rate == "g_corr" ~ "bold(Grazing~(d^-1))",
      rate == "r" ~ "bold(Accumulation~(d^-1))"
    )
  ) %>%
  distinct() %>% 
 mutate(biomass2 = case_when(composite_z <= -1 ~ "bold(Biomass~index~'<'~-1)",
                              composite_z > -1 & composite_z < 1 ~ "bold(-1~'<'~Biomass~index~'<'~1)",
                              composite_z >= 1 ~ "bold(Biomass~index~'>'~1)")) 

p_rates_summary <- rates %>%
  select(exp:phyto, biomass, composite_z, ave_mu, ave_g, ave_r) %>%
  distinct() %>%
  pivot_longer(c(ave_mu, ave_g, ave_r),
               names_to = "rate",
               values_to = "val") %>%
  mutate(
    rate = case_when(
      rate == "ave_mu" ~ "bold(Division~(d^-1))",
      rate == "ave_g" ~ "bold(Grazing~(d^-1))",
      rate == "ave_r" ~ "bold(Accumulation~(d^-1))"
    )
  ) %>%
  distinct() %>%
  left_join(
    .,
    rates %>%
      select(exp:phyto, biomass, composite_z, sd_mu, sd_g, sd_r) %>%
      distinct() %>%
      pivot_longer(
        c(sd_mu, sd_g, sd_r),
        names_to = "rate",
        values_to = "sd"
      ) %>%
      mutate(
        rate = case_when(
          rate == "sd_mu" ~ "bold(Division~(d^-1))",
          rate == "sd_g" ~ "bold(Grazing~(d^-1))",
          rate == "sd_r" ~ "bold(Accumulation~(d^-1))"
        )
      )
  ) %>% 
  mutate(biomass2 = case_when(composite_z <= -1 ~ "bold(Biomass~index~'<'~-1)",
                              composite_z > -1 & composite_z < 1 ~ "bold(-1~'<'~Biomass~index~'<'~1)",
                              composite_z >= 1 ~ "bold(Biomass~index~'>'~1)")) %>% 
  left_join(., p_rates %>% select(phyto, rate) %>% distinct() %>%  mutate(flabel = 1:n()))
```

``` r
dots <- ggplot(p_rates_summary ,
               aes(x = composite_z,
                   y = val,)) +
  geom_errorbar(
    aes(
      color = factor(trt, levels = plot_levels),
      ymin = val - sd,
      ymax = val + sd
    ),
    width = .2
    # position = position_dodge(.3)
  ) +
  stat_smooth(
    aes(fill = factor(amend, levels = plot_levels),
        linetype = factor(amend, levels = plot_levels)),
    color = "black", 
    method = "gam",
    formula = y ~ s(x, bs = "tp"),
    alpha = 0.25,
    size = 2
  ) +
  geom_point(
    aes(color = factor(trt, levels = plot_levels)),
    # position = position_dodge(width = 0.3),
    size = 10,
    alpha = 0.8
  ) +
  facet_grid(
    # ggh4x::facet_nested_wrap(
    factor(rate, levels = plot_levels) ~  factor(phyto, levels = plot_levels),
    labeller = label_parsed,
    scales = "free",
    # nrow = 2
  ) +
  scale_x_continuous(breaks = c(-6, -3, -1, 0 , 1, 3, 6)) +
  scale_fill_manual(values = pal4, guide = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
  scale_color_manual(values = pal4) +
  labs(
    y = expression(""),
    x = expression("Biomass index"),
    fill = "GAM",
    color = ""
  ) +
  theme_linedraw() +
  custom.theme +
  guides(linetype = "none") +
  geom_label(aes(x = 6, y = -0.5, label = flabel), size = 14)
```

``` r
boxes <-
  ggplot(p_rates, aes(x = factor(amend, levels = plot_levels), y = val)) +
  geom_boxplot(
    aes(fill = factor(amend, levels = plot_levels)),
    alpha = 0.7,
    outlier.shape = NA,
    width = 0.5,
    position = position_dodge(0.4)
  ) +
  geom_rect(data = subset(p_rates, biomass2 == "bold(-1~'<'~Biomass~index~'<'~1)"), fill = "light grey",xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.01) +
  geom_jitter(
    aes(color = factor(trt, levels = plot_levels)),
    size = 7,
    position = position_jitterdodge(),
    alpha = 0.5
  ) +
  ggh4x::facet_nested(
    factor(rate, levels = plot_levels) ~  factor(biomass2, levels = plot_levels) + factor(phyto, levels = plot_levels),
    labeller = label_parsed,
    scales = "free"
  ) +
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal4) +
  labs(
    y = expression(""),
    x = expression(""),
    fill = "",
    color = ""
  ) +
  theme_linedraw() +
  custom.theme +
  guides(fill = "none") +
  # guides(color= "none",fill = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.text.x = element_blank()) +
  ggpubr::geom_pwc(
    method = "wilcoxon",
    label = "{p.format}{p.signif}",
    hide.ns = TRUE,
    vjust = 0.2,
    size = 0.6,
    label.size = 7
  )
```

``` r
(
  dots
) /(boxes + guides(color = "none", fill = "none")) + plot_layout(heights = c(1,1)) + plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 45))
```

![](/Users/nicholasbaetge/github/oceprf_smokeonthewater/knitted/6_dilution_exps_files/figure-gfm/Figure4-1.png)<!-- -->
