
library(tidyverse)
library(rstanarm)
source("lib/misc/disturbance_summary.R")
options(mc.cores = 1)
library(patchwork)
library(sf)

getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
colors_regions <- c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377")
colors_scenarios <- c("#BB5566", "#004488", "#DDAA33")
sim_years <- 2019:2050

# Functions ---------------------------------------------------------------

posterior_check <- function(fit, ref = dat_summary_subset) {
  
  post_pred <- posterior_predict(fit) %>%
    as.data.frame() %>%
    mutate(iter = 1:n()) %>%
    gather(key = year, value = disturbance, -iter) %>%
    group_by(iter) %>%
    summarize(mean = mean(disturbance, na.rm = TRUE),
              max = max(disturbance, na.rm = TRUE),
              min = min(disturbance, na.rm = TRUE),
              sd = sd(disturbance, na.rm = TRUE),
              median = median(disturbance, na.rm = TRUE),
              iqr = IQR(disturbance, na.rm = TRUE)) %>%
    ungroup() %>%
    gather(key = metric, value = value, -iter)
  
  true <- ref %>%
    summarize(mean = mean(disturbance, na.rm = TRUE),
              max = max(disturbance, na.rm = TRUE),
              min = min(disturbance, na.rm = TRUE),
              sd = sd(disturbance, na.rm = TRUE),
              median = median(disturbance, na.rm = TRUE),
              iqr = IQR(disturbance, na.rm = TRUE)) %>%
    gather(key = metric, value = value_true)
  
  ggplot(post_pred, aes(x = value)) +
    geom_histogram() +
    geom_vline(data = true, aes(xintercept = value_true), col = "red") +
    facet_wrap(~metric, scales = "free") +
    theme_gray()
  
}

# Get data ----------------------------------------------------------------

names <- list.files("data/timesync/databases_export/", "export.csv") %>%
  strsplit(., "_") %>%
  map(~ .[2]) %>%
  unlist()

country_information <- read_csv("data/countries_master.csv")
country_grouping <- read_delim("data/strata/country_grouping.csv", delim = ";")

weights <- readxl::read_excel("rawdata/forest_area.xlsx") %>%
  mutate(weight_country = forest_area_km2 / sum(forest_area_km2)) %>%
  dplyr::select(country = country_name_short, weight_country, forest_area_km2) %>%
  left_join(country_grouping, by = c("country" = "country_name_short")) %>%
  group_by(euro_region) %>%
  mutate(weight_region = forest_area_km2 / sum(forest_area_km2)) %>%
  ungroup() %>%
  dplyr::select(country, weight_country, weight_region)

country_information$end_year <- 2018
ceur <- c("austria", "czechia", "germany", "poland", "slovakia", "switzerland")
country_information[country_information$country_name_short %in% ceur, "end_year"] <- 2016

dat_summary_all <- list.files("data/timesync/databases_export/", "export.csv", full.names = TRUE) %>%
  map(read_csv) %>%
  map2(.y = arrange(country_information, country_name_short)[country_information$country_name_short %in% names, "end_year"][[1]], 
       ~ disturbance_summary2(.x, 
                              grouping = "year",
                              grouping_range = list(year = 1985:as.integer(.y)))) %>%
  set_names(names) %>%
  bind_rows(.id = "country")

### Repalce first entry of austria to exclude long-term declines, which were still recorded in the first assessment

dat_summary_all[dat_summary_all$country == "austria" & dat_summary_all$year == 1985, "disturbance"] <- NA

# Fit individual trend models ---------------------------------------------

trend <- vector("list", length(names))
trendline <- vector("list", length(names))
annual_rates <- vector("list", length(names))
simulations <- vector("list", length(names))
intercept <- vector("list", length(names))
slope <- vector("list", length(names))

for (k in 1:length(names)) {
  
  cntr <- names[k]
  #cntr <- "estonia"
  
  dat_summary_subset <- dat_summary_all %>%
    filter(country == cntr)
  
  dat_summary_subset <- dat_summary_subset %>%
    mutate(stable = forest - disturbance) %>%
    na.omit(.)
  
  fit_bin <- stan_glmer(cbind(disturbance, stable) ~ year + (1 | year),
                        data = dat_summary_subset,
                        family = binomial(link = "logit"),
                        chains = 4, 
                        iter = 2000)
  
  posterior_checks <- posterior_check(fit_bin)
  ggsave(paste0("results/posterior_checks/posterior-check_", cntr, ".pdf"), posterior_checks, width = 7.5, height = 5)
  
  ### Extract posterior of trend parameter
  
  posterior <- as.matrix(fit_bin)
  
  intercept[[k]] <- posterior[, 1]
  slope[[k]] <- posterior[, 2]
  
  trend[[k]] <- posterior %>%
    as.data.frame() %>%
    gather() %>%
    filter(key == "year") %>%
    mutate(iter = paste0("iter_", 1:4000)) %>%
    mutate(trend = (exp(value) - 1) * 100) %>%
    dplyr::select(-key, -value)
  
  # ggplot(trend, aes(x = trend)) +
  #   geom_histogram()
  
  ### Extract posterior of trendline
  
  pred_trendline <- posterior_linpred(fit_bin, 
                                      newdata = data.frame(year = 1986:2018), 
                                      draws = 4000,
                                      transform = TRUE,
                                      re.form = NA)
  
  trendline[[k]] <- pred_trendline %>% 
    as.data.frame() %>%
    mutate(iter = paste0("iter_", 1:4000)) %>%
    gather(key = year, value = disturbance_rate, -iter) %>%
    mutate(year = rep(1986:2018, each = 4000)) %>%
    mutate(disturbance_rate = disturbance_rate * 100)
  
  # ggplot(trendline[[k]] %>%
  #          filter(iter %in% paste0("iter_", 1:150)),
  #        aes(x = year, y = disturbance_rate)) +
  #   geom_line(aes(group = iter), alpha = 0.2)
  
  ### Extract posterior of annual disturbance rates
  
  pred_rates <- posterior_predict(fit_bin,
                                  newdata = data.frame(year = 1986:2018,
                                                       disturbance = 0,
                                                       stable = unique(dat_summary_subset$forest)),
                                  draws = 4000,
                                  transform = TRUE)
  
  pred_rates <- pred_rates / unique(dat_summary_subset$forest) * 100
  
  annual_rates[[k]] <- pred_rates %>% 
    as.data.frame() %>%
    mutate(iter = paste0("iter_", 1:4000)) %>%
    gather(key = year, value = disturbance_rate, -iter) %>%
    mutate(year = rep(1986:2018, each = 4000))
  
  # ggplot(annual_rates[[k]] %>%
  #          group_by(year) %>%
  #          summarize(est = median(disturbance_rate),
  #                    lower = quantile(disturbance_rate, 0.2),
  #                    upper = quantile(disturbance_rate, 0.8)),
  #        aes(x = year, y = est)) +
  #   geom_point(alpha = 0.2) +
  #   geom_errorbar(aes(ymin = lower, ymax = upper))
  
  ### Simulate future
  
  pred_sim <- posterior_predict(fit_bin, 
                                draws = 4000,
                                newdata = data.frame(year = sim_years,
                                                     disturbance = 0,
                                                     stable = unique(dat_summary_subset$forest)),
                                transform = TRUE)
  
  pred_sim <- pred_sim / unique(dat_summary_subset$forest) * 100
  
  simulations[[k]] <- pred_sim %>% 
    as.data.frame() %>%
    mutate(iter = paste0("iter_", 1:4000)) %>%
    gather(key = year, value = disturbance_rate, -iter) %>%
    mutate(year = rep(sim_years, each = 4000))
  
  # ggplot(simulations[[k]] %>%
  #          filter(iter %in% paste0("iter_", 1:30)),
  #        aes(x = year, y = disturbance_rate)) +
  #   geom_point(alpha = 0.2)
  
}

intercept <- intercept %>%
  map(., ~ data.frame(intercept = .,
                      iter = paste0("iter_", 1:4000))) %>%
  set_names(names) %>%
  bind_rows(.id = "country")

save(intercept, file = "results/intercept.RData")

slope <- slope %>%
  map(., ~ data.frame(slope = .,
                      iter = paste0("iter_", 1:4000))) %>%
  set_names(names) %>%
  bind_rows(.id = "country")

save(slope, file = "results/slope.RData")

trend <- trend %>%
  set_names(names) %>%
  bind_rows(.id = "country")

save(trend, file = "results/trend.RData")

trend <- trend %>%
  set_names(names) %>%
  bind_rows(.id = "country")

save(trend, file = "results/trend.RData")

trendline <- trendline %>%
  set_names(names) %>%
  bind_rows(.id = "country")

save(trendline, file = "results/trendline.RData")

annual_rates <- annual_rates %>%
  set_names(names) %>%
  bind_rows(.id = "country")

save(annual_rates, file = "results/annual_rates.RData")

simulations <- simulations %>%
  set_names(names) %>%
  bind_rows(.id = "country")

save(simulations, file = "results/simulations.RData")

# Calculate summary statistics ---------------------------------------------------

load(file = "results/simulations.RData")
load(file = "results/annual_rates.RData")
load(file = "results/trendline.RData")
load(file = "results/trend.RData")

trend <- trend %>% 
  left_join(country_grouping, by = c("country" = "country_name_short")) %>%
  left_join(weights, by = "country")

trendline <- trendline %>% 
  left_join(country_grouping, by = c("country" = "country_name_short")) %>%
  left_join(weights, by = "country")

annual_rates <- annual_rates %>% 
  left_join(country_grouping, by = c("country" = "country_name_short")) %>%
  left_join(weights, by = "country")

### Country

# Disturbance trend

trend_country_summary <- trend %>% 
  group_by(country) %>%
  summarize(median = median(trend), 
            q05 = quantile(trend, 0.05), 
            q15 = quantile(trend, 0.15), 
            q25 = quantile(trend, 0.25), 
            q75 = quantile(trend, 0.25), 
            q85 = quantile(trend, 0.85), 
            q95 = quantile(trend, 0.95),
            mean = mean(trend),
            sd = sd(trend))

trend_country_summary %>% write_csv("results/trend_country_summary.csv")

# Disturbance trend-line

trendline_country_summary <- trendline %>% 
  group_by(country, year) %>% 
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

# Annual average disturbance rates

annual_rates_country_summary <- annual_rates %>%
  group_by(country, year, iter) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(country, year) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

write_csv(annual_rates_country_summary, "results/annual_rates_country_summary.csv")

annual_rates_country_summary_periods <- annual_rates %>%
  mutate(period = case_when(
    year %in% c(1985:1999) ~ "early",
    year %in% c(2000:2018) ~ "late",
  )) %>%
  group_by(country, period, iter) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(country, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

annual_rates_country_summary_average <- annual_rates %>%
  mutate(period = case_when(
    year %in% c(1985:2018) ~ "reference"
  )) %>%
  group_by(country, period, iter) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(country, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

list(annual_rates_country_summary_periods,
     annual_rates_country_summary_average) %>%
  bind_rows() %>%
  write_csv("results/period_rates_country_summary.csv")

# Annual maximum disturbance rates

annual_rates_max_country_summary <- annual_rates %>%
  group_by(country, year, iter) %>%
  summarize(disturbance_rate = quantile(disturbance_rate, 0.99)) %>%
  group_by(country, year) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

annual_rates_max_country_summary_periods <- annual_rates %>%
  mutate(period = case_when(
    year %in% c(1985:1999) ~ "early",
    year %in% c(2000:2018) ~ "late",
  )) %>%
  group_by(country, period, iter) %>%
  summarize(disturbance_rate = quantile(disturbance_rate, 0.99)) %>%
  group_by(country, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

annual_rates_max_country_summary_average <- annual_rates %>%
  mutate(period = case_when(
    year %in% c(1985:2018) ~ "reference"
  )) %>%
  group_by(country, period, iter) %>%
  summarize(disturbance_rate = quantile(disturbance_rate, 0.99)) %>%
  group_by(country, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

list(annual_rates_max_country_summary_periods,
     annual_rates_max_country_summary_average) %>%
  bind_rows() %>%
  write_csv("results/period_rates_max_country_summary.csv")

# Simulations

simulations_country_summary_periods <- simulations %>%
  mutate(period = case_when(
    year %in% c(2019:2034) ~ "sim_early",
    year %in% c(2035:2050) ~ "sim_late"
  )) %>%
  group_by(country, period, iter) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(country, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

simulations_country_summary_periods %>% write_csv("results/simulations_country_summary.csv")

### Regions

# Disturbance trend

trend_region <- trend %>%
  group_by(region = euro_region, iter) %>%
  summarize(trend = sum(trend * weight_region))

trend_region_summary <- trend_region %>% 
  group_by(region) %>% 
  summarize(median = median(trend), 
            q05 = quantile(trend, 0.05), 
            q15 = quantile(trend, 0.15), 
            q25 = quantile(trend, 0.25), 
            q75 = quantile(trend, 0.25), 
            q85 = quantile(trend, 0.85), 
            q95 = quantile(trend, 0.95),
            mean = mean(trend),
            sd = sd(trend))

trend_region_summary %>% write_csv("results/trend_region_summary.csv")

# Disturbance trend-line

trendline_region <- trendline %>%
  group_by(region = euro_region, year, iter) %>%
  summarize(disturbance_rate = sum(disturbance_rate * weight_region))

trendline_region_summary <- trendline_region %>% 
  group_by(region, year) %>% 
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

# Annual average disturbance rates

annual_rates_region <- annual_rates %>%
  group_by(region = euro_region, year, iter) %>%
  summarize(disturbance_rate = sum(disturbance_rate * weight_region))

annual_rates_region_summary <- annual_rates_region %>% 
  group_by(region, year) %>% 
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

write_csv(annual_rates_region_summary, "results/annual_rates_region_summary.csv")

annual_rates_region_summary_periods <- annual_rates_region %>%
  mutate(period = case_when(
    year %in% c(1985:1999) ~ "early",
    year %in% c(2000:2018) ~ "late",
  )) %>%
  group_by(region, period, iter) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(region, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

annual_rates_region_summary_average <- annual_rates_region %>%
  mutate(period = case_when(
    year %in% c(1985:2018) ~ "reference"
  )) %>%
  group_by(region, period, iter) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(region, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

list(annual_rates_region_summary_periods,
     annual_rates_region_summary_average) %>%
  bind_rows() %>%
  write_csv("results/period_rates_region_summary.csv")

# Annual maximum disturbance rates

annual_rates_region <- annual_rates %>%
  group_by(region = euro_region, year, iter) %>%
  summarize(disturbance_rate = sum(disturbance_rate * weight_region))

annual_rates_max_region_summary <- annual_rates_region %>% 
  group_by(iter, region, year) %>%
  summarize(disturbance_rate = quantile(disturbance_rate, 0.99)) %>%
  group_by(region, year) %>% 
  summarize(median = quantile(disturbance_rate, 0.5), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

write_csv(annual_rates_max_region_summary, "period_rates_max_country_summary.csv")

annual_rates_max_region_summary_periods <- annual_rates_region %>%
  mutate(period = case_when(
    year %in% c(1985:1999) ~ "early",
    year %in% c(2000:2018) ~ "late",
  )) %>%
  group_by(region, period, iter) %>%
  summarize(disturbance_rate = quantile(disturbance_rate, 0.99)) %>%
  group_by(region, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

annual_rates_max_region_summary_average <- annual_rates_region %>%
  mutate(period = case_when(
    year %in% c(1985:2018) ~ "reference"
  )) %>%
  group_by(region, period, iter) %>%
  summarize(disturbance_rate = quantile(disturbance_rate, 0.99)) %>%
  group_by(region, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

list(annual_rates_max_region_summary_periods,
     annual_rates_max_region_summary_average) %>%
  bind_rows() %>%
  write_csv("results/period_rates_max_region_summary.csv")

### Europe

# Disturbance trends

trend_europe <- trend %>%
  group_by(iter) %>%
  summarize(trend = sum(trend * weight_country)) %>%
  mutate(region = "europe")

trend_europe_summary <- trend_europe %>%
  group_by(region) %>%
  summarize(median = median(trend), 
            q05 = quantile(trend, 0.05), 
            q15 = quantile(trend, 0.15), 
            q25 = quantile(trend, 0.25), 
            q75 = quantile(trend, 0.25), 
            q85 = quantile(trend, 0.85), 
            q95 = quantile(trend, 0.95),
            mean = mean(trend),
            sd = sd(trend))

trend_europe_summary %>% write_csv("results/trend_europe_summary.csv")

# Disturbance trend lines

trendline_europe <- trendline %>%
  group_by(year, iter) %>%
  summarize(disturbance_rate = sum(disturbance_rate * weight_country)) %>%
  mutate(region = "europe")

trendline_europe_summary <- trendline_europe %>% 
  group_by(region, year) %>% 
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

# Annual disturbance rates

annual_rates_europe <- annual_rates %>%
  group_by(year, iter) %>%
  summarize(disturbance_rate = sum(disturbance_rate * weight_country)) %>%
  mutate(region = "europe")

annual_rates_europe_summary <- annual_rates_europe %>% 
  group_by(region, year) %>% 
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.75), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

write_csv(annual_rates_europe_summary, "results/annual_rates_europe_summary.csv")

annual_rates_europe_summary_periods <- annual_rates_europe %>%
  mutate(period = case_when(
    year %in% c(1985:1999) ~ "early",
    year %in% c(2000:2018) ~ "late",
  )) %>%
  group_by(region, period, iter) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(region, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.75), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

annual_rates_europe_summary_average <- annual_rates_europe %>%
  mutate(period = case_when(
    year %in% c(1985:2018) ~ "reference"
  )) %>%
  group_by(region, period, iter) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(region, period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.75), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))

list(annual_rates_europe_summary_periods,
     annual_rates_europe_summary_average) %>%
  bind_rows() %>%
  write_csv("results/period_rates_europe_summary.csv")

# Simulations

simulations_europe_summary <- simulations %>%
  left_join(weights, by = "country") %>%
  left_join(country_grouping, by = c("country" = "country_name_short")) %>%
  group_by(year, iter) %>%
  summarize(disturbance_rate = sum(weight_country * disturbance_rate)) %>%
  ungroup(.)

simulations_europe_summary %>% write_csv("results/simulations_europe_summary.csv")

simulations_europe_summary_period <- simulations_europe_summary %>%
  mutate(period = case_when(
    year %in% c(2019:2035) ~ "early",
    year %in% c(2036:2050) ~ "late"
  )) %>%
  group_by(iter, period) %>%
  summarize(disturbance_rate = mean(disturbance_rate)) %>%
  group_by(period) %>%
  summarize(median = median(disturbance_rate), 
            q05 = quantile(disturbance_rate, 0.05), 
            q15 = quantile(disturbance_rate, 0.15), 
            q25 = quantile(disturbance_rate, 0.25), 
            q75 = quantile(disturbance_rate, 0.25), 
            q85 = quantile(disturbance_rate, 0.85), 
            q95 = quantile(disturbance_rate, 0.95),
            mean = mean(disturbance_rate),
            sd = sd(disturbance_rate))


# Intercept and slope (for plotting)

intercept_slope_europe <- intercept %>%
  left_join(slope) %>%
  left_join(weights) %>% 
  group_by(iter) %>%
  summarize(intercept = sum(intercept * weight_country),
            slope = sum(slope * weight_country)) %>%
  mutate(region = "europe")

intercept_slope_europe_summary <- intercept_slope_europe %>%
  gather(key = parameter, value = value, -iter, -region) %>%
  group_by(region, parameter) %>% 
  summarize(median = median(value), 
            q05 = quantile(value, 0.05), 
            q15 = quantile(value, 0.15), 
            q25 = quantile(value, 0.25), 
            q75 = quantile(value, 0.75), 
            q85 = quantile(value, 0.85), 
            q95 = quantile(value, 0.95),
            mean = mean(value),
            sd = sd(value))

# Figures -----------------------------------------------------------------

country_name_changer <- function(x) {
  x <- str_to_title(x)
  x <- ifelse(x == "Bosniaherzegovina", "Bosnia and Herzegovina", x)
  x <- ifelse(x == "Macedonia", "North Macedonia", x)
  x <- ifelse(x == "Unitedkingdom", "United Kingdom", x)
  return(x)
}

### Country disturbance rates and trends

p <- ggplot() +
  # geom_errorbar(data = annual_rates_country_summary%>%
  #                 ungroup() %>%
  #                 mutate(country = country_name_changer(country)),
  #               aes(x = year, ymin = mean - sd, ymax = mean + sd),
  #               col = "grey", width = 0, size = 1.25) +
  # geom_errorbar(data = annual_rates_country_summary%>%
  #                 ungroup() %>%
  #                 mutate(country = country_name_changer(country)),
  #               aes(x = year, ymin = q05, ymax = q95),
  #               col = "grey", width = 0) +
  geom_ribbon(data = annual_rates_country_summary%>%
                  ungroup() %>%
                  mutate(country = country_name_changer(country)),
                aes(x = year, ymin = q05, ymax = q95),
                fill = "#DDDDDD") +
  geom_ribbon(data = annual_rates_country_summary%>%
                ungroup() %>%
                mutate(country = country_name_changer(country)),
              aes(x = year, ymin = mean - sd, ymax = mean + sd),
              fill = "#BBBBBB") +
  geom_point(data = annual_rates_country_summary %>%
               ungroup() %>%
               mutate(country = country_name_changer(country)),
             aes(x = year, y = mean),
             shape = 21, col = "#555555") +
  geom_line(data = trendline_country_summary %>%
              ungroup() %>%
              mutate(country = country_name_changer(country)), aes(x = year, y = mean), col = "#BB5566") + 
  geom_ribbon(data = trendline_country_summary %>%
                ungroup() %>%
                mutate(country = country_name_changer(country)), 
              aes(x = year, ymin = mean - sd, ymax = mean + sd), alpha = 0.3, fill = "#BB5566") +
  theme_minimal() +
  scale_color_brewer(palette = "Paired") +
  labs(x = "Year", 
       y = bquote("Canopy mortality rate (% "*yr^-1*")"),
       col = "Deviance from average (1985-2014)") +
  scale_x_continuous(breaks = c(1985, 2000, 2015)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(size = 6)
  ) +
  facet_wrap(~ country, ncol = 7, scales = "free_y")

ggsave("results/disturbance_trends_country.pdf", p, width = 7.5, height = 6)

### Europe disturbance rates + trends

p <- ggplot() +
  # geom_point(data = annual_rates_europe %>% filter(iter %in% paste0("iter_", 1:150)),
  #            aes(x = year, y = disturbance_rate),
  #            col = "grey",
  #            position = position_jitter(width = 0.075), 
  #            alpha = 0.1, shape = 19, stroke = NA) +
  geom_ribbon(data = annual_rates_europe_summary,
              aes(x = year, ymin = q05, ymax = q95),
              fill = "#DDDDDD") +
  geom_ribbon(data = annual_rates_europe_summary,
              aes(x = year, ymin = mean - sd, ymax = mean + sd),
              fill = "#BBBBBB") +
  # geom_errorbar(data = annual_rates_europe_summary,
  #               aes(x = year, ymin = mean - sd, ymax = mean + sd),
  #               col = "grey", width = 0, size = 1.5) +
  # geom_errorbar(data = annual_rates_europe_summary,
  #               aes(x = year, ymin = q05, ymax = q95),
  #               col = "grey", width = 0) +
  geom_point(data = annual_rates_europe_summary,
             aes(x = year, y = median),
             shape = 21, col = "#555555") +
  # geom_segment(data = annual_rates_europe_summary_periods %>%
  #                mutate(year_start = case_when(period == "early" ~ 1985, period == "late" ~ 2000), 
  #                       year_end = case_when(period == "early" ~ 1999, period == "late" ~ 2018)), 
  #              aes(y = median, yend = median, 
  #                  x = year_start, 
  #                  xend = year_end,
  #                  group = period), 
  #              col = "grey20", linetype = "dashed", alpha = 0.75) +
  geom_line(data = trendline_europe_summary, aes(x = year, y = mean), col = "#BB5566") + 
  geom_ribbon(data = trendline_europe_summary, 
              aes(x = year, ymin = mean - sd, ymax = mean + sd), alpha = 0.3, fill = "#BB5566") +
  theme_minimal() +
  labs(x = "Year", 
       y = bquote("Canopy mortality rate (% "*yr^-1*")"),
       col = "Deviance from average (1986-2018)") +
  ylim(0.5, 1.5) +
  scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010, 2015)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA)) +
  annotate("text", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.1, size = 3,
           label = paste0("Average: ", round(annual_rates_europe_summary_average$mean, 2), " ± ", round(annual_rates_europe_summary_average$sd, 2), " % per year\n",
                          "Trend: ", round(trend_europe_summary$mean, 2), " ± ", round(trend_europe_summary$sd, 2), " % per year"))
  # geom_text(data = annual_rates_europe_summary_periods %>% 
  #             mutate(xpos = case_when(period == "early" ~ 2002, TRUE ~ 1997),
  #                    ypos = case_when(period == "early" ~ 0.725, TRUE ~ 1.075)), 
  #           aes(x = xpos, y = ypos, 
  #               label = paste0("", round(median, 2), "~'%'~yr^-1")),
  #           size = 2.5, parse = TRUE, hjust = 0.5, col = "grey20")

ggsave("results/disturbance_trend_europe.pdf", p, width = 3.5, height = 3.5)
ggsave("results/disturbance_trend_europe_slides.pdf", p + theme(axis.title = element_text(size = 12),
                                                                axis.text.x = element_text(size = 10),
                                                                axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5)), 
       width = 5.5, height = 5.5)

### Trend parameters region

p <- ggplot() +
  geom_errorbar(data = trend_region_summary,
                aes(x = reorder(str_to_title(region), mean), 
                    ymin = mean - sd, ymax = mean + sd, 
                    col = str_to_title(region)), 
                width = 0) +
  geom_point(data = trend_region_summary,
             aes(x = reorder(str_to_title(region), mean), y = mean, 
                 col = str_to_title(region),
                 shape = str_to_title(region))) +
  theme_minimal() +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "darkgrey") + 
  scale_color_manual(values = colors_regions) +
  labs(x = NULL, 
       y = bquote("Trend (% "*yr^-1*")"),
       col = NULL) +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, color = "black"), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank())

ggsave("results/disturbance_trend_parameter_region.pdf", p, width = 1.5, height = 1.5)

# Trend parameter country

p <- ggplot() +
  geom_errorbar(data = trend_country_summary %>%
                  left_join(country_grouping, by = c("country" = "country_name_short")),
                aes(x = reorder(country_name_changer(country), mean), 
                    ymin = mean - sd, ymax = mean + sd, 
                    col = str_to_title(euro_region)), 
                width = 0) +
  geom_point(data = trend_country_summary %>%
               left_join(country_grouping, by = c("country" = "country_name_short")),
             aes(x = reorder(country_name_changer(country), mean), y = mean, 
                 col = str_to_title(euro_region),
                 shape = str_to_title(euro_region))) +
  theme_minimal() +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "darkgrey") + 
  scale_color_manual(values = colors_regions) +
  labs(x = NULL, 
       y = bquote("Trend in canopy mortality (% "*yr^-1*")"),
       col = NULL) +
  ylim(-7, 14) +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, color = "black"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank())

ggsave("results/disturbance_trend_parameter_country.pdf", p, width = 3.5, height = 3.5)

### Correlation between rate and trend

trend_vs_average <- trend_country_summary %>%
  dplyr::select(country, trend_mn = mean, trend_sd = sd) %>%
  left_join(annual_rates_country_summary_average %>%
              dplyr::select(country, rate_mn = mean, rate_sd = sd), 
            by = "country") %>%
  left_join(country_grouping, by = c("country" = "country_name_short"))

p <- ggplot(trend_vs_average, aes(x = trend_mn, y = rate_mn, col = euro_region)) +
  ggrepel::geom_text_repel(data = trend_vs_average %>% filter(country %in% c("slovenia", "portugal")),
                           aes(label = str_to_title(country)),
                           size = 1.5, box.padding = 0.75, segment.colour = "grey", col = "grey") +
  geom_point() +
  geom_errorbarh(aes(xmin = trend_mn - trend_sd, xmax = trend_mn + trend_sd)) +
  geom_errorbar(aes(ymin = rate_mn - rate_sd, ymax = rate_mn + rate_sd)) +
  theme_minimal() +
  scale_color_brewer(palette = "Paired") +
  labs(x = bquote("Trend in canopy mortality (% "*yr^-1*")"),
       y = bquote("Average canopy mortality (% "*yr^-1*")"),
       col = NULL) +
  theme(legend.position = "none",
        axis.title = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 5, color = "black"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank())

ggsave("results/trend_vs_average_country.pdf", p, width = 1.75, height = 1.75)

### Map of regions

regions_map <- read_sf("rawdata/admin/europe/countries_europe.shp") %>%
  left_join(country_grouping, by = c("ISO_CC" = "iso_code")) %>%
  mutate(label = str_to_title(euro_region))

p <- ggplot() +
  geom_sf(data = regions_map, aes(fill = euro_region), col = NA) +
  scale_fill_manual(values = colors_regions) +
  theme_void() +
  theme(legend.position = "none")

ggsave("results/region_map.pdf", p, width = 2.25, height = 3.5 - 2.25)

p <- ggplot() +
  geom_sf(data = regions_map, aes(fill = str_to_title(euro_region)), col = "black") +
  scale_fill_manual(values = colors_regions) +
  theme_void() +
  labs(fill = NULL)

ggsave("results/region_map_large.pdf", p, width = 7.5, height = 5.5)

### Map of trends

countries_shp <- read_sf("rawdata/admin/europe/countries_europe.shp")

countries_shp <- countries_shp %>%
  left_join(country_information, by = c("COUNTRY" = "country_name")) %>%
  left_join(trend_country_summary, by = c("country_name_short" = "country"))

getPalette = colorRampPalette(RColorBrewer::brewer.pal(11, "PiYG"))

pdf("results/trend_map.pdf", width = 4, height = 4.5)
plot(countries_shp["median"],
     main = "",
     key.pos = 1,
     pal = getPalette(17)[12:1][-c(3:5)],
     breaks = c(-4, -2, 0, 2, 4, 6, 8, 10, 12, 14),
     bg = "white")
dev.off()

### Country-trends

cntr = "Ukraine"

p <- ggplot(annual_rates_country_summary %>% filter(country == tolower(cntr)), 
       aes(x = year, y = mean)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.25) + theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  labs(x = "Year", y = bquote("Canopy disturbance rate (% "*yr^-1*")"), title = cntr)

ggsave(filename = paste0("/Users/corneliussenf/Desktop/trend_", cntr, ".png"), plot = p, width = 3, height = 3)

### Future rates 

simulations_stabilize <- simulations %>%
  filter(year %in% 2019)

simulations_stabilize <- list(simulations_stabilize, 
                              simulations_stabilize, 
                              simulations_stabilize, 
                              simulations_stabilize, 
                              simulations_stabilize) %>%
  bind_rows() %>%
  mutate(year = sample(2019:2050, length(year), replace = TRUE))

p <- ggplot() +
  geom_ribbon(data = annual_rates %>%
                ungroup() %>%
                mutate(country = country_name_changer(str_to_title(country))) %>%
                group_by(country, year) %>%
                summarize(mean = mean(disturbance_rate),
                          sd = sd(disturbance_rate)),
              aes(x = year, ymin = mean - sd, ymax = mean + sd),
              alpha = 0.25) +
  geom_line(data = trendline_country_summary %>%
              ungroup() %>%
              mutate(country = country_name_changer(str_to_title(country))),
             aes(x = year, y = mean)) +
  geom_ribbon(data = simulations_stabilize %>%
                ungroup() %>%
                mutate(country = country_name_changer(str_to_title(country))) %>%
                group_by(country, year) %>%
                summarize(mean = mean(disturbance_rate),
                          lower = ifelse(mean - sd(disturbance_rate) < 0, 0, mean - sd(disturbance_rate)),
                          upper = mean + sd(disturbance_rate)),
            aes(x = year, ymin = lower, ymax = upper, 
                fill = "Stabilization at recent mortality rates"),
            alpha = 0.25) +
  geom_line(data = simulations_stabilize %>%
              ungroup() %>%
              mutate(country = country_name_changer(str_to_title(country))) %>%
              group_by(year, country) %>%
              summarize(mean = mean(disturbance_rate),
                        sd = sd(disturbance_rate)),
            aes(x = year, y = mean, col = "Stabilization at recent mortality rates"),
            linetype = "dotted", show.legend = FALSE) +
  geom_ribbon(data = simulations %>%
                ungroup() %>%
                mutate(country = country_name_changer(str_to_title(country))) %>%
                group_by(country, year) %>%
                summarize(mean = mean(disturbance_rate),
                          lower = ifelse(mean - sd(disturbance_rate) < 0, 0, mean - sd(disturbance_rate)),
                          upper = mean + sd(disturbance_rate)),
              aes(x = year, ymin = lower, ymax = upper, fill = "Further increasing canopy mortality rates"),
              alpha = 0.25) +
  geom_line(data = simulations %>%
              ungroup() %>%
              mutate(country = country_name_changer(str_to_title(country))) %>%
              group_by(year, country) %>%
              summarize(mean = mean(disturbance_rate),
                        sd = sd(disturbance_rate)),
            aes(x = year, y = mean, col = "Further increasing canopy mortality rates"),
            linetype = "dotted", show.legend = FALSE) +
  geom_ribbon(data = annual_rates %>%
                ungroup() %>%
                mutate(country = country_name_changer(str_to_title(country))) %>%
                group_by(country) %>%
                mutate(year = sample(2019:2050, length(year), replace = TRUE)) %>%
                group_by(country, year) %>%
                summarize(mean = mean(disturbance_rate),
                          lower = ifelse(mean - sd(disturbance_rate) < 0, 0, mean - sd(disturbance_rate)),
                          upper = mean + sd(disturbance_rate)),
              aes(x = year, ymin = lower, ymax = upper, fill = "Stabilization at past mortality rates"),
              alpha = 0.25) +
  geom_line(data = annual_rates %>%
              ungroup() %>%
              mutate(country = country_name_changer(str_to_title(country))) %>%
              group_by(country) %>%
              mutate(year = sample(2019:2050, length(year), replace = TRUE)) %>%
              group_by(country, year) %>%
              summarize(mean = mean(disturbance_rate),
                        sd = sd(disturbance_rate)),
            aes(x = year, y = mean, col = "Stabilization at past mortality rates"),
            linetype = "dotted", show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Year", 
       y = bquote("Canopy mortality rate (% "*yr^-1*")"),
       linetype = NULL, fill = NULL) +
  scale_x_continuous(breaks = c(1985, 2015, 2045)) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
        strip.text = element_text(size = 6),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA)) +
  facet_wrap(~ country, scales = "free_y") +
  scale_fill_manual(values = colors_scenarios) +
  scale_color_manual(values = colors_scenarios)

ggsave("results/disturbance_trend_country_simulation.pdf", p, width = 7.5, height = 8.5)

# Correlation with growing stock and temperature -----------------------------

#trend_country_summary <- read_csv("results/trend_country_summary.csv")

### Growing stock

growingstock <- readxl::read_excel("data/eea/eea_growing_stock_2015.xlsx") %>%
  mutate(country = tolower(country)) %>%
  mutate(country = case_when(
    country == "bosnia and herzegovina" ~ "bosniaherzegovina",
    country == "czech republic" ~ "czechia",
    country == "the former yugoslav republic of macedonia" ~ "macedonia",
    country == "united kingdom" ~ "unitedkingdom",
    TRUE ~ country
  ))

growingstock <- growingstock %>%
  left_join(trend_country_summary, by = "country") %>%
  left_join(country_grouping, by = c("country" = "country_name_short"))

r2 <- summary(lm(mean ~ growing_stock, data = growingstock))$r.square

p1 <- ggplot(growingstock, aes(x = growing_stock, y = mean)) + 
  geom_point(aes(col = euro_region)) +
  geom_smooth(method = "lm", col = "black", se = TRUE, alpha = 0.5) +
  theme_gray() +
  ggrepel::geom_text_repel(aes(label = country, col = euro_region),
                           box.padding = unit(0.125, "cm"),
                           size = 1.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Paired") +
  labs(x = bquote("Growing stock ("*m^3*ha^-1*")"), 
       y = bquote("Change in Canopy mortality rate (% "*yr^-1*")"),
       col = NULL) +
  theme(legend.position = "none",
        # legend.position = c(1, 0),
        # legend.justification = c(1, 0),
        # legend.text = element_text(size = 5),
        # legend.key.size = unit(0.05, "cm"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank()
        ) +
  annotate("text", x = -Inf, y = Inf, label = paste("R^2:", round(r2, 2)), vjust = 2, hjust = -0.25, parse = TRUE)

### Temperature/precipitation

# Get CRU data

library(data.table)

country_information <- read_csv("data/countries_master.csv")

cntr_cru <- data_frame(country_cru = unique(country_information$country_name),
                       country = country_information$country_name_short) %>%
  mutate(country_cru = case_when(
    country_cru == "Bosnia and Herzegovina" ~ "Bosnia-Herzegovinia",
    country_cru == "The Former Yugoslav Republic of Macedonia" ~ "Macedonia",
    country_cru == "Czech Republic" ~ "Czech_Republic",
    country_cru == "United Kingdom" ~ "United_Kingdom",
    TRUE ~ country_cru
  )) %>%
  filter(!(country_cru %in% c("Russian Federation")))

# Temperature

temp <- vector("list", nrow(cntr_cru))

k <- 0

for (cntr in cntr_cru$country_cru) {
  k <- k + 1
  temp[[k]] <- fread(paste0("https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.03/crucy.1905151143.v4.03/new_countries/tmn/crucy.v4.03.1901.2018.", cntr, ".tmn.per"))
}

temp <- temp %>%
  set_names(cntr_cru$country) %>%
  bind_rows(.id = "country") %>%
  as_tibble() %>%
  dplyr::rename(year = YEAR)

# Precipitation

prec <- vector("list", nrow(cntr_cru))

k <- 0

for (cntr in cntr_cru$country_cru) {
  k <- k + 1
  prec[[k]] <- fread(paste0("https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.03/crucy.1905151143.v4.03/new_countries/pre/crucy.v4.03.1901.2018.", cntr, ".pre.per"))
}

prec <- prec %>%
  set_names(cntr_cru$country) %>%
  bind_rows(.id = "country") %>%
  as_tibble() %>%
  dplyr::rename(year = YEAR)

# Plots

temp_trends <- temp %>%
  #filter(year %in% 1980:2018) %>%
  split(.$country) %>%
  map(., ~ lm(JJA ~ year, data = .)) %>%
  map(., ~ data.frame(temp_trend = coef(.)[2])) %>%
  bind_rows(.id = "country")

prec_trends <- prec %>%
  #filter(year %in% 1980:2018) %>%
  split(.$country) %>%
  map(., ~ lm(JJA ~ year, data = .)) %>%
  map(., ~ data.frame(prec_trend = coef(.)[2])) %>%
  bind_rows(.id = "country")

mod_dat <- growingstock  %>%
  left_join(temp_trends, by = "country") %>%
  left_join(prec_trends, by = "country")

r2 <- summary(lm(mean ~ temp_trend, data = mod_dat))$r.square

p2 <- ggplot(mod_dat %>% mutate(euro_region = as.factor(euro_region)), aes(x = temp_trend, y = mean)) + 
  geom_point(aes(col = euro_region)) +
  geom_smooth(method = "lm", col = "black", se = TRUE, alpha = 0.5) +
  theme_gray() +
  ggrepel::geom_text_repel(aes(label = country, col = euro_region),
                           box.padding = unit(0.125, "cm"),
                           size = 1.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Paired") +
  labs(x = bquote("Change in mean temperature (°C "*yr^-1*")"), 
       y = bquote("Change in Canopy mortality rate (% "*yr^-1*")"),
       col = NULL) +
  theme(legend.position = "none",
        # legend.position = c(1, 0),
        # legend.justification = c(1, 0),
        # legend.text = element_text(size = 5),
        # legend.key.size = unit(0.05, "cm"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank()
  ) +
  annotate("text", x = -Inf, y = Inf, label = paste0("R^2:", round(r2, 2)), vjust = 2, hjust = -0.25, parse = TRUE)

r2 <- summary(lm(mean ~ prec_trend, data = mod_dat))$r.square

p3 <- ggplot(mod_dat %>% mutate(euro_region = as.factor(euro_region)), aes(x = prec_trend, y = mean)) + 
  geom_point(aes(col = euro_region)) +
  geom_smooth(method = "lm", col = "black", se = TRUE, alpha = 0.5) +
  theme_gray() +
  ggrepel::geom_text_repel(aes(label = country, col = euro_region),
                           box.padding = unit(0.125, "cm"),
                           size = 1.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Paired") +
  labs(x = bquote("Change in total precipitation (mm "*yr^-1*")"), 
       y = bquote("Change in Canopy mortality rate (% "*yr^-1*")"),
       col = NULL) +
  theme(legend.position = "none",
        # legend.position = c(1, 0),
        # legend.justification = c(1, 0),
        # legend.text = element_text(size = 5),
        # legend.key.size = unit(0.05, "cm"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank()
  ) +
  annotate("text", x = Inf, y = Inf, label = paste0("R^2:", round(r2, 2)), vjust = 2, hjust = 1.25, parse = TRUE)

p <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave("results/growingstock_temperature_precipitation.pdf", p, width = 7.5, height = 2.5)

ggsave("results/growingstock.pdf", p1, width = 3, height = 3)

# Test for interaction

fit_initial <- lm(median ~ growing_stock * temp_trend * prec_trend, data = mod_dat)

options(na.action = "na.fail")

fits <- MuMIn::dredge(fit_initial)

fits

summary(MuMIn::get.models(fits, subset = delta < 2)[[2]])

# Correlate with FAOSTATS -------------------------------------------------

#annual_rates_country_summary <- read_csv("results/annual_rates_country_summary.csv")
#annual_rates_region_summary <- read_csv("results/annual_rates_region_summary.csv")

faostat <- read_csv("data/faostat/FAOSTAT_data_7-8-2019.csv")

faostat <- faostat %>%
  mutate(country = tolower(gsub(" ", "", gsub(" and ", " ", as.character(Area))))) %>%
  mutate(country = case_when(
    country == "northmacedonia" ~ "macedonia",
    country == "republicofmoldova" ~ "moldova",
    TRUE ~ country)) %>%
  mutate(year = Year) %>%
  mutate(wood_production = Value) %>%
  filter(Item == "Roundwood") %>%
  mutate(wood_production = ifelse(wood_production == 40000 & country == "bosniaherzegovina", NA, wood_production))

### Country

faostat_timesync <- annual_rates_country_summary %>% 
  left_join(faostat, by = c("country", "year")) %>%
  split(.$country) %>%
  map(~ arrange(., year)) %>%
  map(~ mutate(., mean_window = zoo::rollmax(mean, k = 3, fill = NA))) %>%
  bind_rows()

r2 <- faostat_timesync %>%
  split(.$country) %>%
  map(~ lm(mean ~ wood_production, data = na.omit(.))) %>%
  map(~ data.frame(r2 = summary(.)$r.square)) %>%
  bind_rows(.id = "country")

r2_window <- faostat_timesync %>%
  split(.$country) %>%
  map(~ lm(mean_window ~ wood_production, data = na.omit(.))) %>%
  map(~ data.frame(r2 = summary(.)$r.square)) %>%
  bind_rows(.id = "country")

p <- ggplot(faostat_timesync) +
  geom_point(aes(x = mean, y = wood_production / 1000000), shape = 1, size = 0.5) +
  geom_point(aes(x = mean_window, y = wood_production / 1000000), shape = 2, size = 0.5) +
  facet_wrap(~country, scales = "free", ncol = 7) +
  geom_smooth(aes(x = mean, y = wood_production / 1000000), method = "lm", se = TRUE, col = RColorBrewer::brewer.pal(3, "Set1")[1], size = 0.5) +
  geom_smooth(aes(x =  mean_window, y = wood_production / 1000000), method = "lm", se = TRUE, col = RColorBrewer::brewer.pal(3, "Set1")[2], size = 0.5, linetype = "dashed") +
  theme_gray() +
  labs(x = bquote("Canopy mortality rate (% "*yr^-1*")"), 
       y = bquote("Annual roundwood production ("*10^6*m^3*")")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        #strip.text = element_blank()
        strip.text = element_text(size = 6)) +
  #geom_text(aes(x = -Inf, y = Inf, label = str_to_title(country)), vjust = 2, hjust = -0.25, size = 1.5, check_overlap = TRUE) +
  geom_text(data = r2, aes(x = Inf, y = -Inf, label = paste("R^2:", format(round(r2, 2), nsmall = 2))), 
            vjust = -1, hjust = 1.15, parse = TRUE, size = 1.5,
            col = RColorBrewer::brewer.pal(3, "Set1")[1]) +
  geom_text(data = r2_window, aes(x = Inf, y = -Inf, label = paste("R^2:", format(round(r2, 2), nsmall = 2))), 
            vjust = -2.5, hjust = 1.15, parse = TRUE, size = 1.5,
            col = RColorBrewer::brewer.pal(3, "Set1")[2])

ggsave("results/faostat_correlation.pdf", p, width = 7.5, height = 6)

### Europe

faostat_timesync_region <- annual_rates_region_summary %>% 
  left_join(faostat %>% 
              left_join(country_grouping, by = c("country" = "country_name_short")) %>%
              group_by(region = euro_region, year) %>% 
              summarise(wood_production = sum(wood_production, na.rm = TRUE)), by = c("year", "region")) %>%
  arrange(., year) %>%
  mutate(., mean_window = zoo::rollmax(mean, k = 3, fill = NA)) %>%
  ungroup() %>%
  mutate(region = str_to_title(region))

r2 <- faostat_timesync_region %>%
  split(.$region) %>%
  map(~ lm(mean ~ wood_production, data = na.omit(.))) %>%
  map(~ data.frame(r2 = summary(.)$r.square)) %>%
  bind_rows(.id = "region")

r2_window <- faostat_timesync_region %>%
  split(.$region) %>%
  map(~ lm(mean_window ~ wood_production, data = na.omit(.))) %>%
  map(~ data.frame(r2 = summary(.)$r.square)) %>%
  bind_rows(.id = "region")

p <- ggplot(faostat_timesync_region) +
  geom_point(aes(x = mean, y = wood_production / 1000000), shape = 1, size = 0.5) +
  geom_point(aes(x = mean_window, y = wood_production / 1000000), shape = 2, size = 0.5) +
  facet_wrap(~country, scales = "free", ncol = 7) +
  geom_smooth(aes(x = mean, y = wood_production / 1000000), method = "lm", se = TRUE, col = RColorBrewer::brewer.pal(3, "Set1")[1], size = 0.5) +
  geom_smooth(aes(x =  mean_window, y = wood_production / 1000000), method = "lm", se = TRUE, col = RColorBrewer::brewer.pal(3, "Set1")[2], size = 0.5, linetype = "dashed") +
  facet_wrap(~region, scales = "free") +
  theme_gray() +
  labs(x = bquote("Canopy mortality rate (% "*yr^-1*")"), 
       y = bquote("Annual roundwood production ("*10^6*m^3*")")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA)) +
  geom_text(data = r2, aes(x = Inf, y = -Inf, label = paste("R^2:", format(round(r2, 2), nsmall = 2))), 
            vjust = -1, hjust = 1.15, parse = TRUE, size = 2.5, check_overlap = TRUE,
            col = RColorBrewer::brewer.pal(3, "Set1")[1]) +
  geom_text(data = r2_window, aes(x = Inf, y = -Inf, label = paste("R^2:", format(round(r2, 2), nsmall = 2))), 
            vjust = -2.5, hjust = 1.15, parse = TRUE, size = 2.5, check_overlap = TRUE,
            col = RColorBrewer::brewer.pal(3, "Set1")[2])

ggsave("results/faostat_correlation_region.pdf", p, width = 7.5, height = 5.5)


# Correlate with fire -----------------------------------------------------

fire <- read_csv("data/eea/eea_burnt_forest.csv")
names <- strsplit(names(fire), ":") %>% map(~ .[1]) %>% unlist(.)
names[8] <- "Total_EEA39"
names(fire) <- names

p1 <- ggplot(fire %>% filter(Graph == "Burnt Area"), aes(x = Year, y = Total_EEA39)) +
  geom_line() +
  theme_gray() +
  labs(x = "Year", y = "Burnt area (1,000 ha)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))

p2 <- ggplot(fire %>% filter(Graph == "Number of fires"), aes(x = Year, y = Total_EEA39)) +
  geom_line() +
  theme_gray() +
  labs(x = "Year", y = "Number of fires") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))

p <- p1 + p2 + plot_layout(ncol = 2)

ggsave("results/number_fires_bunt_area.pdf", p, width = 5, height = 2.5)

fire <- fire %>%
  filter(Graph == "Burnt Area") %>%
  dplyr::select(-Total, -Total_EEA39, -Percentage, -Graph) %>%
  gather(., key = country, value = area_burnt, -Year) %>%
  dplyr::rename(year = Year) %>%
  mutate(country = tolower(country)) %>%
  left_join(annual_rates_country_summary, by = c("year", "country")) %>%
  mutate(country = str_to_title(country))

p <- ggplot(fire, aes(x = area_burnt, y = mean)) +
  geom_point() +
  facet_wrap(~country, scales = "free", ncol = 3) +
  geom_smooth(method = "lm", col = RColorBrewer::brewer.pal(3, "Set1")[1]) +
  theme_gray() +
  labs(x = "Burnt area (1,000 ha)", y = "Canopy mortality rate") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))

