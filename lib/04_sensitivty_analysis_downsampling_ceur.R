
library(tidyverse)
library(rstanarm)
source("lib/misc/disturbance_summary.R")
options(mc.cores = 1)
library(patchwork)
library(sf)

getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
ceur <- c("austria", "czechia", "germany", "poland", "slovakia", "switzerland")

trend_downsampled <- vector("list", length = 100)

for (j in 1:100) {
  
  ### Get data and downsample
  
  country_information <- read_csv("data/countries_master.csv")
  country_grouping <- read_delim("data/strata/country_grouping.csv", delim = ";")
  country_information$end_year <- 2018
  country_information[country_information$country_name_short %in% ceur, "end_year"] <- 2016
  
  dat_summary_all <- paste0("data/timesync/databases_export/timesync_", ceur, "_export.csv") %>%
    map(read_csv) %>%
    map(~ filter(., plotid %in% sample(unique(.$plotid), 500))) %>% # Downsampling
    map2(.y = arrange(country_information, country_name_short)[country_information$country_name_short %in% ceur, "end_year"][[1]], 
         ~ disturbance_summary2(.x, 
                                grouping = "year",
                                grouping_range = list(year = 1985:as.integer(.y)))) %>%
    set_names(ceur) %>%
    bind_rows(.id = "country")
  
  # Repalce first entry of austria to exclude long-term declines, which were still recorded in the first assessment
  dat_summary_all[dat_summary_all$country == "austria" & dat_summary_all$year == 1985, "disturbance"] <- NA
  
  ### Fit individual trend models
  
  trend <- vector("list", length(ceur))
  
  for (k in 1:length(ceur)) {
    
    cntr <- ceur[k]
    
    dat_summary_subset <- dat_summary_all %>%
      filter(country == cntr) %>%
      mutate(stable = forest - disturbance)
    
    fit_bin <- stan_glmer(cbind(disturbance, stable) ~ year + (1 | year),
                          data = dat_summary_subset,
                          family = binomial(link = "logit"),
                          chains = 4, 
                          iter = 2000)
    
    ## Extract posterior of trend parameter
    
    posterior <- as.matrix(fit_bin)
    
    trend[[k]] <- posterior %>%
      as.data.frame() %>%
      gather() %>%
      filter(key == "year") %>%
      mutate(iter = paste0("iter_", 1:4000)) %>%
      mutate(trend_downsampled = (exp(value) - 1) * 100) %>%
      dplyr::select(-key, -value)
    
  }
  
  trend_downsampled[[j]] <- trend %>%
    set_names(ceur) %>%
    bind_rows(.id = "country")
  
}

trend_downsampled_full <- trend_downsampled %>%
  set_names(1:100) %>%
  bind_rows(.id = "bootstrap")

save(trend_downsampled_full, file = "temp/trend_downsampled_full.RData")
load(file = "temp/trend_downsampled_full.RData")

load(file = "results/trend.RData")

trend <- trend %>%
  filter(country %in% ceur) %>%
  right_join(trend_downsampled_full, by = c("country", "iter"))

trend_summary <- trend %>%
  group_by(country, bootstrap) %>%
  summarise_at(.vars = vars(trend, trend_downsampled), mean)

p <- ggplot(trend_summary) +
  geom_point(aes(x = factor(country, 
                            labels = c("Austria", "Czechia", "Germany", "Poland", "Slovakia", "Switzerland")), 
                 y = trend_downsampled,
                 col = "Sub-sampled"), 
             position = position_jitter(width = 0.1), 
             alpha = 0.25, shape = 19, stroke = NA) +
  geom_point(aes(x = factor(country, 
                            labels = c("Austria", "Czechia", "Germany", "Poland", "Slovakia", "Switzerland")), 
                 y = trend, col = "Full sample"), 
             shape = 3, size = 2) +
  theme_minimal() +
  labs(x = NULL, 
       y = bquote("Canopy mortality rate (% "*yr^-1*")")) +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(size = 6)) +
  scale_color_manual(values = c("red", "grey"))

ggsave("results/down_sampling_ceur.pdf", p, width = 4.5, height = 3.5)


