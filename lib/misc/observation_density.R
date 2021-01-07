
library(tidyverse)
library(sf)

# Load data ---------------------------------------------------------------

# Countries

countries <- read_sf("rawdata/admin/europe/countries_europe.shp")

# Country list

country_master <- read_csv("data/countries_master.csv")

# Observations

observations <- list.files("data/timesync/spectra/", pattern = ".csv$", full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows() %>%
  left_join(country_master, by = "project_id")


# Country statistics ------------------------------------------------------

# Average number of observations per country by year

n_obs_country_year <- observations %>%
  group_by(country_name, plotid, image_year) %>%
  summarize(n = n()) %>%
  group_by(country_name, image_year) %>%
  summarize(n = mean(n)) %>%
  ungroup()

ggplot(n_obs_country_year) +
  geom_line(aes(x = image_year, y = n, col = country_name)) +
  facet_wrap(~country_name) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
        legend.position = "none") +
  labs(x = "Year", y = "Average number of observations") +
  geom_vline(xintercept = c(1999, 2013), linetype = "dashed", alpha = 0.75)

# Average number of observations per country per year

n_obs_country <- n_obs_country_year %>%
  group_by(country_name) %>%
  summarize(n = mean(n)) %>%
  ungroup()

n_obs_country_sf <- countries %>%
  left_join(n_obs_country, by = c("COUNTRY" = "country_name"))

pdf("temp/observation_density.pdf", width = 7, height = 7)
plot(n_obs_country_sf["n"], main = "", key.pos = 1, pal = RColorBrewer::brewer.pal(4, "PuBu"),
     breaks = seq(15, 35, 5))
dev.off()

n_obs_sensor <- observations %>%
  group_by(country_name, plotid, image_year, sensor) %>%
  summarize(n = n()) %>%
  group_by(sensor) %>%
  summarize(n = sum(n)) %>%
  ungroup() %>%
  mutate(prop = (n / sum(n)) * 100)

p <- ggplot(n_obs_sensor %>% mutate(sensor = factor(sensor, levels = c("LT04", "LT05", "LE07", "LC08"))), 
       aes(x = sensor, y = prop)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = NULL, y = "Percent of observations")

ggsave("temp/sensors.pdf", p, width = 2.5, height = 2)

# Total number of images

observations <- observations %>%
  mutate(image_id = paste(sensor, image_year, image_julday, sep = "_"))


