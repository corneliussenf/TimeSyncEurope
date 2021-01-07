
library(tidyverse)
library(raster)
library(landscapemetrics)
library(gdalUtils)
library(rgdal)
library(sf)
library(patchwork)
library(NLMR)

getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
colors_scenarios <- c("#BB5566", "#004488", "#DDAA33")

country_name_changer <- function(x) {
  x <- str_to_title(x)
  x <- ifelse(x == "Bosniaherzegovina", "Bosnia and Herzegovina", x)
  x <- ifelse(x == "Macedonia", "North Macedonia", x)
  x <- ifelse(x == "Unitedkingdom", "United Kingdom", x)
  return(x)
}

# Data --------------------------------------------------------------------

### Canopy tunrover rates

load(file = "results/annual_rates.RData")
load(file = "results/simulations.RData")

### Level off future simulations at values of 2019 (stabilizing scenario)

simulations_stabilize <- simulations %>%
  filter(year %in% 2019:2020)

### Euro-regions

country_grouping <- read_delim("data/strata/country_grouping.csv", delim = ";")

### Forest proportions

weights <- readxl::read_excel("rawdata/forest_area.xlsx") %>%
  mutate(weight_country = forest_area_km2 / sum(forest_area_km2),
         forest_proportion = forest_area_km2 / country_area_km2) %>%
  dplyr::select(country = country_name_short, forest_proportion, weight_country, forest_area_km2) %>%
  left_join(country_grouping, by = c("country" = "country_name_short")) %>%
  group_by(euro_region) %>%
  mutate(weight_region = forest_area_km2 / sum(forest_area_km2)) %>%
  ungroup() %>%
  dplyr::select(country, forest_proportion, weight_country, weight_region)
  
weights_age <- readxl::read_excel("rawdata/forest_area.xlsx") %>%
  dplyr::select(country = country_name_short, forest_area_km2) %>%
  left_join(country_grouping, by = c("country" = "country_name_short")) %>%
  filter(!(country %in% c("bosniaherzegovina", "greece", "serbia", "montenegro", "macedonia", "moldova"))) %>%
  mutate(weight_country = forest_area_km2 / sum(forest_area_km2)) %>%
  group_by(euro_region) %>%
  mutate(weight_region = forest_area_km2 / sum(forest_area_km2)) %>%
  ungroup() %>%
  dplyr::select(country, euro_region, weight_country, weight_region)
  

# Mortality models --------------------------------------------------------

min_age_death <- c(20, 30, 40, 50)
max_age_death <- c(60, 80, 120, 180)

mortality_props <- vector("list", 4)

for (g in 1:4) {
  
  min_age_death_tmp <- min_age_death[g]
  max_age_death_tmp <- max_age_death[g]
  
  mortality_props[[g]] <- data.frame(age = 1:250) %>%
    mutate(mortality_prop = case_when(
      age < min_age_death_tmp ~ 0, 
      age >= min_age_death_tmp & age < max_age_death_tmp ~ (age - min_age_death_tmp) / (max_age_death_tmp - min_age_death_tmp),
      age >= max_age_death_tmp ~ 1)) %>%
    mutate(mortality_prop = mortality_prop / sum(mortality_prop))
  
}

mortality_props <- mortality_props %>%
  set_names(paste(min_age_death, max_age_death, sep = "-")) %>%
  bind_rows(.id = "death_function")

p <- ggplot(mortality_props, aes(x = age, y = mortality_prop)) +
  geom_line(aes(col = death_function)) +
  labs(x = "Age", y = "P(canopy mortality)", fill = NULL, col = NULL) +
  theme_minimal() +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave("results/mortality_prop_function.pdf", p, width = 3.5, height = 2.25)

# Functions ---------------------------------------------------------------

weighted_Random_Sample <- function (data, weights, n) {
  return(data[order(runif(length(data)) ^ (1 / weights), decreasing = TRUE)][1:n])
}

update_nlm <- function(landscape, mortality_prop, dist_rate) {
  
  landscape <- landscape + 1
  
  landscape_df <- data.frame(age = landscape) %>%
    left_join(mortality_prop, by = "age")
  
  deaths <- weighted_Random_Sample(data = 1:nrow(landscape_df), 
                                   weights = landscape_df$mortality_prop, 
                                   n = round(nrow(landscape_df) * dist_rate, 0))
  
  deaths <- na.omit(deaths)
  
  landscape_df$death <- 1
  landscape_df[deaths, "death"] <- 0
  landscape_df$age <- landscape_df$age * landscape_df$death
  
  landscape <- landscape_df$age
  
  return(landscape)
  
}

distance_mature <- function(x) {
  mature_points <- rasterToPoints(x, fun = function(x) x > 30, sp = TRUE)
  regeneration_points <- rasterToPoints(x, fun = function(x) x < 5, sp = TRUE)
  dist_mature <- pointDistance(mature_points, regeneration_points, lonlat = FALSE)
  return(dist_mature)
}

distance_old <- function(x) {
  old_points <- rasterToPoints(x, fun = function(x) x > 120, sp = TRUE)
  dist_old <- pointDistance(old_points, lonlat = FALSE)
  return(dist_old)
}

create_neutral_landscape <- function(cntr, type, grain, extent, weights) {
  
  print(paste0("Creating neutral landscape for ", cntr))
  
  ### Age distribution
  
  age_distribution <- read.table("data/forest_data_2010-2015.csv", 
                                 sep = ";",
                                 dec = ",",
                                 header = TRUE, 
                                 stringsAsFactors = FALSE) %>%
    mutate(country = tolower(NAME)) %>%
    mutate(country = gsub(" ", "", country)) %>%
    mutate(country = gsub("bosniaandherzegovina", "bosniaherzegovina", country)) %>%
    mutate(country = gsub("theformeryugoslavrepublicofmacedonia", "macedonia", country)) %>%
    mutate(country = gsub("czechrepublic", "czechia", country)) %>%
    filter(country == cntr & Description == "B2 reference" & Step == 2015) %>%
    dplyr::select(Area0_20:Area141_) %>%
    gather() %>%
    mutate(value = value / sum(value)) %>%
    slice(rep(1:n(), each = 20))
  
  if (nrow(age_distribution) > 0 & mean(is.na(age_distribution$value)) < 1) {
    
    age_distribution$age <- 1:160
    
    if (type == "landsat") {
      
      if (!file.exists(paste0("data/landscapes/landscape_init_aggregation", grain, "_", cntr, ".tif"))) {
        
        landscape_outline <- read_sf(paste0("data/countries/", cntr, ".shp"))
        landscape_init <- raster(paste0("data/forestmasks/forestmask_", cntr, ".tif")) %>%
          mask(., landscape_outline)
        
        landscape_init <- aggregate(landscape_init, fact = grain, fun = mean)
        landscape_init <- landscape_init > 0.5
        
        forest_n <- cellStats(landscape_init, sum)
        
        ages_init <- sample(age_distribution$age,
                            size = forest_n,
                            prob = age_distribution$value,
                            replace = TRUE)
        
        landscape_init[landscape_init == 1] <- ages_init
        landscape_init[landscape_init == 0] <- NA
        
        writeRaster(landscape_init, paste0("data/nlm/landscapes/landscape_landsat_aggregation", grain, "_", cntr, ".tif"), overwrite = TRUE)
        
      } else {
        
        landscape_init <- raster(paste0("data/nlm/landscapes/landscape_landsat_aggregation", grain, "_", cntr, ".tif"))
        
      } 
    
    } else if (type == "random") {
        
      forest_prop <- weights[weights$country == cntr, "forest_proportion"][[1]]
      
      landscape_init <- nlm_random(nrow = extent / grain,
                                   ncol = extent / grain,
                                   resolution = grain)
      
      landscape_init <- (landscape_init <= forest_prop)
      
      forest_n <- cellStats(landscape_init, sum)
      
      ages_init <- sample(age_distribution$age,
                          size = forest_n,
                          prob = age_distribution$value,
                          replace = TRUE)
      
      landscape_init[landscape_init == 1] <- ages_init
      landscape_init[landscape_init == 0] <- NA
      
      writeRaster(landscape_init, paste0("data/nlm/landscapes/landscape_random_grain", grain, "_extent", extent, "_", cntr, ".tif"), overwrite = TRUE)
      
    } else if (type == "random_cluster") {
      
      forest_prop <- weights[weights$country == cntr, "forest_proportion"][[1]]
      
      landscape_init <- nlm_randomcluster(nrow = extent / grain,
                                          ncol = extent / grain,
                                          resolution = grain,
                                          p = 0.5,
                                          ai = c(1 - forest_prop, forest_prop))
      
      forest_n <- cellStats(landscape_init, sum)
      
      ages_init <- sample(age_distribution$age,
                          size = forest_n,
                          prob = age_distribution$value,
                          replace = TRUE)
      
      landscape_init[landscape_init == 1] <- ages_init
      landscape_init[landscape_init == 0] <- NA
      
      writeRaster(landscape_init, paste0("data/nlm/landscapes/landscape_randomcluster_grain", grain, "_extent", extent, "_", cntr, ".tif"), overwrite = TRUE)
      
    }
   
    return(landscape_init)
     
  }
  
  return(NULL)

} 

run_neutral_landscape_model <- function(landscape_init, cntr, reps, sim_period, mortality_function) {
  
  mortality_prop <- mortality_props %>% filter(death_function == mortality_function)
  
  mortality_prop[mortality_prop$mortality_prop == 0, "mortality_prop"] <- 0.000000001 # Cannot be exact zero because of fast rejection sampling
  
  ### Loop through repetitions
  
  results_historic <- vector("list", length(reps))
  results_future <- vector("list", length(reps))
  results_future_stabilize <- vector("list", length(reps))
  
  for (j in reps) {
    
    print(paste0("Run ", j, " out of ", max(reps)))
    
    landscape_historic <- landscape_init
    landscape_future <- landscape_init
    landscape_future_stabilize <- landscape_init
    
    predictions_tmp <- annual_rates %>% filter(country == cntr)
    simulations_tmp <- simulations %>% filter(country == cntr)
    simulations_stabilize_tmp <- simulations_stabilize %>% filter(country == cntr)
    
    results_historic_tmp <- vector("list", length(sim_period))
    results_future_tmp <- vector("list", length(sim_period))
    results_future_stabilize_tmp <- vector("list", length(sim_period))
    
    ### Loop through simulation years
    
    for (i in sim_period) {
      
      # Grap canopy mortality rates
      
      dist_rate_historic <- predictions_tmp[predictions_tmp$year %in% sample(c(1985:2018), 1), "disturbance_rate"] / 100
      while(length(dist_rate_historic) < 1) dist_rate_historic <- predictions_tmp[predictions_tmp$year %in% sample(c(1985:2018), 1), "disturbance_rate"] / 100
      dist_rate_historic <- sample(dist_rate_historic, 1)
      #dist_rate_future <- sample(simulations_tmp$disturbance_rate / 100, 1)
      dist_rate_future <- sample(simulations_tmp[simulations_tmp$year == sim_years[i], "disturbance_rate"] / 100, 1)
      dist_rate_future_stabilize <- sample(simulations_stabilize_tmp$disturbance_rate / 100, 1)
      
      # Update landscape
      
      updated_historic <- update_nlm(landscape = as.vector(na.omit(values(landscape_historic))),
                                     mortality_prop = mortality_prop, 
                                     dist_rate = dist_rate_historic)
      
      updated_future <- update_nlm(landscape = as.vector(na.omit(values(landscape_future))), 
                                   mortality_prop = mortality_prop, 
                                   dist_rate = dist_rate_future)
      
      updated_future_stabilize <- update_nlm(landscape = as.vector(na.omit(values(landscape_future_stabilize))), 
                                   mortality_prop = mortality_prop, 
                                   dist_rate = dist_rate_future_stabilize)
      
      landscape_historic[!is.na(landscape_historic)] <- updated_historic
      
      landscape_future[!is.na(landscape_future)] <- updated_future
      
      landscape_future_stabilize[!is.na(landscape_future_stabilize)] <- updated_future_stabilize
      
      # Distances 
      
      dist_mature_historic <- distance_mature(landscape_historic)
      dist_mature_future <- distance_mature(landscape_future)
      dist_mature_future_stabilize <- distance_mature(landscape_future_stabilize)
      
      # Calculate indicators
      
      results_historic_tmp[[i]] <- c(age_median = median(updated_historic, na.rm = TRUE),
                                     age_mean = mean(updated_historic, na.rm = TRUE),
                                     age_iqr = IQR(updated_historic, na.rm = TRUE),
                                     age_sd = sd(updated_historic, na.rm = TRUE),
                                     prop_mature = mean(updated_historic > 30, na.rm = TRUE),
                                     prop_old = mean(updated_historic > 120, na.rm = TRUE),
                                     dist_mature_mean = mean(dist_mature_historic, na.rm = TRUE),
                                     dist_mature_sd = sd(dist_mature_historic, na.rm = TRUE),
                                     dist_mature_median = median(dist_mature_historic, na.rm = TRUE),
                                     dist_mature_iqr = IQR(dist_mature_historic, na.rm = TRUE),
                                     dist_mature_min = min(dist_mature_historic, na.rm = TRUE),
                                     dist_mature_max = max(dist_mature_historic, na.rm = TRUE),
                                     shannon_diversity = vegan::diversity(updated_historic, "shannon"),
                                     simpson_diversity = vegan::diversity(updated_historic, "simpson"))
      
      results_future_tmp[[i]] <- c(age_median = median(updated_future, na.rm = TRUE),
                                   age_mean = mean(updated_future, na.rm = TRUE),
                                   age_iqr = IQR(updated_future, na.rm = TRUE),
                                   age_sd = sd(updated_future, na.rm = TRUE),
                                   prop_mature = mean(updated_future > 30, na.rm = TRUE),
                                   prop_old = mean(updated_future > 120, na.rm = TRUE),
                                   dist_mature_mean = mean(dist_mature_future, na.rm = TRUE),
                                   dist_mature_sd = sd(dist_mature_future, na.rm = TRUE),
                                   dist_mature_median = median(dist_mature_future, na.rm = TRUE),
                                   dist_mature_iqr = IQR(dist_mature_future, na.rm = TRUE),
                                   dist_mature_min = min(dist_mature_future, na.rm = TRUE),
                                   dist_mature_max = max(dist_mature_future, na.rm = TRUE),
                                   shannon_diversity = vegan::diversity(updated_future, "shannon"),
                                   simpson_diversity = vegan::diversity(updated_future, "simpson"))
      
      results_future_stabilize_tmp[[i]] <- c(age_median = median(updated_future_stabilize, na.rm = TRUE),
                                             age_mean = mean(updated_future_stabilize, na.rm = TRUE),
                                             age_iqr = IQR(updated_future_stabilize, na.rm = TRUE),
                                             age_sd = sd(updated_future_stabilize, na.rm = TRUE),
                                             prop_mature = mean(updated_future_stabilize > 30, na.rm = TRUE),
                                             prop_old = mean(updated_future_stabilize > 120, na.rm = TRUE),
                                             dist_mature_mean = mean(dist_mature_future_stabilize, na.rm = TRUE),
                                             dist_mature_sd = sd(dist_mature_future_stabilize, na.rm = TRUE),
                                             dist_mature_median = median(dist_mature_future_stabilize, na.rm = TRUE),
                                             dist_mature_iqr = IQR(dist_mature_future_stabilize, na.rm = TRUE),
                                             dist_mature_min = min(dist_mature_future_stabilize, na.rm = TRUE),
                                             dist_mature_max = max(dist_mature_future_stabilize, na.rm = TRUE),
                                             shannon_diversity = vegan::diversity(updated_future_stabilize, "shannon"),
                                             simpson_diversity = vegan::diversity(updated_future_stabilize, "simpson"))
      
    }
    
    results_historic[[j]] <- results_historic_tmp %>%
      do.call("rbind", .) %>%
      as.data.frame(.) %>%
      mutate(sim_period = sim_period)
    
    results_future[[j]] <- results_future_tmp %>%
      do.call("rbind", .) %>%
      as.data.frame(.) %>%
      mutate(sim_period = sim_period)
    
    results_future_stabilize[[j]] <- results_future_stabilize_tmp %>%
      do.call("rbind", .) %>%
      as.data.frame(.) %>%
      mutate(sim_period = sim_period)
    
  }
  
  results_historic <- results_historic %>%
    set_names(reps) %>%
    bind_rows(.id = "reps")
  
  results_future <- results_future %>%
    set_names(reps) %>%
    bind_rows(.id = "reps")
  
  results_future_stabilize <- results_future_stabilize %>%
    set_names(reps) %>%
    bind_rows(.id = "reps")
  
  return(list(results_historic, results_future, results_future_stabilize))
  
}

run_neutral_landscape_model_plotonly <- function(landscape_init, cntr, sim_period, mortality_function) {
  
  mortality_prop <- mortality_props %>% filter(death_function == mortality_function)
  
  mortality_prop[mortality_prop$mortality_prop == 0, "mortality_prop"] <- 0.000000001 # Cannot be exact zero because of fast rejection sampling
  
  predictions_tmp <- annual_rates %>% filter(country == cntr)
  simulations_tmp <- simulations %>% filter(country == cntr)
  
  ### Loop through simulation years
  
  plots <- vector("list", length(sim_period))
  
  for (i in sim_period) {
    
    # Grap canopy mortality rates
    
    dist_rate_future <- sample(simulations_tmp[simulations_tmp$year == sim_years[i], "disturbance_rate"] / 100, 1)
    
    # Update landscape
    
    updated_future <- update_nlm(landscape = as.vector(na.omit(values(landscape_init))), 
                                 mortality_prop = mortality_prop, 
                                 dist_rate = dist_rate_future)
    
    landscape_init[!is.na(landscape_init)] <- updated_future
    
    # plots[[i]] <- ggplot(data = rasterToPoints(landscape_init) %>% as.data.frame(),
    #                      aes(x = x, y = y, fill = clumps)) +
    #   geom_raster() +
    #   theme_void() +
    #   theme(panel.background = element_rect(fill = "grey", color = "black", size = 2),
    #         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    #         legend.key.height = unit(2.2, "cm"),
    #         legend.key.width = unit(0.5, "cm")) +
    #   scale_fill_viridis_c(direction = -1, limits = c(0, 200)) +
    #   scale_x_continuous(expand = c(0.01, 0.01)) +
    #   scale_y_continuous(expand = c(0.01, 0.01)) +
    #   labs(fill = "Age")
    
    plots[[i]] <- rasterToPoints(landscape_init) %>% 
      as.data.frame() %>%
      mutate(distrate = dist_rate_future)
    
  }
  
  return(plots)
  
}

# Run neutral landscape models --------------------------------------------

### Simulation parameters

sim_period <- 1:32
sim_years <- 2019:(2019 + length(sim_period) - 1)
max_age <- (160 + length(sim_period))
reps <- 1:15

### Landscape parameters

grains <- c(100)
extents <- c(10000)

### Countries

countries <- list.files("data/timesync/databases_export/", "export.csv") %>%
  strsplit(., "_") %>%
  map(~ .[2]) %>%
  unlist()

### Create look-up table

simulation_lut <- expand.grid(country = countries,
                              grain = grains,
                              extent = extents,
                              mortality_function = paste(min_age_death, max_age_death, sep = "-"),
                              landscape = c("random", "random_cluster"), 
                              stringsAsFactors = FALSE) %>%
  filter(!(country %in% c("bosniaherzegovina", "greece", "serbia", "montenegro", "macedonia", "moldova")))

write_csv(simulation_lut, "data/nlm/simulation_lut.csv")

# simulation_lut <- simulation_lut %>% 
#   filter(country == "austria" & mortality_function == "30-80" & landscape == "random_cluster")

### Loop through look-up table

# #---
# simulation_lut <- simulation_lut[1:29, ]
# reps <- 1:3
# #---

out_historic <- vector("list", nrow(simulation_lut))
out_future <- vector("list", nrow(simulation_lut))
out_future_stabilize <- vector("list", nrow(simulation_lut))

for (k in 1:nrow(simulation_lut)) {
  
  print(paste0(k, " out of ", nrow(simulation_lut)))
  
  neutral_landscape <- create_neutral_landscape(cntr = simulation_lut[k, "country"], 
                                                type = simulation_lut[k, "landscape"], 
                                                grain = simulation_lut[k, "grain"], 
                                                extent = simulation_lut[k, "extent"], 
                                                weights = weights)
  
  simulation_runs <- run_neutral_landscape_model(landscape_init = neutral_landscape, 
                                                 cntr = simulation_lut[k, "country"],
                                                 reps = reps, 
                                                 sim_period = sim_period, 
                                                 mortality_function = simulation_lut[k, "mortality_function"])
  
  out_historic[[k]] <- simulation_runs[[1]]
  out_future[[k]] <- simulation_runs[[2]]
  out_future_stabilize[[k]] <- simulation_runs[[3]]
  
  removeTmpFiles(h = 0)
  
}

out_historic <- out_historic %>%
  bind_rows()
save(out_historic, file = "temp/out_historic.RData")
out_historic <- cbind(out_historic, 
                      simulation_lut %>% 
                        slice(rep(1:n(), each = length(reps) * length(sim_period))))

out_future <- out_future %>%
  bind_rows()
save(out_future, file = "temp/out_future.RData")
out_future <- cbind(out_future, 
                    simulation_lut %>% 
                      slice(rep(1:n(), each = length(reps) * length(sim_period))))

out_future_stabilize <- out_future_stabilize %>%
  bind_rows()
save(out_future_stabilize, file = "temp/out_future_stabilize.RData")
out_future_stabilize <- cbind(out_future_stabilize, 
                    simulation_lut %>% 
                      slice(rep(1:n(), each = length(reps) * length(sim_period))))

save(out_historic, file = "data/nlm/neutral_landscape_runs_historic.RData")
save(out_future, file = "data/nlm/neutral_landscape_runs_future.RData")
save(out_future_stabilize, file = "data/nlm/neutral_landscape_runs_future_stabilize.RData")

# Figures -----------------------------------------------------------------

load(file = "data/nlm/neutral_landscape_runs_historic.RData")
load(file = "data/nlm/neutral_landscape_runs_future.RData")
load(file = "data/nlm/neutral_landscape_runs_future_stabilize.RData")

dat <- list(out_historic, out_future, out_future_stabilize) %>%
  set_names(c("Historic", "Future", "Future_stabilize")) %>%
  bind_rows(.id = "period") %>%
  mutate(year = sim_period)

dat <- dat %>%
  left_join(country_grouping, by = c("country" = "country_name_short"))

scenario_names <- c("Current", 
                    "Stabilization at past rates", 
                    "Stabilization at current rates", 
                    "Further increasing rates")

### Plot development trajectories

p <- ggplot(dat %>%
         mutate(scenario = factor(period, 
                                  labels = c("Further increasing mortality rates",
                                             "Stabilization at recent mortality rates",
                                             "Stabilization at past mortality rates"))), 
       aes(x = sim_period, y = age_median, 
                group = interaction(reps, landscape, mortality_function, scenario), 
                col = scenario)) + 
  geom_line(alpha = 0.1) + 
  facet_wrap(~country_name_changer(country), scales = "free") +
  scale_color_manual(values = colors_scenarios) +
  labs(x = "Years of forest development", 
       y = "Median forest age",
       col = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
        strip.text = element_text(size = 6),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(0, "cm")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

ggsave("results/nlm_development_trajectories.pdf", p, width = 7.5, height = 8.5)

### Europe age distribution

# Summary

summary1 <- expand.grid(age_class = paste0(c(0, 30, 60), "-", c(30, 60, 120)),
                        period_plot = scenario_names,
                        landscape = unique(dat$landscape),
                        mortality_function = unique(dat$mortality_function))

summary2 <- list(dat %>% filter(period == "Historic" & year == 1),
                 dat %>% filter(period == "Historic" & year == 32),
                 dat %>% filter(period == "Future_stabilize" & year == 32),
                 dat %>% filter(period == "Future" & year == 32)) %>%
  set_names(scenario_names) %>%
  bind_rows(.id = "period_plot") %>%
  mutate(age_class = cut(age_median, c(0, 30, 60, 120), include.lowest = TRUE, labels = paste0(c(0, 30, 60), "-", c(30, 60, 120)))) %>%
  group_by(age_class, period_plot, landscape, mortality_function) %>%
  summarize(n = n())

summary <- summary1 %>% left_join(summary2) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  group_by(period_plot, landscape, mortality_function) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  mutate(period_plot = factor(period_plot, 
                              levels = scenario_names))

summary <- summary %>%
  group_by(age_class, period_plot) %>%
  summarize(n_mean = mean(n),
            n_sd = sd(n, na.rm = TRUE),
            p_mean = mean(p),
            p_sd = sd(p, na.rm = TRUE))

write_csv(summary, "results/nlm_results_summary.csv")

# Flow plot

library(ggalluvial)

plotdat1 <- expand.grid(age_class = paste0(c(0, 30, 60), "-", c(30, 60, 120)),
                        period_plot = scenario_names)

plotdat2 <- list(dat %>% filter(period == "Historic" & year == 1),
                 dat %>% filter(period == "Historic" & year == 32),
                 dat %>% filter(period == "Future_stabilize" & year == 32),
                 dat %>% filter(period == "Future" & year == 32)) %>%
  set_names(scenario_names) %>%
  bind_rows(.id = "period_plot") %>%
  mutate(age_class = cut(age_median, c(0, 30, 60, 120), include.lowest = TRUE, labels = paste0(c(0, 30, 60), "-", c(30, 60, 120)))) %>%
  group_by(age_class, period_plot, mortality_function, landscape) %>%
  summarize(n = n())

plotdat <- plotdat1 %>% left_join(plotdat2) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  group_by(period_plot, mortality_function, landscape) %>%
  mutate(p = n / sum(n)) %>%
  group_by(age_class, period_plot) %>%
  summarize(n = mean(n),
            p = mean(p)) %>%
  ungroup() %>%
  mutate(period_plot = factor(period_plot, 
                              levels = scenario_names[c(2, 1, 3, 4)],
                              c("Stabilization\nat past\nmortality rates",
                                "Current", 
                                "Stabilization\nat recent\nmortality rates",
                                "Further\nincreasing\nmortality rates"))) %>%
  mutate(age_class = factor(age_class, levels = c("60-120", "30-60", "0-30")))

p <- ggplot(plotdat,
       aes(x = period_plot, 
           alluvium = age_class,
           fill = period_plot,
           stratum = age_class,
             label = ifelse(p <= 0.03, "", ifelse(p < 0.1, paste0(age_class, " yr."), paste0(age_class, " yr.\n(", round(p * 100, 0), "%)"))),
           y = p)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(aes(alpha = age_class), width = 0.6, fill = "white", col = "black") +
  geom_stratum(aes(alpha = age_class), width = 0.6) +
  scale_alpha_manual(values = c(1, 0.75, 0.5)) +
  scale_fill_manual(values = c(colors_scenarios[3],
                               "grey",
                               colors_scenarios[2],
                               colors_scenarios[1])) +
  geom_text(stat = "stratum", size = 2.5) +
  labs(x = NULL, y = "Proportion of countries") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(0, "cm"),
        plot.margin = margin(1, 15, 1, 20)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  annotate("segment", x = 1.9, y = -0.03, xend = 1, yend = -0.03, col = "black", arrow = arrow(length = unit(0.05, "cm"))) +
  annotate("segment", x = 2.1, y = -0.03, xend = 4, yend = -0.03, col = "black", arrow = arrow(length = unit(0.05, "cm")))

ggsave("results/age_distribution_v03.pdf", p, width = 3.5, height = 3.5)

p <- ggplot(summary %>% 
              mutate(period_plot = factor(period_plot, 
                                          levels = scenario_names[c(2, 1, 3, 4)],
                                          c("Stabilization\nat past\nmortality rates",
                                            "Current", 
                                            "Stabilization\nat recent\nmortality rates",
                                            "Further\nincreasing\nmortality rates"))) %>%
              mutate(landscape_label = case_when(
                landscape == "random" ~ "Random",
                landscape == "random_cluster" ~ "Clustered")) %>%
              mutate(age_class = factor(age_class, levels = c("60-120", "30-60", "0-30"))),
            aes(x = period_plot, 
                alluvium = age_class,
                fill = period_plot,
                stratum = age_class,
                label = ifelse(p <= 0.03, "", ifelse(p < 0.1, paste0(age_class, " yr."), paste0(age_class, " yr.\n(", round(p * 100, 0), "%)"))),
                y = p)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(aes(alpha = age_class), width = 0.6, fill = "white", col = "black") +
  scale_alpha_manual(values = c(1, 0.75, 0.5)) +
  scale_fill_manual(values = c(colors_scenarios[3],
                               "grey",
                               colors_scenarios[2],
                               colors_scenarios[1])) +
  geom_stratum(aes(alpha = age_class), width = 0.6) +
  geom_text(stat = "stratum", size = 2.15) +
  labs(x = NULL, y = "Proportion of countries") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(0, "cm")) +
  scale_y_continuous(expand = c(0.02, 0)) +
  facet_grid(mortality_function~landscape_label)

ggsave("results/age_distribution_v03_full.pdf", p, width = 3.5 * 2, height = 3.5 * 4)

### Country differences using indicators

# Density plots

dat_country_difference <- dat %>%
  filter(year == 32) %>%
  dplyr::select(period:landscape, euro_region, -sim_period, -extent, -grain) %>%
  gather(key = metric, value = value, -period, -reps, -country, -mortality_function, -landscape, -euro_region) %>%
  spread(key = period, value = value) %>%
  mutate(difference_increase = Future - Historic,
         difference_increase_rel = difference_increase / Historic * 100,
         difference_stabilize = Future_stabilize - Historic,
         difference_stabilize_rel = difference_stabilize / Historic * 100) %>%
  filter(metric %in% c("prop_mature", "dist_mature_mean", "prop_old", "shannon_diversity"))

dat_country_difference_mean <- dat_country_difference %>%
  group_by(country, metric, euro_region) %>%
  summarize_at(.vars = vars("difference_increase", "difference_increase_rel",
                            "difference_stabilize", "difference_stabilize_rel"), 
               .funs = list(mean = mean, sd = sd)) %>%
  ungroup() %>%
  dplyr::select(country, metric, difference_increase_mean, difference_stabilize_mean) %>%
  gather(key = scenario, value = difference, -country, -metric)

# dat_country_difference <- dat_country_difference %>%
#   left_join(dat_country_difference_mean, by = c("country", "metric"))

p1 <- pmap(list(a = dat_country_difference %>%
                 mutate(difference = difference_stabilize) %>%
                 mutate(metric = factor(metric,
                                        levels = c("prop_mature", "dist_mature_mean", "prop_old", "shannon_diversity"),
                                        labels = c("Proportion of\nmature forests", "Distance to\nmature forest", "Proportion of\nold forests", "Age class diversity"))) %>%
                 split(.$metric),
          b = c(1, -1, 1, 1),
          #c = c("a) Resilience", "", "b) Biodiversity", ""),
          d = list(bquote(Delta~"proportion of stands > 30 yr."),
                   bquote(Delta~"distance to stands > 30 yr. (m)"),
                   bquote(Delta~"proportion of stands > 120 yr."),
                   bquote(Delta~"Shannon-Index H'")),
          e = list(element_text(size = 6), element_text(size = 6), element_text(size = 6), element_text(size = 6)),
          # c = c("Dispersal limitations", 
          #       "Seed sources", 
          #       "Old forests", 
          #       "Age class diversity"),
          c = c("A)", "B)", "C)", "D)"),
          label.changer = c(FALSE, TRUE, FALSE, FALSE),
          xlim1 = c(-1.2, -2100, -0.6, -2.3),
          xlim2 = c(0.5, 900, 0.25, 0.94),
          tit = list("A) Mortality rates stabilize at recent rates", NULL, NULL, NULL)),
       function(a, b, c, d, e, label.changer, xlim1, xlim2, tit, subtit) {
         p <- ggplot(a, aes(x = difference * b, group = country, col = difference * b)) + 
           #stat_density(geom = "line", position = "identity", adjust = 3) +
           geom_vline(aes(xintercept = difference * b, col = difference * b), alpha = 0.25, size = 0.5) +
           #geom_vline(aes(xintercept = difference_mean * b, col = difference_mean * b)) +
           geom_vline(xintercept = 0, col = "black") +
           scale_color_gradient2(low = "#a50026", mid = "grey", high = "#313695") +
           theme_minimal() +
           labs(x = d,
                y = NULL,
                fill = NULL,
                col = NULL,
                title = tit) +
           theme(legend.position = "none",
                 plot.title = element_text(size = 8),
                 axis.title = element_text(size = 7),
                 axis.text.x = element_text(size = 7),
                 axis.text.y = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.border = element_rect(fill = NA),
                 strip.background = element_blank(),
                 strip.text = element_blank(),
                 # strip.background = element_blank(),
                 # strip.text.y = element_text(angle = 0, size = 8),
                 plot.margin = unit(x = c(0.5, 0, 0, 0), units = "mm")) +
           facet_wrap(~metric) +
           ylim(0, 1) +
           xlim(xlim1, xlim2)
         if (label.changer) p <- p + scale_x_continuous(labels = function(x) x * -1, limits = c(xlim1, xlim2))
         return(p)
       }) %>%
  wrap_plots(ncol = 1)

p2 <- pmap(list(a = dat_country_difference %>%
                  mutate(difference = difference_increase) %>%
                  mutate(metric = factor(metric,
                                         levels = c("prop_mature", "dist_mature_mean", "prop_old", "shannon_diversity"),
                                         labels = c("Proportion of\nmature forests", "Distance to\nmature forest", "Proportion of\nold forests", "Age class diversity"))) %>%
                  split(.$metric),
                b = c(1, -1, 1, 1),
                #c = c("a) Resilience", "", "b) Biodiversity", ""),
                d = list(bquote(Delta~"proportion of stands > 30 yr."),
                         bquote(Delta~"distance to stands > 30 yr. (m)"),
                         bquote(Delta~"proportion of stands > 120 yr."),
                         bquote(Delta~"Shannon-Index H'")),
                e = list(element_text(size = 6), element_text(size = 6), element_text(size = 6), element_text(size = 6)),
                # c = c("Dispersal limitations", 
                #       "Seed sources", 
                #       "Old forests", 
                #       "Age class diversity"),
                c = c("A)", "B)", "C)", "D)"),
                label.changer = c(FALSE, TRUE, FALSE, FALSE),
                xlim1 = c(-1.2, -2100, -0.6, -2.3),
                xlim2 = c(0.5, 900, 0.25, 0.94),
                tit = list("B) Mortality rates further increase", NULL, NULL, NULL)),
           function(a, b, c, d, e, label.changer, xlim1, xlim2, tit, subtit) {
             p <- ggplot(a, aes(x = difference * b, group = country, col = difference * b)) + 
               #stat_density(geom = "line", position = "identity", adjust = 3) +
               geom_vline(aes(xintercept = difference * b, col = difference * b), alpha = 0.25, size = 0.5) +
               #geom_vline(aes(xintercept = difference_mean * b, col = difference_mean * b)) +
               geom_vline(xintercept = 0, col = "black") +
               scale_color_gradient2(low = "#a50026", mid = "grey", high = "#313695") +
               theme_minimal() +
               labs(x = d,
                    y = NULL,
                    fill = NULL,
                    col = NULL,
                    title = tit) +
               theme(legend.position = "none",
                     plot.title = element_text(size = 8),
                     axis.title = element_text(size = 7),
                     axis.text.x = element_text(size = 7),
                     axis.text.y = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.border = element_rect(fill = NA),
                     strip.background = element_blank(),
                     strip.text.y = element_text(angle = 0, size = 8),
                     plot.margin = unit(x = c(0.5, 0, 0, 0), units = "mm")) +
               facet_wrap(~metric, strip.position = "right") +
               ylim(0, 1) +
               xlim(xlim1, xlim2)
             if (label.changer) p <- p + scale_x_continuous(labels = function(x) x * -1, limits = c(xlim1, xlim2))
             return(p)
           }) %>%
  wrap_plots(ncol = 1)

ggsave("results/indicators_v02_A.pdf", p1, width = 3, height = 3)
ggsave("results/indicators_v02_B.pdf", p2, width = 4, height = 3)

