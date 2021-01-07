
#### Libraries ####

library(raster)
library(tidyverse)
library(sf)
library(velox)
library(gdalUtils)

#### Parameters ####

samplesize_initial <- 750
samplesize_reduced <- 500

countries_already_sampled <- c("austria", "czechia", "germany", "poland", "slovakia", "switzerland")

#### Settings ####

dir.create("temp", showWarnings = FALSE)

rasterOptions(datatype = "FLT4S", 
              progress = "", 
              #tmpdir = "temp", 
              tmptime = 6)

#### Load and clean data ####

countries <- read_sf("rawdata/admin/europe/countries_europe.shp")

# Rename countries with empty spaces

countries <- countries %>%
  mutate(COUNTRY_SHORT = case_when(
    COUNTRY == "Bosnia and Herzegovina" ~ "bosniaherzegovina",
    COUNTRY == "Russian Federation" ~ "russia",
    COUNTRY == "The Former Yugoslav Republic of Macedonia" ~ "macedonia",
    COUNTRY == "United Kingdom" ~ "unitedkingdom",
    COUNTRY == "Czech Republic" ~ "czechia",
    TRUE ~ tolower(COUNTRY)
  ))

# Add project-id (consecutive number from 1 to n)

countries$PROJECT_ID <- 1:length(countries$COUNTRY_SHORT)

# Export country outlines

countries %>%
  split(.$COUNTRY_SHORT) %>%
  map2(.y = sort(unique(countries$COUNTRY_SHORT)), ~ write_sf(.x, paste0("data/countries/", .y, ".shp")))

# Create master-table, which will be updated in following steps

countries_master <- countries %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  dplyr::rename(country_name = COUNTRY, country_name_short = COUNTRY_SHORT, 
                iso_code = ISO_CC, country_area_km2 = AREA_km2,
                project_id = PROJECT_ID)

#### Loop through countries and draw sample (leave out countries already sampled) ####

dir.create("data/samples", showWarnings = FALSE)
dir.create("data/forestmasks", showWarnings = FALSE)

countries_to_sample <- unique(countries$COUNTRY_SHORT)
countries_to_sample <- countries_to_sample[!countries_to_sample %in% countries_already_sampled]


for (cntr in countries_to_sample) {
  
  print(cntr)
  
  countries_tmp <- countries %>% filter(., COUNTRY_SHORT == cntr)
  countries_tmp_extent <- extent(countries_tmp)
  
  gdalbuildvrt(gdalfile = "rawdata/forestcover/europe_lcchange_19852015e_forestmask.bsq", 
               output.vrt = "temp/tmp.vrt", 
               te = st_bbox(countries_tmp))
  
  forest_tmp <- raster("temp/tmp.vrt")
  
  writeRaster(forest_tmp, paste0("data/forestmasks/forestmask_", tolower(cntr), ".tif"), 
              datatype = "INT1U", overwrite = TRUE)
  
  vx <- velox(forest_tmp)
  countries_tmp$id <- 1
  vx$rasterize(countries_tmp, field = "id", band = 1, background = NA)
  countries_tmp_raster <- vx$as.RasterLayer(band = 1)
  forest_tmp <- forest_tmp * countries_tmp_raster
  
  rm(countries_tmp_raster)
  
  forest_tmp[forest_tmp == 0] <- NA
  
  samples <- sampleRandom(forest_tmp, samplesize_initial, na.rm = TRUE, sp = TRUE)
  
  if (length(samples) < samplesize_initial) {
    n_sample <- length(samples)
    while (n_sample < samplesize_initial) {
      samples <- rbind(samples, sampleRandom(forest_tmp, samplesize_initial, na.rm = TRUE, sp = TRUE))
      n_sample <- n_sample + length(samples)
    }  
  }
  
  samples <- samples[1:samplesize_initial, ]
  
  samples$plotid <- 1:nrow(samples)
  samples$layer <- NULL
  
  rgdal::writeOGR(samples, dsn = "data/samples", layer = paste0("samples_", cntr), driver = "ESRI Shapefile")
  
  rm(samples)
  rm(forest_tmp)
  file.remove("temp/tmp.vrt")
  
}

#### Loop through countries already sampled ####

samples_ceur <- read.csv("data/senfetal2018/disturbances_senfetal2018.csv")

samples_ceur <- samples_ceur %>%
  group_by(country, plotid, xcoord, ycoord) %>%
  summarize(disturbance = sum(disturbance == "Yes")) %>%
  ungroup() %>%
  mutate(country = tolower(country))

for(cntr in countries_already_sampled) {
  
  samples <- samples_ceur %>% filter(country == cntr)
  samples <- st_as_sf(samples, coords = c("xcoord", "ycoord"), crs = as.character(st_crs(countries)[[2]]))
  
  samples$country <- NULL
  #samples <- samples[1:samplesize_initial, ]
  samples$plotid2 <- samples$plotid
  samples$plotid <- 1:nrow(samples)
  
  samples %>%
    as.data.frame() %>%
    filter(disturbance > 0) %>%
    write_csv(., paste0("temp/disturbance_plots_", cntr, "_full.csv"))
  
  samples$disturbance <- NULL
  
  write_sf(samples, paste0("data/samples/samples_", cntr, "_full.shp"))
  
}

#### Pixel-outlines and KML files ####

for (cntr in unique(countries$COUNTRY_SHORT)) {
  
  samples <- shapefile(paste0("data/samples/samples_", tolower(cntr), ".shp"))
  
  samples_pixel <- rgeos::gBuffer(samples, byid = TRUE, id = samples$plotid, width = 15, capStyle = "SQUARE")
  rgdal::writeOGR(samples_pixel, dsn = "data/samples", layer = paste0("samples_pixel_", tolower(cntr)), 
                  driver = "ESRI Shapefile", overwrite_layer = TRUE)
  plotKML::plotKML(samples_pixel, file.name = paste0("data/samples/samples_pixel_", tolower(cntr), ".kml"))

}

#### Export samples to CSV ####

dir.create("data/plots", showWarnings = FALSE)

for (cntr in countries_already_sampled) {
  
  samples <- read_sf(paste0("data/samples/samples_", tolower(cntr), "_full.shp"))
  #samples <- samples %>% filter(plotid %in% 1:samplesize_reduced)
  samples_table <- data.frame(project_id = rep(countries_master$project_id[countries_master$country_name_short == cntr], length(samples$plotid)),
                              plotid = samples$plotid,
                              plotid_old = samples$plotid2,
                              x = sf::st_coordinates(samples)[,1],
                              y = sf::st_coordinates(samples)[,2])
  #write.csv(samples_table, paste0("data/plots/plots_n", samplesize_reduced, "_", tolower(cntr), ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(samples_table, paste0("data/plots/plots_full_", tolower(cntr), ".csv"), row.names = FALSE, quote = FALSE)
  
}

#### Add eco/biogeo/political regions ####

samples_files <- list.files("data/samples/", ".shp", full.names = TRUE)
samples_files <- samples_files[-grep("pixel", samples_files)]
samples <- samples_files %>% 
  map(read_sf) %>% 
  map(., ~ dplyr::select(., plotid)) %>%
  map2(.y = gsub(".shp", "", unlist(map(strsplit(basename(samples_files), "_"), ~.[[2]]))), ~ mutate(., country = .y))
samples <- do.call(rbind, samples)

olson <- read_sf("rawdata/biome/terrestrial_ecoregions_olson.shp")
olson <- st_transform(olson, crs = st_crs(samples))

olson_intersect <- st_intersection(samples, olson)

dir.create("data/strata")

write_csv(olson_intersect, "data/strata/strata_olson.csv")

