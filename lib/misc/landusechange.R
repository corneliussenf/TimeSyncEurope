
agents = c("Unknown canopy disturbance", 
           "Biotic", 
           "Fire",
           "Harvest", 
           "Uprooting and breakage", 
           "Gravitational event",
           "Debris", 
           "Decline", 
           "Harvest", 
           "Hydrology", 
           "Other", 
           "Wind")

dat <- list.files("data/timesync/databases_export/", "export.csv", full.names = TRUE) %>%
  map(read_csv) %>%
  map(., ~ mutate(.,
                  agent = lead(change_process),
                  post_disturbance_lc = lead(landcover),
                  post_disturbance_lu = lead(landuse))) %>%
  bind_rows(.)

forestchanges <- dat %>% filter(agent %in% agents)
  
landusechanges <- dat %>% filter(agent %in% agents & post_disturbance_lu != "Forest")

length(unique(paste(landusechanges$project_id, landusechanges$plotid), sep = "-")) / length(unique(paste(forestchanges$project_id, forestchanges$plotid), sep = "-")) * 100
