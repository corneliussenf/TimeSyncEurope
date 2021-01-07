disturbance_summary2 <- function(dat,
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
                                            "Wind"),
                                 grouping,
                                 grouping_range) {
  
  dat <- dplyr::as_tibble(dat)
  
  dat_processed <- dplyr::mutate(dat,
                                 agent = dplyr::lead(change_process),
                                 post_disturbance_lc = dplyr::lead(landcover))
  
  dat_processed <- dplyr::filter(dat_processed, landuse == "Forest" & agent %in% agents)
  
  dat_processed$agent <- as.factor(as.character(dat_processed$agent))
  
  dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, .dots = grouping), 
                                    disturbance = length(agent))
  
  dat_processed <- ungroup(dat_processed)
  
  dat_processed <- dplyr::mutate(dat_processed, year = year + 1)
  
  reference_grouping <- expand.grid(grouping_range)
  names(reference_grouping) <- grouping
  
  dat_processed_full <- reference_grouping %>% 
    left_join(dat_processed, by = grouping) %>%
    mutate(disturbance = ifelse(is.na(disturbance), 0, disturbance))
  
  if (length(grouping) == 1) {
    forest_plots <- sum(dplyr::summarize(dplyr::group_by(dat, plotid),
                                         forest = sum(landuse == "Forest", na.rm = TRUE) > 0)$forest)
    dat_processed_full <- dplyr::mutate(dat_processed_full, forest = forest_plots)
  } else {
    forest_plots <- dat %>%
      dplyr::group_by_(.dots = c("plotid", grouping[2])) %>%
      dplyr::summarize(forest = sum(landuse == "Forest", na.rm = TRUE) > 0) %>%
      group_by_(.dots = grouping[2]) %>%
      summarize(forest = sum(forest)) %>%
      na.omit()
      
    dat_processed_full <- dat_processed_full %>%
      left_join(forest_plots, by = grouping[2])
  }
  
  return(dat_processed_full)
}