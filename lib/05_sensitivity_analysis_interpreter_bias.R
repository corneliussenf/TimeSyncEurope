
library(tidyverse)
library(rstanarm)
source("lib/misc/disturbance_summary_sensitivity.R")
options(mc.cores = 1)
library(patchwork)
library(sf)

### Settings

years <- 1985:2018
plots <- 500
errors <- seq(0.1, 0.8, 0.1)
n_sim <- 30

### Repeat simulation

trends <- vector("list", n_sim)
intercepts <- vector("list", n_sim)

for (sim_num in 1:n_sim) {
  
  #sim_num <- 1
  
  ### Create data
  
  dist_rate_true <- runif(1, 0.01, 0.05)
  dist_rate_true <- ifelse(dist_rate_true < 0, 0, dist_rate_true)
  var_years <- rnorm(1, 0.1, 0.05)
  var_years <- ifelse(var_years < 0, 0, var_years)
  trend_true <- rnorm(1, 0.04, 0.05)
  
  avg <- trend_true * years - 1985 * trend_true + boot::logit(dist_rate_true)
  
  sim_dat <- data.frame(year = years, average = boot::inv.logit(avg))
  
  avg <- avg + rnorm(length(avg), 0, var_years)
  
  sim_dat$disturbance_rate <- boot::inv.logit(avg)
  
  dist_observed <- replicate(plots, sapply(boot::inv.logit(avg), function(x) rbinom(1, 1, x)))
  
  ### Estimate from correct data
  
  disturbance_observed_noerror <- apply(dist_observed, 1, sum)
  dist_observed_noerror <- data.frame(disturbance = disturbance_observed_noerror, stable = forest - disturbance_observed_noerror, year = years)
  
  fit <- stan_glmer(cbind(disturbance, stable) ~ year + (1 | year),
                    data = dist_observed_noerror,
                    family = binomial(link = "logit"),
                    chains = 2, 
                    iter = 2000)
  
  trend_noerror <- as.matrix(fit) %>%
    as.data.frame() %>%
    gather() %>%
    filter(key == "year") %>%
    mutate(iter = paste0("iter_", 1:2000)) %>%
    mutate(trend = (exp(value) - 1) * 100) %>%
    dplyr::select(-key, -value) %>%
    mutate(error = 0)
  
  intercept_noerror <- as.matrix(fit) %>%
    as.data.frame() %>%
    mutate(dist_rate = boot::inv.logit(`(Intercept)` + year * 1985)) %>%
    dplyr::select(dist_rate) %>%
    mutate(iter = paste0("iter_", 1:2000)) %>%
    mutate(error = 0)
  
  ### Introduce error and re-estimate
  
  trend_error_tmp <- vector("list", length(errors))
  intercept_error_tmp <- vector("list", length(errors))
  
  k <- 0
  
  for (error in errors) {
    
    k <- k + 1
    
    dist_observed_error <- dist_observed
    dist_observed_error[sample(which(dist_observed_error == 1), round((length(which(dist_observed_error == 1)) * error), 0))] <- 0
    dist_observed_error <- apply(dist_observed_error, 1, sum)
    dist_observed_error <- data.frame(disturbance = dist_observed_error, stable = forest - dist_observed_error, year = years)
    
    sim_dat[, paste0("disturbance_rate_", error)] <- dist_observed_error$disturbance / 500
    
    fit_error <- stan_glmer(cbind(disturbance, stable) ~ year + (1 | year),
                            data = dist_observed_error,
                            family = binomial(link = "logit"),
                            chains = 2, 
                            iter = 2000)
    
    trend_error_tmp[[k]] <- as.matrix(fit_error) %>%
      as.data.frame() %>%
      gather() %>%
      filter(key == "year") %>%
      mutate(iter = paste0("iter_", 1:2000)) %>%
      mutate(trend = (exp(value) - 1) * 100) %>%
      dplyr::select(-key, -value) %>%
      mutate(error = error)
    
    intercept_error_tmp[[k]] <- as.matrix(fit_error) %>%
      as.data.frame() %>%
      mutate(dist_rate = boot::inv.logit(`(Intercept)` + year * 1985)) %>%
      dplyr::select(dist_rate) %>%
      mutate(iter = paste0("iter_", 1:2000)) %>%
      mutate(error = error)
      
  }
  
  trend <- c(trend_error_tmp, list(trend_noerror)) %>%
    bind_rows() %>%
    mutate(trend_diff = trend - (exp(trend_true) - 1) * 100)
  
  trends[[sim_num]] <- trend
  
  intercept <- c(intercept_error_tmp, list(intercept_noerror)) %>%
    bind_rows() %>%
    mutate(dist_rate_diff = dist_rate - dist_rate_true)
  
  intercepts[[sim_num]] <- intercept
  
  p1 <- trend %>%
    group_by(error) %>%
    summarize(median = median(trend),
              lower = quantile(trend, 0.025),
              upper = quantile(trend, 0.975)) %>%
    ggplot(., aes(x = error, y = median)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    geom_point() +
    geom_hline(yintercept = (exp(trend_true) - 1) * 100, col = "red") +
    theme_minimal() +
    labs(x = "Proportion omitted disturbances", y = bquote("Trend in mortality rate (% "*yr^-1*")"), col = NULL, title = paste0("Simulation #", sim_num)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.border = element_rect(fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(size = 6))
  
  p2 <- intercept %>%
    group_by(error) %>%
    summarize(median = median(dist_rate),
              lower = quantile(dist_rate, 0.025),
              upper = quantile(dist_rate, 0.975)) %>%
    ggplot(., aes(x = error, y = median)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    geom_point() +
    geom_hline(yintercept = dist_rate_true, col = "red") +
    theme_minimal() +
    labs(x = "Year", y = bquote("Average mortality rate (% "*yr^-1*")"), col = NULL, title = "") +
    theme(legend.position = "none",
          plot.title = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.border = element_rect(fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(size = 6))
  
  p <- p1 + p2 + plot_layout(ncol = 2)
  
  ggsave(paste0("results/sensitivitiy_analysis/simulation_results_", sim_num, ".pdf"), plot = p, width = 5, height = 2.5)
  
  p <- sim_dat %>%
    gather(key = error, value = disturbance_rate_error, -year, -average, -disturbance_rate) %>%
    mutate(error_rate = as.double(gsub("disturbance_rate_", "", error))) %>%
    ggplot(.,  aes(x = year)) +
    geom_line(aes(y = average * 100), col = "grey") +
    geom_point(aes(y = disturbance_rate * 100), col = "red") +
    geom_point(aes(y = disturbance_rate_error * 100, col = error_rate)) +
    theme_minimal() +
    labs(x = "Year", y = bquote("Simulated mortality rate (% "*yr^-1*")"), col = NULL, title = paste0("Simulation #", sim_num)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.border = element_rect(fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(size = 6))
  
  ggsave(paste0("results/sensitivitiy_analysis/simulation_", sim_num, ".pdf"), plot = p, width = 2.5, height = 2.5)
  
}

trends <- trends %>%
  set_names(1:n_sim) %>%
  bind_rows(.id = "sim")

intercepts <- intercepts %>%
  set_names(1:n_sim) %>%
  bind_rows(.id = "sim")

p1 <- trends %>%
  group_by(error, sim) %>%
  summarize(median = median(trend_diff),
            lower = quantile(trend_diff, 0.025),
            upper = quantile(trend_diff, 0.975)) %>%
  mutate(sim = factor(sim, levels = 1:30)) %>%
  ggplot(., aes(x = error, y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  theme_minimal() +
  labs(x = "Proportion omitted disturbances", y = bquote("Difference to true value"), col = NULL, title = "Trend in mortality rate") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(size = 6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~sim)

p2 <- intercepts %>%
  group_by(error, sim) %>%
  summarize(median = median(dist_rate_diff),
            lower = quantile(dist_rate_diff, 0.025),
            upper = quantile(dist_rate_diff, 0.975)) %>%
  mutate(sim = factor(sim, levels = 1:30)) %>%
  ggplot(., aes(x = error, y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  theme_minimal() +
  labs(x = "Proportion omitted disturbances", y = bquote("Difference to true value"), col = NULL, title = "Average mortality rate") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(size = 6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~sim)

ggsave("results/sensitivitiy_analysis/sensitivity_analysis_trends.pdf", p1, width = 7.5, height = 7.5 / 6 * 5)
ggsave("results/sensitivitiy_analysis/sensitivity_analysis_rates.pdf", p2, width = 7.5, height = 7.5 / 6 * 5)



