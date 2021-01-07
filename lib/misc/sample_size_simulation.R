
library(tidyverse)

simulation <- function(population_n, mortality_rate, sample_n) {
  
  # Create population
  
  population <- rep(0, population_n)
  population[sample(1:length(population), as.integer(length(population) * mortality_rate))] <- 1
  
  # Sample and estimate 
  
  rate_est <- vector("double", 100)
  
  for (i in 1:100) {
    
    sample <- sample(population, sample_n)
    fit <- glm(cbind(sum(sample), length(sample)) ~ 1, family = "binomial")
    rate_est[i] <- exp(coef(fit)[1])
    
  }
  
  rate_est_mn <- mean(rate_est)
  rate_est_se <- sd(rate_est)
  
  return(data.frame(mean = rate_est_mn, se = rate_est_se))
  
}

# Population size

simdat <- expand.grid(population_n = seq(500, 5000000, length.out = 100), 
                      mortality_rate = 0.1, 
                      sample_n = 500)

sim_out <- vector("list", nrow(simdat))

for (i in 1:nrow(simdat)) {
  
  sim_out[[i]] <- simulation(simdat[i, 1], simdat[i, 2], simdat[i, 3])
  
}

sim_out <- do.call("rbind", sim_out)
sim_out <- cbind(simdat, sim_out)

p <- ggplot(sim_out) +
  geom_errorbar(aes(x = population_n, ymin = mean - se, ymax = mean + se)) +
  geom_point(aes(x = population_n, y = mean)) +
  labs(x = "Population size", y = "Mean and standard error")

ggsave("Desktop/sim01.png", p, width = 5.5, height = 2.5)

# Mortality rate

simdat <- expand.grid(population_n = 5000, 
                      mortality_rate = seq(0.0001, 0.01, 0.0001), 
                      sample_n = 500)

sim_out <- vector("list", nrow(simdat))

for (i in 1:nrow(simdat)) {
  
  sim_out[[i]] <- simulation(simdat[i, 1], simdat[i, 2], simdat[i, 3])
  
}

sim_out <- do.call("rbind", sim_out)
sim_out <- cbind(simdat, sim_out)

p <- ggplot(sim_out) +
  geom_errorbar(aes(x = mortality_rate, ymin = mean - se, ymax = mean + se,
                    col = (mean - se) > 0)) +
  geom_point(aes(x = mortality_rate, y = mean, col = (mean - se) > 0)) +
  labs(x = "Mortality rate", y = "Mean and standard error") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed")

ggsave("Desktop/sim02.png", p, width = 5.5, height = 3.5)

# Sample size

simdat <- expand.grid(population_n = 5000000, 
                      mortality_rate = 0.01, 
                      sample_n = seq(100, 10000, 100))

sim_out <- vector("list", nrow(simdat))

for (i in 1:nrow(simdat)) {
  
  sim_out[[i]] <- simulation(simdat[i, 1], simdat[i, 2], simdat[i, 3])
  
}

sim_out <- do.call("rbind", sim_out)
sim_out <- cbind(simdat, sim_out)

p <- ggplot(sim_out) +
  geom_errorbar(aes(x = sample_n, ymin = mean - se, ymax = mean + se)) +
  geom_point(aes(x = sample_n, y = mean)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Sample size", y = "Mean and standard error")

ggsave("Desktop/sim03.png", p, width = 5.5, height = 2.5)

