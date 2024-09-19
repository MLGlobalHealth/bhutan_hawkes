# Goodness of fits
library(epihawkes)
library(tidyr)
library(ggplot2)
library(ggpubr)

#--------------------------------------------------------------------------------------
# Choose distribution
mu_term <- "sinusoidal_linear"
mu_fn <- mu_sinusoidal_linear

# Falciparum
#--------------------------------------------------------------------------------------
# Loads in real data
real_events <- readRDS(paste0("data/pf_onset_2016.rds"))
parameters <- as.list(readRDS("output/pf_fit_2016.rds"))
parameters$delay <- 15
real_event_times <- as.vector(real_events$time)

# Calculates integral
integral <- vector(length = length(real_event_times))
for (i in 2:length(real_event_times)){
  events_sub <- real_event_times[1:i]
  integral[i] = integral_intensity(events = events_sub, int_kernel = int_ray,
                                   parameters = parameters, mu_fn = mu_sinusoidal_linear,
                                   mu_diff_fn = mu_diff_sinusoidal_linear, 
                                   mu_int_fn = mu_int_sinusoidal_linear,
                                   print_level = 1)
}

taus <- diff(integral)
zks <-  1-exp(-taus)
sorted_zks <- sort(zks)
sorted_zks = sorted_zks[sorted_zks > 0]
n = length(sorted_zks)
k = 1:n
bk = (k-0.5)/n

ks_data_pf_ray <- data.frame(zk = sorted_zks,
                             bk = bk,
                             kernel = rep("ray"), length(sorted_zks))

# KS Plot
p_ks_pf <- ggplot(ks_data_pf_ray) + 
  geom_point(aes(zk, bk)) + 
  geom_abline(intercept = 0, slope = 1, col = 'black') +
  geom_abline(intercept = 0 + 1.36/n^0.5, slope = 1, col = 'black', linetype = "dashed") +
  geom_abline(intercept = 0 - 1.36/n^0.5, slope = 1, col = 'black', linetype = "dashed") +
  theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  xlab("Quantiles") + ylab("Cumulative Distribution Function")

quantiles = qbeta(sorted_zks, shape1 = k, shape2 = n-k+1)
qq_data_pf_ray <- data.frame(zk = sorted_zks,
                                qs = quantiles, 
                                li = sorted_zks + 1.96*(sorted_zks*(1-sorted_zks)/n)^0.5,
                                ui = sorted_zks - 1.96*(sorted_zks*(1-sorted_zks)/n)^0.5,
                                kernel = rep("ray"), length(sorted_zks))
p_qq_pf <- ggplot(qq_data_pf_ray) + 
  geom_point(aes(zk, qs)) + 
  geom_abline(intercept = 0, slope = 1, col = 'black') +
  geom_line(aes(zk, li), col = 'black', linetype = "dashed") +
  geom_line(aes(zk, ui), col = 'black', linetype = "dashed") +
  theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  xlab("Quantiles") + ylab("Empiracle Quantiles")

# Vivax RAY
#-------------------------------------------------------------------------------------------------------------
# Loads in real data
real_events <- readRDS(paste0("data/pv_onset_2016.rds"))
parameters <- as.list(readRDS("output/pv_fit_2016.rds"))
parameters$delay <- 15
real_event_times <- as.vector(real_events$time)

# Calculates integral
integral <- vector(length = length(real_event_times))
for (i in 2:length(real_event_times)){
  events_sub <- real_event_times[1:i]
  integral[i] = integral_intensity(events = events_sub, int_kernel = int_ray,
                                   parameters = parameters, mu_fn = mu_sinusoidal_linear,
                                   mu_diff_fn = mu_diff_sinusoidal_linear, 
                                   mu_int_fn = mu_int_sinusoidal_linear,
                                   print_level = 1)
}

taus <- diff(integral)
zks <-  1-exp(-taus)
sorted_zks <- sort(zks)
sorted_zks = sorted_zks[sorted_zks > 0]
n = length(sorted_zks)
k = 1:n
bk = (k-0.5)/n

ks_data_pv_ray <- data.frame(zk = sorted_zks,
                             bk = bk,
                             kernel = rep("ray"), length(sorted_zks))

# KS Plot
p_ks_pv <- ggplot(ks_data_pv_ray) + 
  geom_point(aes(zk, bk)) + 
  geom_abline(intercept = 0, slope = 1, col = 'black') +
  geom_abline(intercept = 0 + 1.36/n^0.5, slope = 1, col = 'black', linetype = "dashed") +
  geom_abline(intercept = 0 - 1.36/n^0.5, slope = 1, col = 'black', linetype = "dashed") +
  theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  xlab("Quantiles") + ylab("Cumulative Distribution Function")

quantiles = qbeta(sorted_zks, shape1 = k, shape2 = n-k+1)
qq_data_pv_ray <- data.frame(zk = sorted_zks,
                                qs = quantiles, 
                                li = sorted_zks + 1.96*(sorted_zks*(1-sorted_zks)/n)^0.5,
                                ui = sorted_zks - 1.96*(sorted_zks*(1-sorted_zks)/n)^0.5,
                                kernel = rep("ray"), length(sorted_zks))
p_qq_pv <- ggplot(qq_data_pv_ray ) + 
  geom_point(aes(zk, qs)) + 
  geom_abline(intercept = 0, slope = 1, col = 'black') +
  geom_line(aes(zk, li), col = 'black', linetype = "dashed") +
  geom_line(aes(zk, ui), col = 'black', linetype = "dashed") +
  theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  xlab("Quantiles") + ylab("Empiracle Quantiles")


#----------------------------------------------------------------------------------------------------


p2 <- ggarrange(p_ks_pf, p_ks_pv,
                p_qq_pf, p_qq_pv,                
                labels = "AUTO",
                common.legend = T, 
                legend = "bottom")
ggsave("figures/fig_s1_qof_2016.pdf", p2, width = 8, height = 8)
