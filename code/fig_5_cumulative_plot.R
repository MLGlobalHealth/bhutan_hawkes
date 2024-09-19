rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
print(gc())

library(epihawkes)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

#--------------------------------------------------------------------------------------
# Choose distribution
mu_term <- "sinusoidal_linear"
mu_fn <- mu_sinusoidal_linear

#--------------------------------------------------------------------------------------
# Pf plots
#--------------------------------------------------------------------------------------
# Loads in real data
real_events <- readRDS("data/pf_onset_2016.rds")
output  <- readRDS("output/pf_fit_2016.rds")
parameters <- output
parameters$delay <- 15

real_event_times <- real_events$time
imported_events <- readRDS("output/pf_imported_2016.rds")
  
T_max <- max(real_event_times)

real_intensity <- compute_intensity_function(events = real_event_times, 
                                             kernel = ray_kernel,
                                             parameters = parameters, 
                                             N = 5000, T_max = T_max, mu_fn = mu_fn)

parameters_mu <- parameters
parameters_mu$alpha <- 0
parameters_mu$delta <- 0
real_intensity_imported <- compute_intensity_function(events = imported_events, 
                                                      kernel = ray_kernel,
                                                      parameters = parameters_mu, 
                                                      N = 5000, T_max = T_max, mu_fn = mu_fn)

real_counts <- compute_count_function(events = real_event_times, N = 5000, T_max = T_max)
real_counts_imported <- compute_count_function(events = imported_events, N = 5000, T_max = T_max)

#---------------------------------------------------------------------------------------
num_procs <- 10

intensities_list <- list(length = num_procs)
counts_list <- list(length = num_procs)
intensities_list_imp <- list(length = num_procs)
counts_list_imp <- list(length = num_procs)

for (i in 1:num_procs){
  load(paste0("output/pf_sims_", num_procs, "_2016.Rdata"))
  intensities_tmp <- t(intensities)
  counts_tmp <- t(counts)
  
  intensities_tmp_imp <- t(intensities_imp)
  counts_tmp_imp <- t(counts_imp)
  
  intensities_list[[i]] <- intensities_tmp
  counts_list[[i]] <- counts_tmp
  
  intensities_list_imp[[i]] <- intensities_tmp_imp
  counts_list_imp[[i]] <- counts_tmp_imp
}

intensities <- do.call(cbind, intensities_list)
counts <- do.call(cbind, counts_list)

intensities_imp <- do.call(cbind, intensities_list_imp)
counts_imp <- do.call(cbind, counts_list_imp)

intensity_df <- data.frame(intensities)
intensity_df$times <- time_vec
intensity_df_long <- gather(intensity_df, -times, "key" = key, "value" = values)
intensity_df_long$rowid <- 1:length(intensity_df$intensities)

counts_df <- data.frame(counts)
counts_df$times <- time_vec
counts_df_long <- data.frame(gather(counts_df, -times, "key" = key, "value" = values))
counts_df_long$rowid <- 1:length(counts_df$counts)

intensity_df_imp <- data.frame(intensities_imp)
intensity_df_imp$times <- time_vec
intensity_df_long_imp <- gather(intensity_df_imp, -times, "key" = key, "value" = values)
intensity_df_long_imp$rowid <- 1:length(intensity_df_imp$intensities_imp)

counts_df_imp <- data.frame(counts_imp)
counts_df_imp$times <- time_vec
counts_df_long_imp <- data.frame(gather(counts_df_imp, -times, "key" = key, "value" = values))
counts_df_long_imp$rowid <- 1:length(counts_df_imp$counts_imp)

counts_df_long = counts_df_long[which(counts_df_long$key %in% c("X1", "X2")),]
counts_df_long_imp = counts_df_long_imp[which(counts_df_long_imp$key %in% c("X1", "X2")),] 

p1 <- ggplot(counts_df_long) + geom_line(aes(times, values, group=factor(key)), alpha = 0.1) +
  geom_line(data = counts_df_long_imp, aes(times, values, group=factor(key)), alpha = 0.1, col = 'blue') +
  geom_line(data = real_counts, aes(x = t, y = counts), col='red') + 
  geom_line(data = real_counts_imported, aes(x = t, y = counts), col='green') +
  xlab("Time (days)") + ylab("Case count") + theme_bw() +
  scale_x_continuous(expand = c(0, 0))

counts_end = unlist(counts_df[5001,1:(1000*num_procs)])
importations_end = rev(unlist(counts_df_imp[5001,1:(1000*num_procs)]))
ratio = importations_end/counts_end
print(sprintf("Prop importations falciparum: %f [95 CI %f - %f]",
              mean(ratio)*100, quantile(ratio, probs = 0.025)*100, 
              quantile(ratio, probs = 0.975)*100))

ggsave("figures/fig_5_cumulative_counts_2016_1.pdf", p1, 
       width = 12, height = 5.9)

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
print(gc())

#--------------------------------------------------------------------------------------
# Pv plots
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
library(epihawkes)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)


# Choose distribution
mu_term <- "sinusoidal_linear"
mu_fn <- mu_sinusoidal_linear
num_procs <- 10

# Loads in real data
real_events <- readRDS("data/pv_onset_2016.rds")
output  <- readRDS("output/pv_fit_2016.rds")
parameters <- output
parameters$delay <- 15

real_event_times <- real_events$time
imported_events <- readRDS("output/pv_imported_2016.rds")

T_max <- max(real_event_times)

real_intensity <- compute_intensity_function(events = real_event_times,
                                             kernel = ray_kernel,
                                             parameters = parameters,
                                             N = 5000, T_max = T_max, mu_fn = mu_fn)

parameters_mu <- parameters
parameters_mu$alpha <- 0
parameters_mu$delta <- 0
real_intensity_imported <- compute_intensity_function(events = imported_events,
                                                      kernel = ray_kernel,
                                                      parameters = parameters_mu,
                                                      N = 5000, T_max = T_max, mu_fn = mu_fn)

real_counts <- compute_count_function(events = real_event_times, N = 5000, T_max = T_max)
real_counts_imported <- compute_count_function(events = imported_events, N = 5000, T_max = T_max)

#---------------------------------------------------------------------------------------
intensities_list <- list(length = num_procs)
counts_list <- list(length = num_procs)
intensities_list_imp <- list(length = num_procs)
counts_list_imp <- list(length = num_procs)


for (i in 1:num_procs){
  load(paste0("output/pv_sims_", num_procs, "_2016.Rdata"))
  intensities_tmp <- t(intensities)
  counts_tmp <- t(counts)

  intensities_tmp_imp <- t(intensities_imp)
  counts_tmp_imp <- t(counts_imp)

  intensities_list[[i]] <- intensities_tmp
  counts_list[[i]] <- counts_tmp

  intensities_list_imp[[i]] <- intensities_tmp_imp
  counts_list_imp[[i]] <- counts_tmp_imp
}

intensities <- do.call(cbind, intensities_list)
counts <- do.call(cbind, counts_list)

intensities_imp <- do.call(cbind, intensities_list_imp)
counts_imp <- do.call(cbind, counts_list_imp)

intensity_df <- data.frame(intensities)
intensity_df$times <- time_vec
intensity_df_long <- data.frame(gather(intensity_df, -times, "key" = key, "value" = values))
intensity_df_long$rowid <- 1:length(intensity_df$intensities)

counts_df <- data.frame(counts)
counts_df$times <- time_vec
counts_df_long <- gather(counts_df, -times, "key" = key, "value" = values)
counts_df_long$rowid <- 1:length(counts_df$counts)

intensity_df_imp <- data.frame(intensities_imp)
intensity_df_imp$times <- time_vec
intensity_df_long_imp <- gather(intensity_df_imp, -times, "key" = key, "value" = values)
intensity_df_long_imp$rowid <- 1:length(intensity_df_imp$intensities_imp)

counts_df_imp <- data.frame(counts_imp)
counts_df_imp$times <- time_vec
counts_df_long_imp <- data.frame(gather(counts_df_imp, -times, "key" = key, "value" = values))
counts_df_long_imp$rowid <- 1:length(counts_df_imp$counts_imp)

counts_df_long = counts_df_long[which(counts_df_long$key %in% c("X1", "X2")),]
counts_df_long_imp = counts_df_long_imp[which(counts_df_long_imp$key %in% c("X1", "X2")),] 

p2 <- ggplot(counts_df_long) + geom_line(aes(times, values, group=factor(key)), alpha = 0.1) +
  geom_line(data = counts_df_long_imp, aes(times, values, group=factor(key)), alpha = 0.1, col = 'blue') +
  geom_line(data = real_counts, aes(x = t, y = counts), col='red') +
  geom_line(data = real_counts_imported, aes(x = t, y = counts), col='green') +
  xlab("Time (days)") + ylab("Case count") + theme_bw() +
  scale_x_continuous(expand = c(0, 0))

counts_end = unlist(counts_df[5001,1:(1000*num_procs)])
importations_end = rev(unlist(counts_df_imp[5001,1:(1000*num_procs)]))
ratio = importations_end/counts_end
print(sprintf("Prop importations vivax: %f [95 CI %f - %f]",
              mean(ratio)*100, quantile(ratio, probs = 0.025)*100,
              quantile(ratio, probs = 0.975)*100))

#p <- ggarrange(p1, p2,
#               common.legend = T,
#               legend = "bottom",
#               labels = "AUTO")

ggsave("figures/fig_5_cumulative_counts_2016_2.pdf", p2,
       width = 12, height = 5.9)
