library(epihawkes)
library(tidyr)
library(ggplot2)
library(matrixStats)
library(dplyr)
library(varhandle)

#--------------------------------------------------------------------------------------
# Choose distribution
mu_term <- "sinusoidal_linear"
mu_fn <- mu_sinusoidal_linear

#--------------------------------------------------------------------------------------
# Loads in real data
real_events <- readRDS("data/pf_onset_2016.rds")
output  <- readRDS("output/pf_fit_2016.rds")
parameters <- output
parameters$delay <- 15

real_event_times <- real_events$time
imported_events <- readRDS("output/pf_imported_2016.rds")

T_max <- max(real_event_times)

real_event_times = floor(real_event_times)
tab_days <- data.frame(table(real_event_times))

tab_days$real_event_times <- unfactor(tab_days$real_event_times)
real_days <- left_join(data.frame(real_event_times =  0:ceiling(T_max)), tab_days)
real_days[is.na(real_days)] <- 0

real_imports_days <- floor(imported_events)
tab_days_imports <- data.frame(table(real_imports_days))
tab_days_imports$real_imports_days <- unfactor(tab_days_imports$real_imports_days)
real_days_imports <- left_join(data.frame(real_imports_days =  0:ceiling(T_max)), tab_days_imports)
real_days_imports[is.na(real_days_imports)] <- 0

#---------------------------------------------------------------------------------------
num_procs <- 10

events <- list(length = num_procs)
counts_list <- list(length = num_procs)

for (i in 1:num_procs){
  load(paste0("output/pf_sims_", num_procs, "_2016.Rdata"))
  days <- data.frame(events_days =  0:ceiling(T_max))
  for (j in 1:length(event_times)){
    events_days <- floor(event_times[[j]])
    events_days <- events_days[events_days < ceiling(T_max)]
    tab_days <- data.frame(table(events_days))
    tab_days$events_days <- unfactor(tab_days$events_days)
    days <- left_join(days, tab_days, by = "events_days")
  }
  events[[i]] <- days[,2:ncol(days)]
}


events_all <- do.call(cbind, events)
events_all[is.na(events_all)] <- 0
events_df <- data.frame(events_all)
events_df$times <- 0:ceiling(T_max)
events_long <- gather(events_df, -times, "key" = key, "value" = values)

data <- data.frame("t" = 0:ceiling(T_max),
                   "real" = real_days$Freq)

#---------------------------------------------------------------------------------------
# Gathers mu data

num_procs <- 10

importations <- list(length = num_procs)
for (i in 1:num_procs){
  days <- data.frame(events_days =  0:ceiling(T_max))
  for (j in 1:length(event_times_imp)){
    events_days <- floor(event_times[[j]])
    tab_days <- data.frame(table(events_days))
    tab_days$events_days <- unfactor(tab_days$events_days)
    days <- left_join(days, tab_days, by = "events_days")
  }
  importations[[i]] <- days[,2:ncol(days)]
}

importations_all <- do.call(cbind, importations)
importations_all[is.na(importations_all)] <- 0
importations_df <- data.frame(importations_all)
importations_df$times <- 0:ceiling(T_max)
importations_long <- gather(importations_df, -times, "key" = key, "value" = values)


data_importations <- data.frame("t" = 0:ceiling(T_max),
                                "real" = real_days_imports$Freq)
#--------------------------------------------------------------------------------------
# Make plots
p1 <- ggplot(events_long) +  geom_line(aes(times, values, group = key), col = "black", alpha = 0.1) + 
  geom_col(data = data, aes(t, real), col = "red") + 
  xlab("Time (days)") + ylab("Daily cases") + theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(0,ceiling(T_max))) + scale_y_continuous(expand = c(0, 0)) + 
  theme(plot.margin=unit(c(.2,.5,.2,.2),"cm"))

p2 <- ggplot(importations_long) +  geom_line(aes(times, values, group = key), col = "black", alpha = 0.1) + 
  geom_col(data = data_importations, aes(t, real), col = "red") + 
  xlab("Time (days)") + ylab("Daily importations") + theme_bw() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,ceiling(T_max))) + scale_y_continuous(expand = c(0, 0)) + 
  theme(plot.margin=unit(c(.2,.5,.2,.2),"cm"))


#--------------------------------------------------------------------------------------
# Loads in real data
real_events <- readRDS("data/pv_onset_2016.rds")
output  <- readRDS("output/pv_fit_2016.rds")
parameters <- output
parameters$delay <- 15

real_event_times <- real_events$time
imported_events <- readRDS("output/pv_imported_2016.rds")

T_max <- max(real_event_times)

real_event_times = floor(real_event_times)
tab_days <- data.frame(table(real_event_times))

tab_days$real_event_times <- unfactor(tab_days$real_event_times)
real_days <- left_join(data.frame(real_event_times =  0:ceiling(T_max)), tab_days)
real_days[is.na(real_days)] <- 0

real_imports_days <- floor(imported_events)
tab_days_imports <- data.frame(table(real_imports_days))
tab_days_imports$real_imports_days <- unfactor(tab_days_imports$real_imports_days)
real_days_imports <- left_join(data.frame(real_imports_days =  0:ceiling(T_max)), tab_days_imports)
real_days_imports[is.na(real_days_imports)] <- 0

#---------------------------------------------------------------------------------------
num_procs <- 10

events <- list(length = num_procs)
counts_list <- list(length = num_procs)

for (i in 1:num_procs){
  load(paste0("output/pv_sims_", num_procs, "_2016.Rdata"))
  days <- data.frame(events_days =  0:ceiling(T_max))
  for (j in 1:length(event_times)){
    events_days <- floor(event_times[[j]])
    events_days <- events_days[events_days < ceiling(T_max)]
    tab_days <- data.frame(table(events_days))
    tab_days$events_days <- unfactor(tab_days$events_days)
    days <- left_join(days, tab_days, by = "events_days")
  }
  events[[i]] <- days[,2:ncol(days)]
}


events_all <- do.call(cbind, events)
events_all[is.na(events_all)] <- 0
events_df <- data.frame(events_all)
events_df$times <- 0:ceiling(T_max)
events_long <- gather(events_df, -times, "key" = key, "value" = values)

data <- data.frame("t" = 0:ceiling(T_max),
                   "real" = real_days$Freq)

#---------------------------------------------------------------------------------------
# Gathers mu data

num_procs <- 1

importations <- list(length = num_procs)
for (i in 1:num_procs){
  days <- data.frame(events_days =  0:ceiling(T_max))
  for (j in 1:length(event_times_imp)){
    events_days <- floor(event_times[[j]])
    tab_days <- data.frame(table(events_days))
    tab_days$events_days <- unfactor(tab_days$events_days)
    days <- left_join(days, tab_days, by = "events_days")
  }
  importations[[i]] <- days[,2:ncol(days)]
}

importations_all <- do.call(cbind, importations)
importations_all[is.na(importations_all)] <- 0
importations_df <- data.frame(importations_all)
importations_df$times <- 0:ceiling(T_max)
importations_long <- gather(importations_df, -times, "key" = key, "value" = values)


data_importations <- data.frame("t" = 0:ceiling(T_max),
                                "real" = real_days_imports$Freq)
#--------------------------------------------------------------------------------------
# Make plots
p3 <- ggplot(events_long) +  geom_line(aes(times, values, group = key), col = "black", alpha = 0.1) + 
  geom_col(data = data, aes(t, real), col = "red") + 
  xlab("Time (days)") + ylab("Daily cases") + theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(0,ceiling(T_max))) + scale_y_continuous(expand = c(0, 0)) + 
  theme(plot.margin=unit(c(.2,.5,.2,.2),"cm"))

p4 <- ggplot(importations_long) +  geom_line(aes(times, values, group = key), col = "black", alpha = 0.1) + 
  geom_col(data = data_importations, aes(t, real), col = "red") + 
  xlab("Time (days)") + ylab("Daily importations") + theme_bw() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,ceiling(T_max))) + scale_y_continuous(expand = c(0, 0)) + 
  theme(plot.margin=unit(c(.2,.5,.2,.2),"cm"))

p <- ggarrange(p1, p2, p3, p4,                 
               common.legend = T, 
               legend = "bottom",
               labels = c("AUTO"))
ggsave("figures/fig_s2_daily_counts.tiff", p, width = 10, height = 10)
