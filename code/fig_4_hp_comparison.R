library(epihawkes)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)

pf_parameters <- as.list(readRDS("output/pf_fit_2016.rds"))
pf_parameters$delay = 15

pv_parameters <- as.list(readRDS("output/pv_fit_2016.rds"))
pv_parameters$delay = 15

t <- seq(0, 40, length.out=151)
y_pf <- unlist(lapply(t, ray_kernel, parameters = pf_parameters))
y_pv <- unlist(lapply(t, ray_kernel, parameters = pv_parameters))

df <- data.frame("t" = t, "Falciparum" = y_pf, "Vivax" = y_pv)
df_long <- gather(data = df, key = "Species", value = "value", -t)
p1 <- ggplot(df_long, aes(t, value, group = Species)) + 
  geom_line(aes(col = Species, linetype = Species)) + theme_bw() + 
  ylab("Kernel intensity") + xlab("Time since symptom onset") + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme(legend.position = "bottom", plot.margin=unit(c(.2,.5,.2,.2),"cm"))

ts <- 1:1000
mus_pf <- mu_sinusoidal_linear(ts, parameters = pf_parameters)
mus_pv <- mu_sinusoidal_linear(ts, parameters = pv_parameters)
df <- data.frame("t" = ts, "Falciparum" = mus_pf, "Vivax" = mus_pv)
df_long <- gather(data = df, key = "Species", value = "value", -t)
p2 <- ggplot(df_long, aes(t, value, group = Species)) + 
  geom_line(aes(col = Species, linetype = Species)) + theme_bw() + 
  ylab("Importation intensity") + xlab("Time") + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme(legend.position = "bottom", plot.margin=unit(c(.2,.5,.2,.2),"cm"))

print(sprintf("Max importations falciparum: %f", max((df_long %>% filter(Species == "Falciparum"))$value)))
print(sprintf("Max importations vivax: %f", max((df_long %>% filter(Species == "Vivax"))$value)))

p <- ggarrange(p1, p2, ncol = 2,
               common.legend = T, 
               legend = "bottom",
               labels = "AUTO")
ggsave("figures/fig_4_fig_fits.pdf", p, width = 12, height = 5)
