# August 2023
# Author: Jess K Hopf
# Toy_model_plotting_v0

# plots model outputs from MA_Toy_model_v0.m

# Load packages
library(tidyverse); library(viridisLite); library(scales); library(Ecfun)
library(ggridges); library(ggpubr)

# Clear environment
rm(list=ls())

# Single data ------------

# load single data
Results_in <- read_csv("model_outputs/2023-10-03_open_FisherySqueeze_no_NoiseCorr_none_decline_1_V0results_5pn2000_BR.csv") %>%
  mutate(logRR = log(RR))

# Results_in <- read_csv("model_outputs/2023-11-13_closed_FisherySqueeze_no_NoiseCorr_space_decline_1_V0results_5pn2000_BR.csv") %>%
#   mutate(logRR = log(RR))


# get summary
Results_sum <- Results_in %>% group_by(var, measure, time) %>% 
  summarise(mean = mean(logRR), median = median(logRR), max = max(logRR), min = min(logRR),
            sd = sd(logRR), se=sd(logRR)/sqrt(n()), n = n()) %>% 
  mutate

# get 95% lines
line95 <- left_join(Results_sum,
                         Results_sum %>% summarise(mm = max(mean))) %>% 
               mutate(r25 = first(time[mean > mm*0.25]),
                      r50 = first(time[mean > mm*0.5]),
                      r75 = first(time[mean > mm*0.75]),
                      r95 = first(time[mean > mm*0.95])) %>% 
               filter(time == 1)



# plots

# colours:
col_vec = c("#462B77", "#489D89", "#BADC3B")
leg_vals = c("Before-after", "Inside-outside", "BACI")

# simulated RRs over time
ggplot(Results_sum, aes(time, mean)) +
  geom_line(aes(color = measure), size = 1.5) +
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = measure), alpha = 0.2) +
  # geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = measure)) +
  # geom_vline(data = line95, aes(xintercept = r95), colour = 'grey80') +
  # geom_vline(data = line95, aes(xintercept = r75), colour = 'grey60') +
  # geom_vline(data = line95, aes(xintercept = r50), colour = 'grey40') +
  # geom_vline(data = line95, aes(xintercept = r25), colour = 'grey20') +
  facet_grid(rows = vars(var)) +
  scale_color_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  scale_fill_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  theme_bw() +
  theme(legend.position="bottom") +
  ylab('logRR') + xlab('Years afer MPA implimentation')


# ggplot(Results_sum, aes(time, median)) +
#   geom_line(aes(color = measure)) +
#   geom_ribbon(aes(ymin = median-se, ymax = median+se, fill = measure),
#               alpha = 0.2) +
#   facet_grid(rows = vars(var), scales = 'free') +
#   theme_bw()+
#   ylab('median(logRR)')
# 
# 
# ggplot(Results_in %>% filter(time == 10)) +
#   geom_density(aes(x = logRR, fill = measure), alpha = 0.5) +
#   facet_grid(rows = vars(var), scales = 'free') +
#   theme_bw()


# multiple data ------------

# load data
Results_inM <- rbind(read_csv("model_outputs/2023-10-03_open_FisherySqueeze_no_NoiseCorr_none_decline_1_V0results_5pn2000_BR.csv") %>% mutate(connect = "open", pinknoise_lvl = 5, noise_corr = "none", scenario = 'baseline'),
                    read_csv("model_outputs/2023-10-03_open_FisherySqueeze_yes_NoiseCorr_none_decline_1_V0results_5pn2000_BR.csv") %>% mutate(connect = "open", pinknoise_lvl = 5, noise_corr = "none", scenario = 'Fsqueeze'),
                    read_csv("model_outputs/2023-10-03_open_FisherySqueeze_no_NoiseCorr_none_decline_0.5_V0results_5pn2000_BR.csv") %>% mutate(connect = "open", pinknoise_lvl = 5, noise_corr = "none", scenario = 'decline'),
                    read_csv("model_outputs/2023-10-03_closed_FisherySqueeze_no_NoiseCorr_none_decline_1_V0results_5pn2000_BR.csv") %>% mutate(connect = "closed", pinknoise_lvl = 5, noise_corr = "none", scenario = 'baseline'),
                    read_csv("model_outputs/2023-10-03_closed_FisherySqueeze_yes_NoiseCorr_none_decline_1_V0results_5pn2000_BR.csv") %>% mutate(connect = "closed", pinknoise_lvl = 5, noise_corr = "none", scenario = 'Fsqueeze'),
                    read_csv("model_outputs/2023-10-03_closed_FisherySqueeze_no_NoiseCorr_none_decline_0.5_V0results_5pn2000_BR.csv") %>% mutate(connect = "closed", pinknoise_lvl = 5, noise_corr = "none", scenario = 'decline'),
                    read_csv("model_outputs/2023-11-13_open_FisherySqueeze_no_NoiseCorr_space_decline_1_V0results_5pn2000_BR.csv") %>% mutate(connect = "open", pinknoise_lvl = 5, noise_corr = "space", scenario = 'baseline'),
                    read_csv("model_outputs/2023-11-13_open_FisherySqueeze_yes_NoiseCorr_space_decline_1_V0results_5pn2000_BR.csv") %>% mutate(connect = "open", pinknoise_lvl = 5, noise_corr = "space", scenario = 'Fsqueeze'),
                    read_csv("model_outputs/2023-11-13_open_FisherySqueeze_no_NoiseCorr_space_decline_0.5_V0results_5pn2000_BR.csv") %>% mutate(connect = "open", pinknoise_lvl = 5, noise_corr = "space", scenario = 'decline'),
                    read_csv("model_outputs/2023-11-13_closed_FisherySqueeze_no_NoiseCorr_space_decline_1_V0results_5pn2000_BR.csv") %>% mutate(connect = "closed", pinknoise_lvl = 5, noise_corr = "space", scenario = 'baseline'),
                    read_csv("model_outputs/2023-11-13_closed_FisherySqueeze_yes_NoiseCorr_space_decline_1_V0results_5pn2000_BR.csv") %>% mutate(connect = "closed", pinknoise_lvl = 5, noise_corr = "space", scenario = 'Fsqueeze'),
                    read_csv("model_outputs/2023-11-13_closed_FisherySqueeze_no_NoiseCorr_space_decline_0.5_V0results_5pn2000_BR.csv") %>% mutate(connect = "closed", pinknoise_lvl = 5, noise_corr = "space", scenario = 'decline')) %>%
      mutate(logRR = log(RR))


# get summary
Results_sumM <- Results_inM %>% group_by(var, measure, time, connect, pinknoise_lvl, noise_corr, scenario) %>% 
  summarise(mean = mean(logRR), median = median(logRR), max = max(logRR), min = min(logRR),
            sd = sd(logRR), se=sd(logRR)/sqrt(n()), n = n()) %>% 
  mutate(scenario = factor(scenario, levels = c('baseline', 'Fsqueeze','decline')),
         connect = factor(connect, levels = c('open', 'closed')))


# plots
# colours:
col_vec = c("#462B77", "#489D89", "#BADC3B")
leg_vals = c("Before-after", "Inside-outside", "BACI")


# simulated RRs over time
ggplot(Results_sumM %>% filter(var == 'biomass' && noise_corr == 'none'), 
       aes(time, mean)) +
  geom_line(aes(color = measure)) +
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = measure),
              alpha = 0.2) +
  facet_grid(rows = vars(connect), cols = vars(scenario)) +
  scale_color_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  scale_fill_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  theme_bw() +
  theme(legend.position="bottom") +
  ylab('logRR') + xlab('Years afer MPA implimentation') +
  ggtitle('Biomass: No noise correlation')

ggplot(Results_sumM %>% filter(var == 'biomass' && noise_corr == 'space'), 
       aes(time, mean)) +
  geom_line(aes(color = measure)) +
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = measure),
              alpha = 0.2) +
  facet_grid(rows = vars(connect), cols = vars(scenario)) +
  scale_color_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  scale_fill_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  theme_bw() +
  theme(legend.position="bottom") +
  ylab('logRR') + xlab('Years afer MPA implimentation') +
ggtitle('Biomass: Spatial noise correlation')




ggplot(Results_sumM %>% filter(var == 'density' && noise_corr == 'none'), 
       aes(time, mean)) +
  geom_line(aes(color = measure)) +
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = measure),
              alpha = 0.2) +
  facet_grid(rows = vars(connect), cols = vars(scenario)) +
  scale_color_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  scale_fill_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  theme_bw() +
  theme(legend.position="bottom") +
  ylab('logRR') + xlab('Years afer MPA implimentation') +
  ggtitle('Abundance: No noise correlation')

ggplot(Results_sumM %>% filter(var == 'density' && noise_corr == 'space'), 
       aes(time, mean)) +
  geom_line(aes(color = measure)) +
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = measure),
              alpha = 0.2) +
  facet_grid(rows = vars(connect), cols = vars(scenario)) +
  scale_color_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  scale_fill_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  theme_bw() +
  theme(legend.position="bottom") +
  ylab('logRR') + xlab('Years afer MPA implimentation') +
  ggtitle('Abundance: Spatial noise correlation')




ggplot(Results_sumM %>% filter(var == 'size' && noise_corr == 'none'), 
       aes(time, mean)) +
  geom_line(aes(color = measure)) +
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = measure),
              alpha = 0.2) +
  facet_grid(rows = vars(connect), cols = vars(scenario)) +
  scale_color_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  scale_fill_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  theme_bw() +
  theme(legend.position="bottom") +
  ylab('logRR') + xlab('Years afer MPA implimentation') +
  ggtitle('Size: No noise correlation')

ggplot(Results_sumM %>% filter(var == 'size' && noise_corr == 'space'), 
       aes(time, mean)) +
  geom_line(aes(color = measure)) +
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = measure),
              alpha = 0.2) +
  facet_grid(rows = vars(connect), cols = vars(scenario)) +
  scale_color_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  scale_fill_manual("Comparison metric", values = col_vec, labels = leg_vals) + 
  theme_bw() +
  theme(legend.position="bottom") +
  ylab('logRR') + xlab('Years afer MPA implimentation') +
  ggtitle('Size: Spatial noise correlation')

























