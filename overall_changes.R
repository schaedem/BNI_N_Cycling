#setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/02_Station_N_Cycling/final_data")
library(tidyverse)
library(cowplot)
source("univariate_time_loc.R")
source("enzyme_mean_by_time_plots.R")

delta_calc <- function(loc_df) {
  sum_df <- loc_df %>%
    filter(timepoint == 1 | timepoint == 6) %>%
    mutate(bio_rep = str_sub(sample, 2)) %>%
    dplyr::select(c(NP, log_dea, mg_kg_NO3N, mg_kg_NH4N, DON, NH4N_net, 
                    treatment, block, timepoint, location, bio_rep, sample)) %>%
    mutate(mg_kg_N2O = (exp(log_dea)/1000000)*(44.013*1000)) %>% #convert DEA to same scale as NP
    mutate(log_N2O = log(mg_kg_N2O)) %>%
    mutate(log_N2O_1 = ifelse(timepoint ==1, log_N2O, NA),
           log_N2O_6 = ifelse(timepoint ==6, log_N2O, NA),
           NP_1 = ifelse(timepoint==1, NP, NA),
           NP_6 = ifelse(timepoint ==6, NP, NA),
           DEA_1 = ifelse(timepoint ==1, log_dea, NA),
           DEA_6 = ifelse(timepoint == 6, log_dea, NA),
           NO3_1 = ifelse(timepoint==1, mg_kg_NO3N, NA),
           NO3_6 = ifelse(timepoint==6, mg_kg_NO3N, NA),
           NH4_1 = ifelse(timepoint==1, mg_kg_NH4N, NA),
           NH4_6 = ifelse(timepoint==6, mg_kg_NH4N, NA),
           PMN_1 = ifelse(timepoint == 1, NH4N_net, NA),
           PMN_6 = ifelse(timepoint ==6, NH4N_net, NA),
           DON_1 = ifelse(timepoint ==1, DON, NA),
           DON_6 = ifelse(timepoint == 6, DON, NA)) %>%
    dplyr::select(c(treatment, block, bio_rep, log_N2O_1, log_N2O_6,
                    NP_1, NP_6, DEA_1, DEA_6, NO3_1, NO3_6,
                    NH4_1, NH4_6, PMN_1, PMN_6, DON_1, DON_6, location)) %>%
     group_by(bio_rep) %>%
     fill(log_N2O_1, log_N2O_6, NP_1, NP_6, DEA_1, DEA_6, NO3_1, NO3_6,
          NH4_1, NH4_6, PMN_1, PMN_6, DON_1, DON_6, .direction="updown") %>%
    ungroup() %>%
    unique() %>%
    # mutate(mean_NP1 = mean(NP_1),
    #        mean_DEA1 = mean(DEA_1),
    #        mean_NO31 = mean(NO3_1),
    #        mean_NH41 = mean(NH4_1),
    #        mean_NP6 = mean(NP_6),
    #        mean_DEA6 = mean(DEA_6),
    #        mean_NO36 = mean(NO3_6),
    #        mean_NH46 = mean(NH4_6))
      mutate(
           delta_N2O = log_N2O_6-log_N2O_1,
           delta_NP = NP_6 - NP_1,
           delta_DEA = DEA_6 - DEA_1, 
           delta_NO3 = NO3_6 - NO3_1,
           delta_NH4 = NH4_6 - NH4_1,
           delta_PMN = PMN_6 - PMN_1,
           delta_DON = DON_6 - DON_1) %>%
    dplyr::select(c(treatment, block, bio_rep, delta_NP, delta_DEA, delta_NO3, 
                    delta_NH4, delta_PMN, delta_DON, delta_N2O, location
                    ))
  
  sum_df$treatment<-factor(sum_df$treatment,
                       levels=c("Napier grass", "Napier + Desmodium", 
                                "Brachiaria cv. Mulato II", "Brachiaria + Desmodium",
                                "Desmodium distortum", 
                                "maize monoculture", "Maize + Desmodium"), ordered=TRUE)
  return(sum_df)
  
}

#calculate total changes 
delta_karama <- delta_calc(karama_df) 
delta_rubona <- delta_calc(rubona_df)
all_delta <- rbind(delta_karama, delta_rubona)

sum_delta <- all_delta %>%
  group_by(location, treatment) %>%
  summarise(mean_delta_NP = mean(delta_NP),
            mean_delta_DEA = mean(delta_DEA),
            mean_delta_NH4 = mean(delta_NH4),
            mean_delta_NO3 = mean(delta_NO3),
            mean_delta_N2O = mean(delta_N2O)) %>%
  #convert DEA to same scale as NP
  mutate(mean_mg_N2O_kg = (mean_delta_DEA/1000000)*44.013*1000)
            # mean_NH41 = mean(mean_NH41),
            # mean_NH46 = mean(mean_NH46),
            # mean_NO31 = mean(mean_NO31),
            # mean_NO36 = mean(mean_NO36),
            # mean_NP1 = mean(mean_NP1),
            # mean_NP6 = mean(mean_NP6),
            # mean_DEA1 = mean(mean_DEA1),
            # mean_DEA6 = mean(mean_DEA6))

# fold_change <- all_delta %>%
#   group_by(location, treatment) %>%
#   summarise(fold_NP = abs(mean_NP1 - mean_NP6)/((abs(mean_NP1)+abs(mean_NP6))/2),
#             fold_DEA = abs(mean_DEA1 - mean_DEA6)/((abs(mean_DEA1)+abs(mean_DEA6))/2),
#             fold_NO3 = abs(mean_NO31 - mean_NO36)/((abs(mean_NO31)+abs(mean_NO36))/2),
#             fold_NH4 = abs(mean_NH41 - mean_NH46)/((abs(mean_NH41)+abs(mean_NH46))/2)) %>%
#   ungroup()

#### graphing ####
#create df for graphing significance
trts <- unique(all_delta$treatment)
np_sig <- data.frame(location = c("Karama", "Rubona"), 
                     treatment = "maize monoculture", 
                     label = c("*", ""), 
                     p.value = c(0.007, 0.07),
                     ylab = c(0.37, 0.5))

#visualize abs change
#custom trt palette
nap_np <- all_delta %>%
  filter(grepl("Napier grass", treatment)) %>%
  group_by(location) %>%
  summarize(mean = mean(delta_NP))

pal <- c("#B00000", "#FFB4B3",  "#065F9B", "#C0E8FD", "#6A3D9A", "#006D00", "#CCF9A4")

#### NP delta plot ####
np_delta <- ggplot(data=all_delta, aes(x=treatment, y=delta_NP, fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~location) +
  scale_fill_manual(values = pal) +
  coord_flip() +
  theme_classic() +
  xlab("") +
  ylab(expression("mg"~NO[3]^-~""~"-N"~"kg"~soil^-1~hr^-1)) +
  geom_hline(yintercept = 0, linetype=2, color = "darkgrey", size=1) + 
  geom_hline(aes(yintercept = mean), data=nap_np, color = "#B00000", size=1) +
  geom_text(data=np_sig, aes(x=treatment, y=ylab, label=label), 
            size=8,
            color = "#006D00",
            fontface="bold",
            inherit.aes=FALSE) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title=element_text(hjust=0.5, size=15),
        strip.text = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  ggtitle(Delta~"NP") 
np_delta

####DEA delta plot ####
nap_dea <- all_delta %>%
  filter(grepl("Napier grass", treatment)) %>%
  group_by(location) %>%
  summarize(mean = mean(delta_N2O))

dea_delta <- ggplot(data=all_delta, aes(x=treatment, y=delta_N2O, fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~location) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = pal) +
  xlab("") +
  ylab(expression("log"~"mg"~N[2]~"O"~"kg"~soil^-1~hr^-1)) +
  geom_hline(yintercept = 0, linetype=2, color = "darkgrey", size=1) + 
  geom_hline(aes(yintercept = mean), data=nap_dea, color = "#B00000", size=1) +
  theme(legend.position = "bottom",
        legend.title = element_blank(), 
        plot.title=element_text(hjust=0.5, size=15),
        strip.text = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  ggtitle(Delta~"DEA")
dea_delta

#### Nitrate delta plot ####
nap_nit <- all_delta %>%
  filter(grepl("Napier grass", treatment)) %>%
  group_by(location) %>%
  summarize(mean = mean(delta_NO3))

nit_sig <- data.frame(location = c("Karama"), 
                      treatment = "maize monoculture", 
                      label = c("*"), 
                      p.value = 0.047,
                      ylab = -24)

no3_delta <- ggplot(data=all_delta, aes(x=treatment, y=delta_NO3, fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~location) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = pal) +
  xlab("") +
  ylab(expression("mg"~NO[3]^-~""~"-N"~"kg"~soil^-1))+
  geom_hline(yintercept = 0, linetype=2, color = "darkgrey", size=1) + 
  geom_hline(aes(yintercept = mean), data=nap_nit, color = "#B00000", size=1) +
  geom_text(data=nit_sig, aes(x=treatment, y=ylab, label=label), 
            size=8,
            color = "#006D00",
            fontface="bold",
            inherit.aes=FALSE) +
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5, size=15),
        strip.text = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  ggtitle(Delta~NO[3]^-~""~"-N")
no3_delta

#### NH4 delta plot ####
nap_amm <- all_delta %>%
  filter(grepl("Napier grass", treatment)) %>%
  group_by(location) %>%
  summarize(mean = mean(delta_NH4))


nh4_delta <- ggplot(data=all_delta, aes(x=treatment, y=delta_NH4, fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~location) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = pal) +
  xlab("") +
  ylab(expression("mg"~NH[4]^+~""~"-N"~"kg"~soil^-1)) +
  geom_hline(yintercept = 0, linetype=2, color="darkgrey", size=1) + 
  geom_hline(aes(yintercept = mean), data=nap_amm, color = "#B00000", size=1) +
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5, size=15),
        strip.text = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+ 
  ggtitle(Delta~NH[4]^+~""~"-N") 
nh4_delta #nap line overlaps zero (mean = -0.00358)

#### arrange plots ####
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Thesis/Rwanda/Lab work/02_Station_N_Cycling/final_data/Figures")
legend <- cowplot::get_legend(np_delta)

delta_enzyme <- cowplot::plot_grid(np_delta + theme(legend.position="none"), 
                   dea_delta + theme(legend.position="none"), 
                   labels = "AUTO",
                   ncol=1)
delta_enzyme_final <- plot_grid(delta_enzyme, legend,
                                ncol = 1,
                                rel_heights = c(6,0.5))

delta_enzyme_final
ggsave("Fig3_enzyme_delta_new.png", dpi=300, height=6.5, width=8, unit="in", bg="white")

delta_minN <- cowplot::plot_grid(no3_delta,
                                 nh4_delta,
                                 labels= "AUTO",
                                 ncol=1)
delta_minN_final <- plot_grid(delta_minN, legend,
                              ncol = 1,
                              rel_heights = c(6,0.5))
delta_minN_final
ggsave("Fig4_minN_delta_new.png", dpi=300, height=6.5, width=8, unit="in", bg="white")

####Karama stats ####
k_delta_NP <- univar_trt(delta_NP, delta_karama) 
k_delta_NP$tuk %>% filter(adj.p.value < 0.06)
#maize had higher overall increase in NP than Napier over sampling period
# A tibble: 3 × 7
# term      contrast                              null.value estimate conf.low conf.high adj.p.value
# <chr>     <chr>                                      <dbl>    <dbl>    <dbl>     <dbl>       <dbl>
#   1 treatment maize monoculture-Desmodium distortum          0    0.273   0.0607     0.486     0.00676
# 2 treatment maize monoculture-Napier grass                 0    0.224   0.0285     0.419     0.0183 
# 3 treatment maize monoculture-Napier + Desmodium           0    0.220   0.0247     0.415     0.0210 

k_delta_dea <- univar_trt(delta_DEA, delta_karama)
k_delta_dea$tuk%>% filter(grepl("Napier grass", contrast)) #nothing

k_delta_NO3 <- univar_trt(delta_NO3, delta_karama)
k_delta_NO3$tuk %>% filter(grepl("Napier grass", contrast)) #maize is sig
#maize had higher overall DECREASE in NO3 than Napier over sampling period
# term      contrast                                   null.value estimate conf.low conf.high adj.p.value
# <chr>     <chr>                                           <dbl>    <dbl>    <dbl>     <dbl>       <dbl>
#   1 treatment Napier grass-maize monoculture                      0     11.7    0.127      23.3      0.0465
# 2 treatment Brachiaria cv. Mulato II-maize monoculture          0     12.8    1.16       24.3      0.0253
# 3 treatment Napier + Desmodium-maize monoculture                0     12.9    1.31       24.5      0.0231
# 4 treatment Desmodium distortum-maize monoculture               0     13.5    0.859      26.1      0.0315
# 5 treatment Maize + Desmodium-maize monoculture                 0     13.9    2.29       25.5      0.0127

k_delta_NH4 <- univar_trt(delta_NH4, delta_karama)
k_delta_NH4$tuk %>% filter(grepl("Napier grass", contrast)) #nothing

###Rubona stats ####
r_delta_NP <- univar_trt(delta_NP, delta_rubona)
r_delta_NP$tuk%>% filter(grepl("Napier grass", contrast))
#maize had marginally higher overall increase in NP than Napier over sampling period
# <chr>     <chr>                                      <dbl>    <dbl>    <dbl>     <dbl>     <dbl>
#   1 treatment Napier grass-Brachiaria + Desmodium            0  0.0809   -0.147      0.308    0.900 
# 2 treatment Napier grass-Napier + Desmodium                0  0.0710   -0.157      0.299    0.943 
# 3 treatment Napier grass-Brachiaria cv. Mulato II          0  0.0395   -0.188      0.267    0.997 
# 4 treatment Napier grass-Desmodium distortum               0  0.00652  -0.221      0.234    1.00  
# 5 treatment Maize + Desmodium-Napier grass                 0  0.167    -0.0601     0.395    0.248 
# 6 treatment maize monoculture-Napier grass                 0  0.213    -0.0148     0.440    0.0765

r_delta_dea <- univar_trt(delta_DEA, delta_rubona)
r_delta_dea$tuk %>% filter(grepl("Napier grass", contrast)) #nothing

r_delta_NO3 <- univar_trt(delta_NO3, delta_rubona)
r_delta_NO3$tuk %>% filter(grepl("Napier grass", contrast)) #nothing

r_delta_NH4 <- univar_trt(delta_NH4, delta_rubona)
r_delta_NH4$tuk %>% filter(grepl("Napier grass", contrast)) #nothing

