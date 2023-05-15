library(tidyverse)
library(SNCAnalysis)

full_dat <- SNCAnalysis::full_dat %>%
  mutate(DON = TN - mg_kg_NO3N - mg_kg_NH4N,
         treatment = as.factor(treatment)) 
head(full_dat)

karama_df <- full_dat %>% filter(location=="Karama") %>%
  mutate(date = lubridate::mdy(date),
         )

unique(karama_df$date)
#### Karama, Brachiaria ####
k_graph_df_brach <- karama_df %>%
  dplyr::select(NP, timepoint, date, log_dea, treatment) %>%
  dplyr::filter(grepl("Brachiaria", treatment)) %>%
  group_by(date) %>%
  mutate(mean_NP = mean(NP),
         se_NP = sd(NP)/sqrt(4),
         mean_dea = mean(log_dea),
         se_dea = sd(log_dea)/ sqrt(4)) %>%
  ungroup() %>%
  dplyr::select(-c(NP, log_dea)) %>%
  unique()

k_brach_np_plot <- 
  ggplot(data = k_graph_df_brach, aes(x=date, y=mean_NP, group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NP - se_NP, ymax= mean_NP + se_NP), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#A6CEE3","#1F78B4")) +
  ylab("NP") +
  xlab("") +
  ylim(-0.14,0.3)

k_brach_np_plot 

k_brach_dea_plot <-   
  ggplot(data = k_graph_df_brach, aes(x=date, y=mean_dea, group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_dea - se_dea, ymax= mean_dea + se_dea), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#A6CEE3","#1F78B4")) +
  ylim(0, 3.3) +
  ylab("log(DEA)") +
  xlab("") 
k_brach_dea_plot

brach_legend <- get_legend(
  # create some space to the left of the legend
  k_brach_dea_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

k_brach <- 
  cowplot::plot_grid(
    k_brach_np_plot + theme(legend.position = "none"),
    k_brach_dea_plot + theme(legend.position = "none"),
    labels = c("C","D"),
    label_size = 10
  )

k_brach_final <-
  cowplot::plot_grid(
    brach_legend,
    k_brach,
    rel_heights = c(0.3,4),
    ncol=1
  )
k_brach_final

# 
# RColorBrewer::brewer.pal(7, "Paired")
# RColorBrewer::display.brewer.pal(7, "Paired")
####

#### Karama, Napier ####
k_graph_df_nap <- karama_df %>%
  dplyr::select(NP, timepoint, date, log_dea, treatment) %>%
  dplyr::filter(grepl("Napier", treatment)) %>%
  group_by(date) %>%
  mutate(mean_NP = mean(NP),
         se_NP = sd(NP)/sqrt(4),
         mean_dea = mean(log_dea),
         se_dea = sd(log_dea)/ sqrt(4)) %>%
  ungroup() %>%
  dplyr::select(-c(NP, log_dea)) %>%
  unique()

k_nap_np_plot <- 
  ggplot(data = k_graph_df_nap, aes(x=date, y=mean_NP, 
                                    group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NP - se_NP, ymax= mean_NP + se_NP), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#FB9A99", "#E31A1C")) +
  ylab("NP") +
  xlab("") +
  ylim(-0.14,0.3)
k_nap_np_plot

k_nap_dea_plot <-   
  ggplot(data = k_graph_df_nap, aes(x=date, y=mean_dea, group=treatment, 
                                    color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_dea - se_dea, ymax= mean_dea + se_dea), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#FB9A99","#E31A1C")) +
  ylim(0, 3.3) +
  ylab("log(DEA)") +
  xlab("") 
k_nap_dea_plot

nap_legend <- get_legend(
  # create some space to the left of the legend
  k_nap_dea_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

k_nap <- 
  cowplot::plot_grid(
    k_nap_np_plot + theme(legend.position = "none"),
    k_nap_dea_plot + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 10
  )

k_nap_final <-
  cowplot::plot_grid(
    nap_legend,
    k_nap,
    rel_heights = c(0.3,4),
    ncol=1
  )
k_nap_final

#### Karama, maize####
pal <- c("#B2DF8A","#33A02C")

k_graph_df_maize <- karama_df %>%
  dplyr::select(NP, timepoint, date, log_dea, treatment) %>%
  dplyr::filter(grepl("aize", treatment)) %>%
  group_by(date) %>%
  mutate(mean_NP = mean(NP),
         se_NP = sd(NP)/sqrt(4),
         mean_dea = mean(log_dea),
         se_dea = sd(log_dea)/ sqrt(4)) %>%
  ungroup() %>%
  dplyr::select(-c(NP, log_dea)) %>%
  unique()

k_maize_np_plot <- 
  ggplot(data = k_graph_df_maize, aes(x=date, y=mean_NP, 
                                    group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NP - se_NP, ymax= mean_NP + se_NP), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=pal) +
  annotate("text", label = ".", size=10, color ="#33A02C", x= lubridate::mdy("9/17/20"), y=0.06) +
  ylab("NP") +
  xlab("") +
  ylim(-0.14,0.3)
k_maize_np_plot

k_maize_dea_plot <-   
  ggplot(data = k_graph_df_maize, aes(x=date, y=mean_dea, group=treatment, 
                                    color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_dea - se_dea, ymax= mean_dea + se_dea), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=pal) +
  ylim(0, 3.3) +
  ylab("log(DEA)") +
  xlab("") 
k_maize_dea_plot

maize_legend <- get_legend(
  # create some space to the left of the legend
  k_maize_dea_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

k_maize <- 
  cowplot::plot_grid(
    k_maize_np_plot + theme(legend.position = "none"),
    k_maize_dea_plot + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 10
  )

k_maize_final <-
  cowplot::plot_grid(
    maize_legend,
    k_maize,
    rel_heights = c(0.3,4),
    ncol=1
  )
k_maize_final




#####Karama all enzyme plots ####
all_k_plots <- 
  cowplot::plot_grid(
  k_nap_np_plot + theme(legend.position = "none"),
  k_brach_np_plot + theme(legend.position = "none") + ylab(""),
  k_maize_np_plot + theme(legend.position = "none") + ylab(""),
  k_nap_dea_plot + theme(legend.position = "none"),
  k_brach_dea_plot + theme(legend.position = "none") + ylab(""),
  k_maize_dea_plot + theme(legend.position = "none") + ylab(""),
  ncol=3,
  labels="AUTO",
  label_size = 10
)

all_k_enzyme_final <-
  cowplot::plot_grid(
    all_k_plots,
    nap_legend,
    brach_legend,
    maize_legend,
    ncol=1,
    rel_heights = c(5,0.3,0.3,0.3)
  )
all_k_enzyme_final
