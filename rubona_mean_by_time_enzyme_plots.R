library(tidyverse)
library(SNCAnalysis)

full_dat <- SNCAnalysis::full_dat %>%
  mutate(DON = TN - mg_kg_NO3N - mg_kg_NH4N,
         treatment = as.factor(treatment)) 
head(full_dat)

rubona_df <- full_dat %>% 
  filter(location=="Rubona") %>%
  mutate(date = ifelse(date == "9/17/20", "9/16/20",
                       ifelse(date == "10/7/20", "10/6/20",
                              ifelse(date == "11/18/20", "11/19/20",
                                     ifelse(date == "12/7/20", "12/8/20",
                                            ifelse(date=="1/19/21", "1/20/21",
                                                   ifelse(date == "2/10/21", "2/11/21",
                                                          date)))))),
         date = lubridate::mdy(date),
  )

unique(rubona_df$date)
#### Karama, Brachiaria ####
r_graph_df_brach <- rubona_df %>%
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

r_brach_np_plot <- 
  ggplot(data = r_graph_df_brach, aes(x=date, y=mean_NP, group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NP - se_NP, ymax= mean_NP + se_NP), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#A6CEE3","#1F78B4")) +
  ylab("NP") +
  xlab("") +
  ylim(-0.35,0.25)

r_brach_np_plot 

r_brach_dea_plot <-   
  ggplot(data = r_graph_df_brach, aes(x=date, y=mean_dea, group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_dea - se_dea, ymax= mean_dea + se_dea), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#A6CEE3","#1F78B4")) +
  ylim(0, 2.8) +
  ylab("log(DEA)") +
  xlab("") 
r_brach_dea_plot

brach_legend <- get_legend(
  # create some space to the left of the legend
  r_brach_dea_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

r_brach <- 
  cowplot::plot_grid(
    r_brach_np_plot + theme(legend.position = "none"),
    r_brach_dea_plot + theme(legend.position = "none"),
    labels = c("C","D"),
    label_size = 10
  )

r_brach_final <-
  cowplot::plot_grid(
    brach_legend,
    r_brach,
    rel_heights = c(0.3,4),
    ncol=1
  )
r_brach_final

# 
# RColorBrewer::brewer.pal(7, "Paired")
# RColorBrewer::display.brewer.pal(7, "Paired")
####

#### Karama, Napier ####
r_graph_df_nap <- rubona_df %>%
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

r_nap_np_plot <- 
  ggplot(data = r_graph_df_nap, aes(x=date, y=mean_NP, 
                                    group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NP - se_NP, ymax= mean_NP + se_NP), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#FB9A99", "#E31A1C")) +
  ylab("NP") +
  xlab("") +
  ylim(-0.35,0.25)
r_nap_np_plot

r_nap_dea_plot <-   
  ggplot(data = r_graph_df_nap, aes(x=date, y=mean_dea, group=treatment, 
                                    color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_dea - se_dea, ymax= mean_dea + se_dea), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#FB9A99","#E31A1C")) +
  ylim(0, 2.8) +
  ylab("log(DEA)") +
  xlab("") 
r_nap_dea_plot

nap_legend <- get_legend(
  # create some space to the left of the legend
  r_nap_dea_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

r_nap <- 
  cowplot::plot_grid(
    r_nap_np_plot + theme(legend.position = "none"),
    r_nap_dea_plot + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 10
  )

r_nap_final <-
  cowplot::plot_grid(
    nap_legend,
    r_nap,
    rel_heights = c(0.3,4),
    ncol=1
  )
r_nap_final

#### Karama, maize####
pal <- c("#B2DF8A","#33A02C")

r_graph_df_maize <- rubona_df %>%
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

r_maize_np_plot <- 
  ggplot(data = r_graph_df_maize, aes(x=date, y=mean_NP, 
                                      group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NP - se_NP, ymax= mean_NP + se_NP), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=pal) +
  annotate("text", label = "*", size=8, color ="#33A02C", x= lubridate::mdy("10/07/20"), y=-0.1) +
  ylab("NP") +
  xlab("") +
  ylim(-0.35,0.25)
r_maize_np_plot

r_maize_dea_plot <-   
  ggplot(data = r_graph_df_maize, aes(x=date, y=mean_dea, group=treatment, 
                                      color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_dea - se_dea, ymax= mean_dea + se_dea), width=0.3, 
                position = position_dodge((width=10))) +
  annotate("text", label = ".", size=10, color ="#33A02C", x= lubridate::mdy("10/07/20"), y=1.8) +
  scale_color_manual(values=pal) +
  ylim(0, 2.8) +
  ylab("log(DEA)") +
  xlab("") 
r_maize_dea_plot

maize_legend <- get_legend(
  # create some space to the left of the legend
  r_maize_dea_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

r_maize <- 
  cowplot::plot_grid(
    r_maize_np_plot + theme(legend.position = "none"),
    r_maize_dea_plot + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 10
  )

r_maize_final <-
  cowplot::plot_grid(
    maize_legend,
    r_maize,
    rel_heights = c(0.3,4),
    ncol=1
  )
r_maize_final




#####Karama all enzyme plots ####
all_r_plots <- 
  cowplot::plot_grid(
    r_nap_np_plot + theme(legend.position = "none"),
    r_brach_np_plot + theme(legend.position = "none") + ylab(""),
    r_maize_np_plot + theme(legend.position = "none") + ylab(""),
    r_nap_dea_plot + theme(legend.position = "none"),
    r_brach_dea_plot + theme(legend.position = "none") + ylab(""),
    r_maize_dea_plot + theme(legend.position = "none") + ylab(""),
    ncol=3,
    labels="AUTO",
    label_size = 10
  )

all_r_enzyme_final <-
  cowplot::plot_grid(
    all_r_plots,
    nap_legend,
    brach_legend,
    maize_legend,
    ncol=1,
    rel_heights = c(5,0.3,0.3,0.3)
  )

all_r_enzyme_final
