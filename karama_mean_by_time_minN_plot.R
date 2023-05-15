library(tidyverse)
library(SNCAnalysis)

full_dat <- SNCAnalysis::full_dat %>%
  mutate(DON = TN - mg_kg_NO3N - mg_kg_NH4N,
         treatment = as.factor(treatment)) 
head(full_dat)

karama_df <- full_dat %>% filter(location=="Karama") %>%
  mutate(date = lubridate::mdy(date),
  )

#### Karama, Brachiaria ####
k_graph_df_brach <- karama_df %>%
  dplyr::select(mg_kg_NO3N, mg_kg_NH4N, timepoint, date, treatment) %>%
  dplyr::filter(grepl("Brachiaria", treatment)) %>%
  group_by(date) %>%
  mutate(mean_NO3N = mean(mg_kg_NO3N),
         se_NO3N = sd(mg_kg_NO3N)/sqrt(4),
         mean_NH4N = mean(mg_kg_NH4N),
         se_NH4N = sd(mg_kg_NH4N)/ sqrt(4)) %>%
  ungroup() %>%
  dplyr::select(-c(mg_kg_NO3N, mg_kg_NH4N)) %>%
  unique()

k_brach_no3_plot <- 
  ggplot(data = k_graph_df_brach, aes(x=date, y=mean_NO3N, group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NO3N - se_NO3N, ymax= mean_NO3N + se_NO3N), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#A6CEE3","#1F78B4")) +
  ylab(bquote('mg'~N0[3]^-~"" ~"-N"~g^-1~"soil")) +
  xlab("") +
  ylim(5, 30)

k_brach_no3_plot 

k_brach_nh4_plot <-   
  ggplot(data = k_graph_df_brach, aes(x=date, y=mean_NH4N, group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NH4N - se_NH4N, ymax= mean_NH4N + se_NH4N), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#A6CEE3","#1F78B4")) +
  annotate("text", label = "*", size=10, color ="#A6CEE3", x= lubridate::mdy("01/20/21"), y=1.3) +
  ylim(0, 4.8) +
  ylab(bquote('mg'~NH[4]^+~"" ~"-N"~g^-1~"soil")) +
  xlab("") 
k_brach_nh4_plot

brach_legend <- get_legend(
  # create some space to the left of the legend
  k_brach_nh4_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

k_brach <- 
  cowplot::plot_grid(
    k_brach_no3_plot + theme(legend.position = "none"),
    k_brach_nh4_plot + theme(legend.position = "none"),
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
  dplyr::select(mg_kg_NO3N, timepoint, date, mg_kg_NH4N, treatment) %>%
  dplyr::filter(grepl("Napier", treatment)) %>%
  group_by(date) %>%
  mutate(mean_NO3N = mean(mg_kg_NO3N),
         se_NO3N = sd(mg_kg_NO3N)/sqrt(4),
         mean_NH4N = mean(mg_kg_NH4N),
         se_NH4N = sd(mg_kg_NH4N)/ sqrt(4)) %>%
  ungroup() %>%
  dplyr::select(-c(mg_kg_NH4N, mg_kg_NO3N)) %>%
  unique()

k_nap_no3_plot <- 
  ggplot(data = k_graph_df_nap, aes(x=date, y=mean_NO3N, 
                                    group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NO3N - se_NO3N, ymax= mean_NO3N + se_NO3N), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#FB9A99", "#E31A1C")) +
  ylab(bquote('mg'~N0[3]^-~"" ~"-N"~g^-1~"soil")) +
  xlab("") +
  ylim(5,30)
k_nap_no3_plot

k_nap_nh4_plot <-   
  ggplot(data = k_graph_df_nap, aes(x=date, y=mean_NH4N, group=treatment, 
                                    color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NH4N - se_NH4N, ymax= mean_NH4N + se_NH4N), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=c("#FB9A99","#E31A1C")) +
  ylim(0, 4.8) +
  ylab(bquote('mg'~NH[4]^+~"" ~"-N"~g^-1~"soil")) +
  xlab("") 
k_nap_nh4_plot

nap_legend <- get_legend(
  # create some space to the left of the legend
  k_nap_nh4_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

k_nap <- 
  cowplot::plot_grid(
    k_nap_no3_plot + theme(legend.position = "none"),
    k_nap_nh4_plot + theme(legend.position = "none"),
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
  dplyr::select(mg_kg_NO3N, timepoint, date, mg_kg_NH4N, treatment) %>%
  dplyr::filter(grepl("aize", treatment)) %>%
  group_by(date) %>%
  mutate(mean_NO3N = mean(mg_kg_NO3N),
         se_NO3N = sd(mg_kg_NO3N)/sqrt(4),
         mean_NH4N = mean(mg_kg_NH4N),
         se_NH4N = sd(mg_kg_NH4N)/ sqrt(4)) %>%
  ungroup() %>%
  dplyr::select(-c(mg_kg_NO3N, mg_kg_NH4N)) %>%
  unique()

k_maize_no3_plot <- 
  ggplot(data = k_graph_df_maize, aes(x=date, y=mean_NO3N, 
                                      group=treatment, color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NO3N - se_NO3N, ymax= mean_NO3N + se_NO3N), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=pal) +
  annotate("text", label = "*", size=10, color ="#33A02C", x= lubridate::mdy("10/15/20"), y=25) +
  annotate("text", label = "*", size=10, color = "#B2DF8A", x=lubridate::mdy("10/05/20"), y=25) +
  ylab(bquote('mg'~N0[3]^-~"" ~"-N"~g^-1~"soil")) +
  xlab("") +
  ylim(5,30)
k_maize_no3_plot

k_maize_nh4_plot <-   
  ggplot(data = k_graph_df_maize, aes(x=date, y=mean_NH4N, group=treatment, 
                                      color=treatment)) +
  geom_point(position=position_dodge(width=10)) +
  geom_line(position=position_dodge(width=10)) + 
  theme_classic() +
  geom_errorbar(aes(ymin = mean_NH4N - se_NH4N, ymax= mean_NH4N + se_NH4N), width=0.3, 
                position = position_dodge((width=10))) +
  scale_color_manual(values=pal) +
  annotate("text", label = "*", size=10, color ="#B2DF8A", x= lubridate::mdy("01/20/21"), y=1.3) +
  annotate("text", label = "*", size =10, color = "#33A02C", x=lubridate::mdy("1/28/21"), y=1.3) +
  ylim(0, 4.8) +
  ylab(bquote('mg'~NH[4]^+~"" ~"-N"~g^-1~"soil")) +
  xlab("") 
k_maize_nh4_plot

maize_legend <- get_legend(
  # create some space to the left of the legend
  k_maize_nh4_plot + theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin = margin(0, 0, 0, 12))
)

k_maize <- 
  cowplot::plot_grid(
    k_maize_no3_plot + theme(legend.position = "none"),
    k_maize_nh4_plot + theme(legend.position = "none"),
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
    k_nap_no3_plot + theme(legend.position = "none"),
    k_brach_no3_plot + theme(legend.position = "none") + ylab(""),
    k_maize_no3_plot + theme(legend.position = "none") + ylab(""),
    k_nap_nh4_plot + theme(legend.position = "none"),
    k_brach_nh4_plot + theme(legend.position = "none") + ylab(""),
    k_maize_nh4_plot + theme(legend.position = "none") + ylab(""),
    ncol=3,
    labels="AUTO",
    label_size = 10
  )

all_k_minN_final <-
  cowplot::plot_grid(
    all_k_plots,
    nap_legend,
    brach_legend,
    maize_legend,
    ncol=1,
    rel_heights = c(5,0.3,0.3,0.3)
  )
all_k_minN_final
