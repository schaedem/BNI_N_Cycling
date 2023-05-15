library(tidyverse)
library(lubridate)
library(SNCAnalysis)
library(cowplot)

full_dat <- SNCAnalysis::full_dat %>%
  mutate(DON = TN - mg_kg_NO3N - mg_kg_NH4N,
         treatment = as.factor(treatment)) 
head(full_dat)

karama_df <- full_dat %>% filter(location=="Karama")

rubona_df <- full_dat %>% filter(location=="Rubona") %>%
  mutate(date = case_when(
    timepoint == 1 ~ "9/17/20",
    timepoint == 2 ~ "10/7/20",
    timepoint == 3 ~ "11/19/20",
    timepoint ==4 ~ "12/8/20",
    timepoint ==5 ~ "1/20/21",
    timepoint ==6 ~ "2/11/21"
  ))
  
## set palettes
brach_pal <- c("#A6CEE3","#1F78B4")
nap_pal <- c("#FB9A99", "#E31A1C")
maize_pal <- c("#B2DF8A","#33A02C")


graph_fun <- function(trtname, pal, df) {
  
  graph_df <- df %>%
    mutate(date = mdy(date)) %>%
    dplyr::select(NP, log_dea, timepoint, date, treatment) %>%
    dplyr::filter(grepl(trtname, treatment)) %>%
    group_by(timepoint) %>%
    mutate(mean_NP = mean(NP),
          se_NP = sd(NP)/sqrt(48),
          mean_dea = mean(log_dea),
          se_dea = sd(log_dea)/ sqrt(48)) %>%
    ungroup() %>%
    dplyr::select(-c(NP, log_dea)) %>%
    unique()
    
    np_plot <-
      ggplot(data = graph_df, aes(x=date, y=mean_NP, 
                                          group=treatment, color=treatment)) +
      geom_point(position=position_dodge(width=10)) +
      geom_line(position=position_dodge(width=10)) + 
      theme_classic() +
      geom_errorbar(aes(ymin = mean_NP - se_NP, ymax= mean_NP + se_NP), width=0.3, 
                    position = position_dodge((width=10))) +
      scale_color_manual(values=pal) +
      ylab("NP") +
      xlab("") +
      ylim(-0.3,0.25)
      
    
    dea_plot <- 
      ggplot(data = graph_df, aes(x=date, y=mean_dea, 
                                  group=treatment, color=treatment)) +
      geom_point(position=position_dodge(width=10)) +
      geom_line(position=position_dodge(width=10)) + 
      theme_classic() +
      geom_errorbar(aes(ymin = mean_dea - se_dea, ymax= mean_dea + se_dea), width=0.3, 
                    position = position_dodge((width=10))) +
      scale_color_manual(values=pal) +
      ylim(0, 2.8) +
      ylab("log(DEA)") +
      xlab("") 
    
    legend <- get_legend(
      # create some space to the left of the legend
        dea_plot + theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.box.margin = margin(0, 0, 0, 12))
    )
    
    done <- list("np_plot" = np_plot, "dea_plot"=dea_plot, "legend" = legend, "graph_df" = graph_df)
    return(done)
}

k_brach <- graph_fun("Brachiaria", brach_pal, karama_df)
k_nap <- graph_fun("apier", nap_pal, karama_df)
k_maize<- graph_fun("aize", maize_pal, karama_df)

r_brach <- graph_fun("Brachiaria", brach_pal, rubona_df)
r_nap <- graph_fun("apier", nap_pal, rubona_df)
r_maize<- graph_fun("aize", maize_pal, rubona_df)

#####Karama all enzyme plots ####
all_k_plots <- 
  cowplot::plot_grid(
  k_nap$np_plot + theme(legend.position = "none"),
  k_brach$np_plot + theme(legend.position = "none") + ylab(""),
  k_maize$np_plot + theme(legend.position = "none") + ylab("") + 
                    annotate("text", label =".", size=8, x= as.Date("2020-09-20"), 
                             y=0.04, color = "#33A02C", fontface="bold"),
  k_nap$dea_plot + theme(legend.position = "none"),
  k_brach$dea_plot + theme(legend.position = "none") + ylab(""),
  k_maize$dea_plot + theme(legend.position = "none") + ylab(""),
  ncol=3,
  labels="AUTO",
  label_size = 10
)

all_k_enzyme_final <-
  cowplot::plot_grid(
    all_k_plots,
    k_nap$legend,
    k_brach$legend,
    k_maize$legend,
    ncol=1,
    rel_heights = c(5,0.3,0.3,0.3)
  )

ggsave(all_k_enzyme_final, filename = "karama_enzyme_time_trt_plot.png", dpi=300)

#### Rubona all enzyme plots####
all_r_plots <- 
  cowplot::plot_grid(
    r_nap$np_plot + theme(legend.position = "none"),
    r_brach$np_plot + theme(legend.position = "none") + ylab(""),
    r_maize$np_plot + theme(legend.position = "none") + ylab("") + 
                        annotate(geom="text", x=as.Date("2020-10-10"), 
                                 y=-0.18, label = "*", size = 6, fontface="bold", color= "#33A02C"),
    r_nap$dea_plot + theme(legend.position = "none"),
    r_brach$dea_plot + theme(legend.position = "none") + ylab(""),
    r_maize$dea_plot + theme(legend.position = "none") + ylab("") +
      annotate(geom="text", x=as.Date("2020-10-10"), 
               y=1.5, label = ".", size = 8, fontface="bold", color= "#33A02C") ,
    ncol=3,
    labels="AUTO",
    label_size = 10
  )

all_r_enzyme_final <-
  cowplot::plot_grid(
    all_r_plots,
    r_nap$legend,
    r_brach$legend,
    r_maize$legend,
    ncol=1,
    rel_heights = c(5,0.3,0.3,0.3)
  )
all_r_enzyme_final
ggsave(all_r_enzyme_final, filename = "rubona_enzyme_time_trt_plot.png", dpi=300)
