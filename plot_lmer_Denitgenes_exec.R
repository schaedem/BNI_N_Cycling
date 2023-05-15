library(SNCAnalysis)
library(tidyverse)
library(emmeans)
library(ggforce)
library(ggforestplot)
library(lmerTest)
library(multcomp)
library(plotrix) #for std error
library(broom.mixed) #tidy mixed effects model results
library(ggpubr) #for ggarange
library(gridExtra) #for grid.arrange
library(cowplot) #for plot_grid

setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/02_Station_N_Cycling/final_data")
source("genes_lmer_plot_fun.R")

#### Denitrification genes, Karama ####
## nirK ##
knirK_fun <- gene_lmer_fun(nirK_Fungi, karama_df)
knirK_flacu <- gene_lmer_fun(nirK_FlaCu, karama_df)
knirK_876 <- gene_lmer_fun(nirK_876, karama_df)
## norB ##
knorB_2 <- gene_lmer_fun(norB_2, karama_df)
kqnorB <- gene_lmer_fun(`qnorB_2F-5R`, karama_df)
## nosZ ##
knosZ_I <- gene_lmer_fun(nosZ_1F, karama_df)
knosZ_II <- gene_lmer_fun(nosZ_912F, karama_df)

#### Denitrification genes, Rubona ####
## nirK ##
rnirK_fun <- gene_lmer_fun(nirK_Fungi, rubona_df)
rnirK_flacu <- gene_lmer_fun(nirK_FlaCu, rubona_df)
rnirK_876 <- gene_lmer_fun(nirK_876, rubona_df)
## norB ##
rnorB_2 <- gene_lmer_fun(norB_2, rubona_df)
rqnorB <- gene_lmer_fun(`qnorB_2F-5R`, rubona_df)
## nosZ ##
rnosZ_I <- gene_lmer_fun(nosZ_1F, rubona_df)
rnosZ_II <- gene_lmer_fun(nosZ_912F, rubona_df)

#### nirK plot layout ####
kplots_nirk <- 
  cowplot::plot_grid(
    knirK_fun$plot + xlab(""),
    knirK_876$plot+ xlab(bquote('copies'~g^-1~"soil")),
    knirK_flacu$plot+ xlab(bquote('copies'~g^-1~"soil")),
    labels = "AUTO",
    label_size = 10
  )
kplots_nirk

ktitle <- ggdraw() + 
  draw_label(expression
             (underline("                             Karama                             ")), 
             fontface = "bold", vjust=0.5, hjust=0.5, size=10) +
  theme(plot.margin = margin(0,0,0,7))

k_final_nirk <- cowplot::plot_grid(ktitle,kplots_nirk, ncol=1, rel_heights=c(0.5,10))

rplots_nirk <-
  cowplot::plot_grid(
      rnirK_fun$plot + xlab(""),
      rnirK_876$plot+ xlab(bquote('copies'~g^-1~"soil")),
      rnirK_flacu$plot+ xlab(bquote('copies'~g^-1~"soil")),
      labels = c("D", "E","F"),
      label_size = 10
  )
  
rplots_nirk

rtitle <- ggdraw() + 
  draw_label(expression
             (underline("                             Rubona                             ")), 
             fontface = "bold", vjust=0.5, hjust=0.5, size=10) +
  theme(plot.margin = margin(0,0,0,7))

r_final_nirK <- cowplot::plot_grid(rtitle, rplots_nirk, ncol = 1, rel_heights = c(0.5,10))

bigtitle <- ggdraw() + 
  draw_label(label = "Nitrite reductase (nirK)", fontface = "bold", x=0.5, hjust=0.5)

both <- cowplot::plot_grid(k_final_nirk, r_final_nirK, ncol = 2)

#### get legend #####
get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}
legend <- get_only_legend(knirK_876$plot + theme(legend.position="bottom", 
                                                 legend.title = element_blank()))

#final nirK plot
nirk_final <- cowplot::plot_grid(
  bigtitle, both, legend, ncol=1, rel_heights = c(0.3,10, 0.8))

#### norB plot playout ####
#karama plots
kplots_norB <- 
  cowplot::plot_grid(
    knorB_2$plot + xlab(""),
    kqnorB$plot + xlab(""),
    labels = "AUTO",
    label_size = 10
  )
kplots_norB

k_final_norb <- cowplot::plot_grid(ktitle, kplots_norB, ncol = 1, rel_heights = c(0.5,10))

#rubona plots
rplots_norB <- 
  cowplot::plot_grid(
    rnorB_2$plot+ xlab(bquote('copies'~g^-1~"soil")),
    rnorB_2$plot+ xlab(bquote('copies'~g^-1~"soil")), 
    labels = c("C", "D"),
    label_size = 10
  )

r_final_norb <- cowplot::plot_grid(rtitle, rplots_norB, ncol = 1, rel_heights = c(0.5,10))

both_norb <- cowplot::plot_grid(k_final_norb, r_final_norb, ncol = 1)

norb_title <- ggdraw() + 
  draw_label(label = "Nitric oxide reductase (norB)", fontface = "bold", x=0.5, hjust=0.5)

norb_final <- cowplot::plot_grid(norb_title, both_norb, legend, ncol = 1, rel_heights = c(0.3,10, 0.8))

####nosZ plot layout ####

#karama plots
kplots_nosZ <- 
  cowplot::plot_grid(
    knosZ_I$plot + xlab(""),
    knosZ_II$plot + xlab(""),
    labels = "AUTO",
    label_size = 10
  )
kplots_nosZ

k_final_nosZ <- cowplot::plot_grid(ktitle, kplots_nosZ, ncol = 1, rel_heights = c(0.5,10))

#rubona plots
rplots_nosZ <- 
  cowplot::plot_grid(
    rnosZ_I$plot+ xlab(bquote('copies'~g^-1~"soil")),
    rnosZ_II$plot+ xlab(bquote('copies'~g^-1~"soil")),
    labels = c("C", "D"),
    label_size = 10
  )
rplots_nosZ

r_final_nosZ <- cowplot::plot_grid(rtitle, rplots_nosZ, ncol = 1, rel_heights = c(0.5,10))

nosZ_title <- ggdraw() + 
  draw_label(label = "Nitrous oxide reductase (nosZ)", fontface = "bold", x=0.5, hjust=0.5)

both_nosZ <- cowplot::plot_grid(k_final_nosZ, r_final_nosZ, ncol = 1)

nosZ_final <- cowplot::plot_grid(nosZ_title, both_nosZ, legend, ncol = 1, rel_heights = c(0.3,10, 0.8))




