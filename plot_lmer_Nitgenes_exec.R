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

##### Nitrification genes, Karama ####
#amoA
kArch_amoA_FA <- gene_lmer_fun(Arch_amoA_FA, karama_df) #sig (nap/desmodium) 
kArch_amoA_FB <- gene_lmer_fun(Arch_amoA_FB, karama_df)
kArch_amoA_F <- gene_lmer_fun(Arch_amoA_F, karama_df) 
kGamma_amoA_F1R1 <- gene_lmer_fun(Gamma_amoA_F1R1, karama_df)
kGamma_amoA_F2R1 <- gene_lmer_fun(Gamma_amoA_F2R1, karama_df)
kGamma_amoA_F1R2<- gene_lmer_fun(Gamma_amoA_F1R2, karama_df)
kBeta_amoA <- gene_lmer_fun(Beta_amoA, karama_df)
#nxrB
kNitros_nxrB <- gene_lmer_fun(Nitrospira_nxrB, karama_df)
kNitrob_nxrB <- gene_lmer_fun(Nitrobacter_nxrB, karama_df)
#hao
kProteo_hao <- gene_lmer_fun(Proteo_hao, karama_df)
#Comammox, anammox
kcomaB <- gene_lmer_fun(comaB, karama_df)
khzocl <- gene_lmer_fun(Annamox_hzocl, karama_df)

##### Nitrification genes, Rubona ####
#amoA
rArch_amoA_FA <- gene_lmer_fun(Arch_amoA_FA, rubona_df) #sig (maize mono) 
rArch_amoA_FB <- gene_lmer_fun(Arch_amoA_FB, rubona_df)
rArch_amoA_F <- gene_lmer_fun(Arch_amoA_F, rubona_df) 
rGamma_amoA_F1R1 <- gene_lmer_fun(Gamma_amoA_F1R1, rubona_df)
rGamma_amoA_F2R1 <- gene_lmer_fun(Gamma_amoA_F2R1, rubona_df)
rGamma_amoA_F1R2 <- gene_lmer_fun(Gamma_amoA_F1R2, rubona_df)
rBeta_amoA <- gene_lmer_fun(Beta_amoA, rubona_df)
#nxrB
rNitros_nxrB <- gene_lmer_fun(Nitrospira_nxrB, rubona_df)
rNitrob_nxrB <- gene_lmer_fun(Nitrobacter_nxrB, rubona_df)
#hao
rProteo_hao <- gene_lmer_fun(Proteo_hao, rubona_df)
#comammox, anammox
rcomaB <- gene_lmer_fun(comaB, rubona_df)
rhzocl <- gene_lmer_fun(Annamox_hzocl, rubona_df)

##### amoA plot layout ####

kplots_amoA <- 
  cowplot::plot_grid(
  kArch_amoA_FA$plot + xlab(""), 
  # kArch_amoA_FB$plot+ xlab(""),
  # kArch_amoA_F$plot+ xlab(""),
  # kGamma_amoA_F1R1$plot+ xlab(bquote('copies'~g^-1~"soil")),
  # kGamma_amoA_F2R1$plot+ xlab(bquote('copies'~g^-1~"soil")),
  # kGamma_amoA_F1R2$plot+ xlab(bquote('copies'~g^-1~"soil")),
  kBeta_amoA$plot + xlab(""),
  labels = "AUTO",
  label_size = 10
)
kplots_amoA

ktitle <- ggdraw() + 
  draw_label(expression
             (underline("                             Karama                             ")), 
             fontface = "bold", vjust=0.5, hjust=0.5, size=11) +
  theme(plot.margin = margin(0,0,0,7))

k_final_amoA <- cowplot::plot_grid(ktitle,kplots_amoA, ncol=1, rel_heights=c(0.5,10)) +
  theme(plot.margin = margin(0,0,0,0))

rplots_amoA <-
  cowplot::plot_grid(
    rArch_amoA_FA$plot + xlab(bquote('copies'~g^-1~"soil")), 
    # rArch_amoA_FB$plot+ xlab(""),
    # rArch_amoA_F$plot+ xlab(""),
    # rGamma_amoA_F1R1$plot+ xlab(bquote('copies'~g^-1~"soil")),
    # rGamma_amoA_F2R1$plot+ xlab(bquote('copies'~g^-1~"soil")),
    # rGamma_amoA_F1R2$plot+ xlab(bquote('copies'~g^-1~"soil")),
    rBeta_amoA$plot + xlab(bquote('copies'~g^-1~"soil")),
    labels = c("C", "D"),
    label_size = 10
  )
rplots_amoA 

rtitle <- ggdraw() + 
  draw_label(expression
             (underline("                             Rubona                             ")), 
             fontface = "bold", vjust=0.5, hjust=0.5, size=11) +
  theme(plot.margin = margin(0,0,0,7))

r_final_amoA <- cowplot::plot_grid(rtitle, rplots_amoA, ncol = 1, rel_heights = c(0.5,10)) + 
  theme(plot.margin = margin(0,0,0,0))

bigtitle <- ggdraw() + 
  draw_label(label = "Ammonia monooxygenase (amoA)", fontface = "bold", x=0.5, hjust=0.5) +
  theme(plot.margin = margin(0,0,7,0))

both <- cowplot::plot_grid(k_final_amoA, r_final_amoA, nrow = 2)
both
#### get legend #####

get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}
legend <- get_only_legend(khzocl$plot + 
                            theme(legend.position="bottom", 
                                  legend.title = element_blank()))

#### final amoA plot ####
amoa_final <- cowplot::plot_grid(
  bigtitle, both, legend, ncol=1, rel_heights = c(0.35,10, 0.8))
amoa_final

ggsave("Fig2_amoA_lmer_new.png", dpi=300, height=10, width=8, unit="in", bg="white")

#### nxrB and hao plot layout ####
hao_label <- ggdraw() +
  draw_label(label = "Hydroxylamine oxidoreductase (hao)", fontface = "bold", x=0.5, hjust=0.5)

nxrB_label <- ggdraw() + 
  draw_label(label = "Nitrite oxidoreductase (nxrB) & \n Hydroxylamine oxidoreductase (hao)", 
             fontface = "bold", x=0.5, hjust=0.5)


###Karama###

kplots_nxrB_hao <- 
  cowplot::plot_grid(kNitros_nxrB$plot + xlab(""), 
                     kNitrob_nxrB$plot + xlab(bquote('copies'~g^-1~"soil")), 
                     kProteo_hao$plot + xlab(bquote('copies'~g^-1~"soil")),
                     nrow=2, labels="AUTO", label_size = 10)
kplots_nxrB_hao_final <-
  cowplot::plot_grid(ktitle,kplots_nxrB_hao, ncol=1, rel_heights = c(0.5,10))

###Rubona###

rplots_nxrB_hao <- 
  cowplot::plot_grid(rNitros_nxrB$plot + xlab(""), 
                     rNitrob_nxrB$plot + xlab(bquote('copies'~g^-1~"soil")), 
                     rProteo_hao$plot + xlab(bquote('copies'~g^-1~"soil")),
                     nrow=2, labels=c("D", "E", "F"), label_size=10)

rplots_nxrB_hao_final <-
  cowplot::plot_grid(rtitle, rplots_nxrB_hao, ncol=1, rel_heights = c(0.5,10))

both_hao_nxrB <- 
  cowplot::plot_grid(
  kplots_nxrB_hao_final, rplots_nxrB_hao_final, nrow=1
)

#### final nxrb/hao plot ####
final_hao_nxrB <-
  cowplot::plot_grid(
    nxrB_label, both_hao_nxrB, legend,ncol=1, rel_heights = c(1,10, 0.8)
  )

#### comammox and anammox layout ####
kplots_mox <- 
  cowplot::plot_grid(
    kcomaB$plot + xlab(""),
    khzocl$plot + xlab(""),
    nrow=1,
    labels = "AUTO",
    label_size = 10
  )

rplots_mox <- 
  cowplot::plot_grid(
    rcomaB$plot+ xlab(bquote('copies'~g^-1~"soil")),
    rhzocl$plot+ xlab(bquote('copies'~g^-1~"soil")),
    nrow=1,
    labels = c("C", "D"),
    label_size = 10
  )

kplots_mox_final <- 
  cowplot::plot_grid(
    ktitle, kplots_mox, ncol=1, rel_heights = c(0.5,10)
  )
  
rplots_mox_final <- 
  cowplot::plot_grid(
    rtitle, rplots_mox, ncol=1, rel_heights = c(0.5,10)
  )

moxtitle <- ggdraw() + 
  draw_label(label = "Comammox (comaB) & Hydrazine oxidoreductase (hzocl)", fontface = "bold", x=0.5, hjust=0.5)


both_mox <- 
  cowplot::plot_grid(
    kplots_mox_final, rplots_mox_final, ncol=1
  )
  
#### final comammox/anammox plot ####
mox_final <- cowplot::plot_grid(
  moxtitle, both_mox, legend, ncol=1, rel_heights = c(0.3,4, 0.8))
mox_final

#### nifH ####
knifh <- gene_lmer_fun(nifH_IGK3, karama_df)
knifh_plot <- knifh$plot + xlab(bquote('copies'~g^-1~"soil")) 
knifh_final <- cowplot::plot_grid(ktitle, knifh_plot, ncol=1, rel_heights = c(0.5,10))

rnifh <- gene_lmer_fun(nifH_IGK3, rubona_df)
rnifh_plot <- rnifh$plot + xlab(bquote('copies'~g^-1~"soil"))
rnifh_final <- cowplot::plot_grid(rtitle, rnifh_plot, ncol=1, rel_heights = c(0.5,10))

nifh_title <-  ggdraw() + 
  draw_label(label = "Nitrogen fixation (nifH)", fontface = "bold", x=0.5, hjust=0.5)

both_nifh <-   cowplot::plot_grid(
  knifh_final, rnifh_final, ncol=2, labels="AUTO", label_size=10
)

nifh_final <- cowplot::plot_grid(nifh_title, both_nifh, legend, ncol=1, rel_heights=c(0.3,4))

#Figure 1
kplots_amoA_fig2 <- 
  cowplot::plot_grid(
    kArch_amoA_FA$plot + xlab(""), 
    kBeta_amoA$plot + xlab(""),
    labels = "AUTO",
    label_size = 10
  )
kplots_amoA_fig2

rplots_amoA_fig2 <- 
  cowplot::plot_grid(
    rArch_amoA_FA$plot + xlab(bquote('copies'~g^-1~"soil")),
    rBeta_amoA$plot + xlab(bquote('copies'~g^-1~"soil")),
    labels =c ("C", "D"),
    label_zie = 10
  )
rplots_amoA_fig2

cowplot::plot_grid(kplots_amoA_fig2, rplots_amoA_fig2, nrow=2,
                   rel_heights = c(6, 15))
