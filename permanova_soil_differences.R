setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/02_Station_N_Cycling/final_data")
source("permanova_community_differences.R")

#### format soil matrices ####
soil_mat <- full_meta
k_meta <- full_meta %>% filter(location=="Karama")
r_meta <- full_meta %>% filter(location=="Rubona") %>% drop_na()

all_mat <- soil_mat %>%
  select(-c(timepoint, location, treatment, block, date, truesoil, TN, 
            growth, mmol_kg_hr)) %>%
  #filter(sample != 717) %>% # sample 717 prevented convergence due to high NPOC - go into dataset and check why
  column_to_rownames("sample")

all_mat_std <- decostand(all_mat, method = "standardize")
all_dist_mat <- as.matrix(all_mat_std)
all_dist_euc <- vegdist(all_dist_mat, method="euclidean", na.rm=TRUE)

k_soil_mat <- soil_mat %>%
  filter(location == "Karama") %>%
  select(-c(timepoint, location, treatment, block, date, truesoil, TN, 
            growth, mmol_kg_hr)) %>%
  column_to_rownames("sample")
k_mat_std <- decostand(k_soil_mat, method = "standardize")
k_dist_mat <- as.matrix(k_mat_std)
k_dist_euc <- vegdist(k_dist_mat, method = "euclidean", na.rm=TRUE)

r_soil_mat <- soil_mat %>%
  filter(location == "Rubona") %>%
  select(-c(timepoint, location, treatment, block, date, truesoil, TN, 
            growth, mmol_kg_hr)) %>%
   column_to_rownames("sample")
r_mat_std <- decostand(r_soil_mat, method = "standardize")
r_dist_mat <- as.matrix(r_mat_std)
r_dist_euc <- vegdist(r_dist_mat, method="euclidean", na.rm=TRUE)

#### run NMDS ####

set.seed(123)
all_soil_nmds <- metaMDS(all_dist_mat, distance = "euclidean", na.rm=TRUE,
                         autotransform = FALSE)

nmds_all_soil_scores <- scores(all_soil_nmds$points) %>%
  as_tibble(rownames = "sample")

all_soil_plot_df <- nmds_all_soil_scores %>%
  merge(full_meta, by="sample")

all_env_df <- nmds_all_soil_scores %>% merge(full_meta, by="sample")

all_env_results <- envfit(all_soil_nmds, all_env_df, na.rm=T)

all_env_coord <- all_env_results$vectors$arrows %>%
  as_tibble(rownames = "vector")

all_env_coord <- all_env_coord[c(4,5,7,9:13,15,17), ]

#Karama
set.seed(123)
k_soil_nmds <- metaMDS(k_dist_mat, distance = "euclidean", na.rm = TRUE,
                       autotransform = TRUE)
plot(k_soil_nmds)

k_nmds_soil_scores <- scores(k_soil_nmds$points) %>%
  as_tibble(rownames = "sample")

k_soil_plot_df <- k_nmds_soil_scores %>%
  merge(k_meta, by="sample")

k_env_df <- k_nmds_soil_scores %>% merge(k_meta, by="sample")

k_env_results <- envfit(k_soil_nmds, k_soil_plot_df, na.rm=T)

k_env_coord <- k_env_results$vectors$arrows %>%
  as_tibble(rownames = "vector")

k_env_coord <- k_env_coord[c(4,5,7,9:13,15,17), ]

#Rubona
set.seed(123)
r_soil_nmds <- metaMDS(r_dist_mat, distance = "euclidean", na.rm = TRUE,
                       autotransform = TRUE)
plot(r_soil_nmds)

r_nmds_soil_scores <- scores(r_soil_nmds$points) %>%
  as_tibble(rownames = "sample")

r_soil_plot_df <- r_nmds_soil_scores %>%
  merge(full_meta, by="sample")

r_env_df <- r_nmds_soil_scores %>% merge(r_meta, by="sample")

r_env_results <- envfit(r_soil_nmds, r_soil_plot_df, na.rm=T)

r_env_coord <- r_env_results$vectors$arrows %>%
  as_tibble(rownames = "vector") 
r_env_coord <- r_env_coord[c(4,5,7,9:13,15,17), ]
# NMDS1    NMDS2     r2 Pr(>r)    
# MDS1        1.00000  0.00000 1.0000  0.001 ***
#   MDS2        0.00000  1.00000 1.0000  0.001 ***
#   timepoint   0.98556 -0.16932 0.4604  0.001 ***
#   gwc         0.98741 -0.15815 0.4967  0.001 ***
#   pH          0.19167 -0.98146 0.0885  0.003 ** 
#   block       0.42969  0.90298 0.0193  0.199    
# NP          0.99171 -0.12853 0.6564  0.001 ***
#   mmol_kg_hr  0.61334 -0.78982 0.3595  0.001 ***
#   log_dea     0.70580 -0.70842 0.4516  0.001 ***
#   mg_kg_NO3N -0.90583 -0.42364 0.5931  0.001 ***
#   PMN        -0.03459 -0.99940 0.3813  0.001 ***
#   mg_kg_NH4N -0.99749  0.07085 0.4790  0.001 ***
#   mg_kg_POXC  0.35570 -0.93460 0.2222  0.001 ***
#   truesoil   -0.99240  0.12304 0.5008  0.001 ***
#   NPOC       -0.93386  0.35763 0.2134  0.001 ***
#   TN         -0.70533 -0.70888 0.5267  0.001 ***
#   DON        -0.61639 -0.78744 0.4717  0.001 ***

#### Plot NDMS soil features ####
#timepoint as color factor
soil_nmds_plot_time <-
  ggplot() +
  geom_point(data=all_soil_plot_df, aes(x=MDS1, y=MDS2, color=as.factor(timepoint)), 
             alpha=0.2, size=3) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") + 
  stat_ellipse(geom="polygon", aes(fill=as.factor(timepoint),x=MDS1, y=MDS2,), 
               data=all_soil_plot_df, alpha=0.2) +
  viridis::scale_fill_viridis(discrete="TRUE", name="Timepoint") +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=all_env_coord, size=1,
               alpha=0.6, color = "grey30") +
  geom_text_repel(data=all_env_coord, aes(x=NMDS1, y=NMDS2), color="grey30",
            fontface = "bold", label=all_env_coord$vector, 
            size=3.5, max.overlaps=15, min.segment.length = 0) +
  xlab("NMDS1") +
  ylab("NMDS2")
soil_nmds_plot_time

k_soil_nmds_plot_time <-
  ggplot() +
  geom_point(data=k_soil_plot_df, aes(x=MDS1, y=MDS2, 
                                      color=as.factor(timepoint)), alpha=0.2, size=3) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") + 
  stat_ellipse(geom="polygon", aes(fill=as.factor(timepoint),x=MDS1, y=MDS2,), 
               data=k_soil_plot_df, alpha=0.2) +
  viridis::scale_fill_viridis(discrete="TRUE", name="Timepoint") +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=k_env_coord, size=1,
               alpha=0.6, color = "grey30") +
  geom_text_repel(data=k_env_coord, aes(x=NMDS1, y=NMDS2), color="grey30",
                  fontface = "bold", label=k_env_coord$vector, 
                  size=3.5, max.overlaps=15, min.segment.length = 0) +
  xlab("NMDS1") +
  ylab("NMDS2")
k_soil_nmds_plot_time

r_soil_nmds_plot_time <- 
  ggplot() +
  geom_point(alpha=0.2, size=3, data=r_soil_plot_df, 
             aes(x=MDS1, y=MDS2, color=as.factor(timepoint))) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") + 
  stat_ellipse(geom="polygon", aes(fill=as.factor(timepoint),x=MDS1, y=MDS2,), 
               data=r_soil_plot_df, alpha=0.2) +
  viridis::scale_fill_viridis(discrete="TRUE", name="Timepoint")  +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=r_env_coord, size=1,
               alpha=0.6, color = "grey30") +
  geom_text_repel(data=r_env_coord, aes(x=NMDS1, y=NMDS2), color="grey30",
                  fontface = "bold", label=all_env_coord$vector, 
                  size=3.5, max.overlaps=15, min.segment.length = 0) +
  xlab("NMDS1") +
  ylab("NMDS2")
r_soil_nmds_plot_time

time_legend <- cowplot::get_legend(
  # create some space to the left of the legend
    soil_nmds_plot_time + theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 12))
)

#final time plot for soil nmds

cowplot::plot_grid(
  soil_nmds_plot_time + theme(legend.position = "none"),
  k_soil_nmds_plot_time + theme(legend.position = "none"),
  r_soil_nmds_plot_time + theme(legend.position = "none"),
  time_legend,
  labels = c("A", "B", "C", ""),
  label_size = 10,
  rel_widths = c(4,4,4,0.3)
)
ggsave("fig5_nmds.png", dpi=300, height = 8, width=12, unit="in", bg="white")
#location as color factor
all_soil_nmds_plot_loc <- 
ggplot(data=all_soil_plot_df, aes(x=MDS1, y=MDS2, color=location)) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(legend.title=element_blank()) +
  stat_ellipse() +
  xlab("NMDS1") +
  ylab("NMDS2")

all_soil_nmds_plot_loc
ggsave("soil_nmds_location.png")

#forage treatment as color factor
all_soil_nmds_plot_trt <- 
  ggplot(data=all_soil_plot_df, aes(x=MDS1, y=MDS2, color=as.factor(treatment))) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  scale_color_brewer(palette="Dark2") +
  theme(legend.title=element_blank()) +
  stat_ellipse() +
  xlab("NMDS1") +
  ylab("NMDS2")

k_soil_nmds_plot_trt <-
  ggplot(data=k_soil_plot_df, aes(x=MDS1, y=MDS2, color=as.factor(treatment))) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  scale_color_brewer(palette="Dark2") +
  theme(legend.title=element_blank()) +
  stat_ellipse() +
  theme(legend.position = "none") +
  xlab("NMDS1") +
  ylab("NMDS2")
k_soil_nmds_plot_trt

r_soil_nmds_plot_trt <-
  ggplot(data=k_soil_plot_df, aes(x=MDS1, y=MDS2, color=as.factor(treatment))) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  scale_color_brewer(palette="Dark2") +
  theme(legend.title=element_blank()) +
  stat_ellipse() +
  theme(legend.position = "none") +
  xlab("NMDS1") +
  ylab("NMDS2")
r_soil_nmds_plot_trt


plant_legend <- get_legend(
  # create some space to the left of the legend
  all_soil_nmds_plot_trt + theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 12))
)

all_trt_soil_plots <-
  cowplot::plot_grid(all_soil_nmds_plot_trt + theme(legend.position="none"),
            k_soil_nmds_plot_trt,
            r_soil_nmds_plot_trt,
            plant_legend,
            labels=c("A", "B", "C", ""),
            label_size = 10)
  
all_trt_soil_plots

####Significance testing: PERMANOVA ####
#all data combined
#Location, timepoint, block
set.seed(123)
# adonis2(all_dist_euc ~ location*block*timepoint, method = "euclidean", data=soil_mat, na.rm=TRUE)
# Df SumOfSqs      R2       F Pr(>F)    
# location                   1    341.8 0.10269 47.2004  0.001 ***
#   block                      1     33.5 0.01006  4.6252  0.001 ***
#   timepoint                  1    329.0 0.09886 45.4373  0.001 ***
#   location:block             1     38.4 0.01154  5.3050  0.001 ***
#   location:timepoint         1    197.3 0.05928 27.2481  0.001 ***
#   block:timepoint            1      9.9 0.00298  1.3682  0.230    
# location:block:timepoint   1     10.4 0.00313  1.4364  0.190    

##Plant treatment
set.seed(123)
# adonis2(all_dist_euc ~ treatment, method = "euclidean", data=soil_mat, na.rm=TRUE)
# Df SumOfSqs      R2      F Pr(>F)  
# treatment   6     87.6 0.02633 1.4784  0.033 *

#Karama only
set.seed(123)
# adonis2(k_dist_euc ~ timepoint*block, data=k_meta, na.rm=TRUE)
# adonis2(formula = k_dist_euc ~ timepoint * block, data = k_meta, na.rm = TRUE)
# Df SumOfSqs      R2       F Pr(>F)    
# timepoint         1    99.24 0.05978 10.8692  0.001 ***
#   block             1    59.63 0.03592  6.5304  0.001 ***
#   timepoint:block   1    12.85 0.00774  1.4077  0.187    
# Residual        163  1488.28 0.89655                   
# Total           166  1660.00 1.00000      
set.seed(123)
# adonis2(k_dist_euc ~ treatment, data=k_meta, na.rm=TRUE)
# adonis2(formula = k_dist_euc ~ treatment, data = k_meta, na.rm = TRUE)
# Df SumOfSqs     R2      F Pr(>F)   
# treatment   6   103.09 0.0621 1.7656  0.002 **
#   Residual  160  1556.91 0.9379                 
# Total     166  1660.00 1.0000  

#Rubona only
set.seed(123)
# adonis2(r_dist_euc ~ treatment, data=r_meta, na.rm=TRUE)
# Df SumOfSqs      R2      F Pr(>F)    
# treatment   6   127.09 0.07664 2.2272  0.001 ***
#   Residual  161  1531.16 0.92336                  
# Total     167  1658.25 1.00000  
set.seed(123)
adonis2(r_dist_euc ~ timepoint*block, data=r_meta, na.rm=TRUE)

#### Envfit: soil variables ####
#Tells us which variables correlate with the NMDS ordination axes




