setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/02_Station_N_Cycling/final_data")
source("permanova_community_differences.R")
library(ggplot2)
library(cowplot)
library(rlang)
assay_df <- read_csv("assay_list_2.csv") %>%
  select(assay, acronym, pathway)
head(assay_df)
#Run NMDS for each site
#Karama
set.seed(123)
karama_nice_nmds <- metaMDS(karama_dist_mat, distance="bray")
plot(karama_nice_nmds)
#Rubona
set.seed(963)
rubona_nice_nmds <- metaMDS(rubona_dist_mat, distance = "bray")
plot(rubona_nice_nmds)

#### Extract scores - Karama ####
k_nmds_plot <- scores(karama_nice_nmds) 

k_species_scores <- k_nmds_plot$species %>%
  as_tibble(rownames = "assay") %>%
  merge(assay_df, by="assay") %>%
  select(assay, acronym, NMDS1, NMDS2, pathway) %>%
  mutate(assay = as.numeric(assay))

k_plot_df <- k_nmds_plot$sites %>%
  as_tibble(rownames="sample") %>%
  merge(full_meta, by="sample")

#### Extract scores - Rubona ####
r_nmds_plot <- scores(rubona_nice_nmds) 

r_species_scores <- r_nmds_plot$species %>%
  as_tibble(rownames = "assay") %>%
  merge(assay_df, by="assay") %>%
  select(assay, acronym, NMDS1, NMDS2, pathway)

r_plot_df <- r_nmds_plot$sites %>%
  as_tibble(rownames="sample") %>%
  merge(full_meta, by="sample")

####plot results ####

#Site plots
k_abund_nmds_time <-
  ggplot(data=k_plot_df, aes(x=NMDS1, y=NMDS2)) +
 # geom_point(aes(color=as.factor(timepoint)), alpha=0.4, size=3) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") + 
  stat_ellipse(aes(fill = as.factor(timepoint)), geom = "polygon",
               alpha=0.2) + 
  viridis::scale_fill_viridis(discrete="TRUE", name="Timepoint")
k_abund_nmds_time

r_abund_nmds_time <-
  ggplot(data=r_plot_df, aes(x=NMDS1, y=NMDS2)) +
  #geom_point(aes(color=as.factor(timepoint)), alpha=0.4, size=3) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") + 
  stat_ellipse(aes(fill = as.factor(timepoint)), geom = "polygon",
               alpha=0.2) + 
  viridis::scale_fill_viridis(discrete="TRUE", name="Timepoint")
r_abund_nmds_time

#### envfit arrows ####

#Karama

#find which variables are sig correlated with the axes
#removed main effects (location, date, timepoint, etc to identify important env variables)
k_env_dat <- full_meta %>%
  filter(location == "Karama") %>%
  select(-c(location, treatment,date, growth, block, 
            sample, timepoint, truesoil,mmol_kg_hr, TN))
k_en_results <- envfit(karama_nice_nmds, k_env_dat, na.rm=TRUE) #gwc, ph,poxc not sig

#extract envfit results
k_env_coord <- as.data.frame(scores(k_en_results, "vectors")) * ordiArrowMul(en_results)
#drop non-sig parameters
k_env_coord <- k_env_coord[c(3:7,9:10),]

#Rubona
r_env_dat <- full_meta %>%
  filter(location == "Rubona") %>%
  select(-c(location, treatment,date, growth, block, 
            sample, timepoint, truesoil,mmol_kg_hr, TN))
r_en_results <- envfit(rubona_nice_nmds, r_env_dat, na.rm=TRUE) #pH,poxc, npoc not sig

#extract envfit results
r_env_coord <- as.data.frame(scores(r_en_results, "vectors")) * ordiArrowMul(en_results)
#drop non-sig parameters
r_env_coord <- r_env_coord[c(1, 3:7,10),]

##### plot results: Karama, site ####

#final plot - both locations
k_nice_nmds_time <-
  ggplot(data=k_plot_df, aes(x=NMDS1, y=NMDS2, color=as.factor(timepoint))) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=k_env_coord, size=1,
               alpha=0.6, color = "grey30") +
  geom_text(data=k_env_coord, aes(x=NMDS1 +0.004, y=NMDS2), color="grey30",
            fontface = "bold", label=rownames(k_env_coord)) +
  geom_point(data=k_species_scores, aes(x=NMDS1, y=NMDS2), position="identity",
             inherit.aes=FALSE, shape=18, size=4, alpha=0.5,
             color = "royalblue")
ggsave(dpi=300, "karama_nice_abund_nmds_envfit.png")

#### plot results: Rubona, site ####
r_nice_nmds_time <-
  ggplot(data=r_plot_df, aes(x=NMDS1, y=NMDS2, color=as.factor(timepoint))) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=r_env_coord, size=1,
               alpha=0.6, color = "grey30") +
  geom_text(data=r_env_coord, aes(x=NMDS1 +0.004, y=NMDS2), color="grey30",
            fontface = "bold", label=rownames(r_env_coord)) +
  geom_point(data=r_species_scores, aes(x=NMDS1, y=NMDS2), position="identity",
             inherit.aes=FALSE, shape=18, size=4, alpha=0.5,
             color = "royalblue")
ggsave(dpi=300, "rubona_nice_abund_nmds_envfit.png")


### #arrange all plots together ####
library(cowplot)
nmds_envfit <- 
  cowplot::plot_grid(
    final_nice_abund + theme(legend.position="none"),
    k_nice_nmds_time+ theme(legend.position="none"),
    r_nice_nmds_time+ theme(legend.position="none"),
    labels = "AUTO",
    label_size = 10
  )
nmds_envfit
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  final_nice_abund + theme(
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(nmds_envfit, legend, nrow = 2, rel_widths = c(3, .4),
          rel_heights = c(3, .4))

#### Species plots ####
k_species_plot <- 
  k_abund_nmds_time +
  geom_point(data=k_species_scores, aes(x=NMDS1, y=NMDS2, shape=pathway), alpha=1, 
             size=4, color="darkblue", inherit.aes=FALSE) +
  # geom_text_repel(data=k_species_scores, aes(x=NMDS1, y=NMDS2), color="grey30",
  #                 label=species_scores$acronym) +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=k_env_coord, size=1,
             alpha=0.6, color = "grey30") +
  geom_text(data=k_env_coord, aes(x=NMDS1, y=NMDS2), color="grey30",
            fontface = "bold", label=rownames(k_env_coord)) +
  geom_text_repel(data=k_species_scores, aes(x=NMDS1, y=NMDS2), color="grey30",
                                  label=k_species_scores$assay) +
  guides(fill=guide_legend("Timepoint"),
         shape = guide_legend("Functional Pathway"))

r_species_plot <- 
  r_abund_nmds_time +
  geom_point(data=r_species_scores, aes(x=NMDS1, y=NMDS2, shape=pathway), alpha=1, 
             size=4, color = "darkblue", inherit.aes=FALSE) +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=r_env_coord, size=1,
               alpha=0.6, color = "grey30") +
  geom_text_repel(data=r_env_coord, aes(x=NMDS1, y=NMDS2+0.0005), color="grey30",
            fontface = "bold", label=rownames(r_env_coord)) +
  geom_text_repel(data=r_species_scores, aes(x=NMDS1, y=NMDS2), color="grey30",
                  label=r_species_scores$assay) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend("Timepoint"), 
         shape = guide_legend("Pathway"))

#### Arrange species plots ####
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  r_species_plot + theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 12))
)

all_species_plots <-
  cowplot::plot_grid(
    nice_species_plot + theme(legend.position="none"),
    k_species_plot + theme(legend.position = "none"),
    r_species_plot + theme(legend.position = "none"),
  labels = "AUTO",
  label_size = 10,
  nrow=2, ncol=2
)

final_species_plots <-
  plot_grid(
    all_species_plots, 
    legend,
     rel_widths = c(3,.4),
    # rel_hieghts = c(3,.4),
    ncol=2
  )
final_species_plots

ggsave("all_species_plots.png")
