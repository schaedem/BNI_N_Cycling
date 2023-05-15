library(tidyverse)
library(vegan)
library(reshape2)
library(bioDist)
library(SNCAnalysis)
#source("cooccurrence_permanova.R")

setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/02_Station_N_Cycling/final_data")

nice_mat <- read.table(file="nice_matrix.Rdata") 
meta <- read_csv("SNC_metadata_2020.csv") %>%
  mutate(growth = case_when(timepoint==1 ~ "anthesis",
                            timepoint ==2 ~ "regrowth",
                            timepoint == 3 ~ "anthesis",
                            timepoint ==4 ~ "regrowth",
                            timepoint == 5 ~ "anthesis",
                            timepoint == 6 ~ "regrowth"))

soil_dat <- SNCAnalysis::full_dat

full_meta <- meta %>%
  merge(soil_dat, by=c("sample", "timepoint", "location","treatment")) %>%
  dplyr::select(-c(gwc.y, pH.y, date.y, auc_dea, log_auc_dea, auc_np, block.y)) %>%
  rename(gwc = gwc.x,
         pH = pH.x,
         date = date.x,
         block = block.x,
         PMN = NH4N_net) %>%
  mutate(DON = TN - mg_kg_NH4N, mg_kg_NO3N)

env_dat <- full_meta %>%
  dplyr::select(-c(location, treatment,date, growth, block, 
            sample, timepoint, truesoil,mmol_kg_hr, TN))


nice_dist_mat <- nice_mat %>%
  arrange(sample) %>%
  pivot_wider(names_from=assay, values_from=log_conc_truesoil) %>%
  as.data.frame()

#separate out into locations
loc_dist_mat <- nice_dist_mat%>%
  merge(meta, by="sample") 

karama_dist_mat <- loc_dist_mat %>%
  filter(location=="Karama")%>%
  dplyr::select(-c(location, timepoint, pH, gwc, date, treatment, block, growth))

rubona_dist_mat <- loc_dist_mat %>%
  filter(location=="Rubona")%>%
  dplyr::select(-c(location, timepoint, pH, gwc, date, treatment, block, growth))

####format gene abundanace distance matrix ####
rownames(nice_dist_mat) <- nice_dist_mat$sample
nice_dist_mat <- nice_dist_mat[,2:21]
nice_dist_mat <- as.matrix(nice_dist_mat)

dist <- vegdist(nice_dist_mat, method="bray")

#Karama distance matrix
rownames(karama_dist_mat) <- karama_dist_mat$sample
karama_dist_mat <- karama_dist_mat[,2:21]
karama_dist_mat <- as.matrix(karama_dist_mat)

#Rubona distance matrix
rownames(rubona_dist_mat) <- rubona_dist_mat$sample
rubona_dist_mat <- rubona_dist_mat[,2:21]
rubona_dist_mat <- as.matrix(rubona_dist_mat)


####permanova testing: gene abundance/community structure ####
adonis_df <- nice_dist_mat %>%
  as.data.frame() %>%
  rownames_to_column("sample")

adonis_merge <- adonis_df %>% 
  as_tibble() %>%
  merge(meta, by="sample")

set.seed(123)
treatment_adonis <- adonis2(dist~treatment, data=adonis_merge, permutations=999)

set.seed(123)
#timepoint_adonis <- adonis2(dist ~ timepoint, data=adonis_merge, permutations=999)
#timepoint highly sig, R2=0.48792, p=0.0001

set.seed(123)
#takes a LONG time to run
#timepoint_loc_adonis <- adonis2(dist ~ timepoint*location*block*treatment, 
#                                data=adonis_merge, permutations=9999)
#timepoint_loc_adonis
#timepoint, location, and timepoint:location are sig
#timepoint: R2 =0.49 , p=0.0001
#location: R2 = 0.03, p = 0.0001
#timepoint:location: R2 = 0.01, p=0.0026

set.seed(123)
#growth_all_adonis <- adonis2(dist~as.factor(timepoint)*location*gwc, 
 #                            data=adonis_merge, permutations=9999)
#growth_all_adonis

#### pairwise comparisons to find sig differences between timepoints ####
pairwise_df <- numeric() #initialize empty df to put p-values

# pairwise_adonis <- function(val1, val2, adonis_merge) {
#   val1 <- val1
#   val2 <- val2
# 
#   filter_merge <- adonis_merge %>%
#     filter(timepoint == val1 | timepoint == val2)
# 
#   dist_filter <- filter_merge %>%
#     select(all_of(.[["sample"]])) %>%
#     as.dist()
# 
#   set.seed(123)
#   result <- adonis(dist_filter ~ timepoint, data=filter_merge)
# 
#   result_p <- result[["aov.tab"]][["Pr(>F)"]][1]
# 
#   return(result_p)
# 
# }
# # 
# pairwise_1_2 <- pairwise_adonis(1,2, adonis_merge) #0.03
# pairwise_1_3 <- pairwise_adonis(1,3, adonis_merge) #0.001
# pairwise_1_4 <- pairwise_adonis(1,4, adonis_merge) #0.001
# pairwise_1_5 <- pairwise_adonis(1,5, adonis_merge) #0.001
# pairwise_1_6 <- pairwise_adonis(1,6, adonis_merge) #0.001
# pairwise_2_3 <- pairwise_adonis(2,3, adonis_merge) #0.001
# pairwise_2_4 <- pairwise_adonis(2,4,adonis_merge) #0.001
# pairwise_2_5 <- pairwise_adonis(2,5, adonis_merge) #0.001
# pairwise_2_6 <- pairwise_adonis(2,6, adonis_merge) #0.001
# pairwise_3_4 <- pairwise_adonis(3,4, adonis_merge) #0.001
# pairwise_3_5 <- pairwise_adonis(3,5, adonis_merge) #0.009
# pairwise_3_6 <- pairwise_adonis(3,6, adonis_merge) #0.001
# pairwise_4_5 <- pairwise_adonis(4,5, adonis_merge) #0.006
# pairwise_4_6 <- pairwise_adonis(4,6,adonis_merge) #0.001
# pairwise_5_6 <- pairwise_adonis(5,6,adonis_merge) #0.001
# 
# pairwise_df[c("1-2", "1-3", "1-4", "1-5", "1-6", "2-3", "2-4", "2-5", "2-6",
#               "3-4", "3-5", "3-6", "4-5", "4-6", "5-6")] <-
#    c(pairwise_1_2, pairwise_1_3, pairwise_1_4,pairwise_1_5,pairwise_1_6,
#      pairwise_2_3, pairwise_2_4 ,pairwise_2_5 ,pairwise_2_6 ,pairwise_3_4 ,pairwise_3_5 ,
#      pairwise_3_6 ,pairwise_4_5 , pairwise_4_6 ,pairwise_5_6  )
# #
#  pairwise_p <- p.adjust(pairwise_df, method="BH") #adjust for multiple comparisons
#  all(pairwise_p < 0.05)
# 

# set.seed(222)
# loc_adonis <- adonis2(dist ~ location, data=adonis_merge)
# #location sig, R2=0.028, p=0.001
# set.seed(222)
# gwc_adonis <- adonis2(dist~gwc, data=adonis_merge)
# #gwc sig, R2=0.14665, p=0.001
# set.seed(222)
# treat_adonis <- adonis2(dist~treatment, data=adonis_merge)
# #treatment not sig, R2=0.0047, P=0.993
# set.seed(222)
# growth_adonis <- adonis2(dist~growth, data=adonis_merge)
# #growth sig, R2=0.03028, P=0.002; likely conflated with other env variables
# 

###NMDS plots - gene structure ####
 #set.seed(100192) 
 set.seed(123)
 nmds <- metaMDS(nice_dist_mat, distance="bray")

plot(nmds)

nmds_plot <- scores(nmds) 
assay_df <- read_csv("assay_list_2.csv") %>%
  select(assay, acronym, pathway)
head(assay_df)

species_scores <- nmds_plot$species %>%
  as_tibble(rownames = "assay") %>%
  merge(assay_df, by="assay")

plot_df <- nmds_plot$sites %>%
  as_tibble(rownames="sample") %>%
  merge(full_meta, by="sample")
 
# abund_nmds_growth <- 
#   ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=growth)) + 
#   geom_point(alpha=0.8) + 
#   theme_classic() +
#   scale_color_brewer(palette="Set2") +
#   theme(legend.title=element_blank())
# abund_nmds_growth

# site scores
abund_nmds_time <-
  ggplot(data=plot_df, aes(x=NMDS1, y=NMDS2, color=as.factor(timepoint))) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") +
  stat_ellipse()
abund_nmds_time

ellipse_nmds_time <- 
  ggplot(data=plot_df, aes(x=NMDS1, y=NMDS2)) +
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") + 
  stat_ellipse(aes(fill = as.factor(timepoint)), geom = "polygon",
               alpha=0.2) + 
  viridis::scale_fill_viridis(discrete="TRUE", name="Timepoint")
ellipse_nmds_time

# species scores
library(ggrepel)

#find which variables are sig correlated with the axes
#removed main effects (location, date, timepoint, etc to identify important env variables)
en_results <- envfit(nmds, env_dat, na.rm=TRUE) #everything is sig except poxc and pH

#extract envfit results
env_coord <- as.data.frame(scores(en_results, "vectors")) * ordiArrowMul(en_results)
#drop non-sig parameters
env_coord <- env_coord[c(1,3:7, 9:10), ]

nice_species_plot <- 
  ellipse_nmds_time + 
  geom_point(data=species_scores, aes(x=NMDS1, y=NMDS2, shape=pathway), 
             alpha = 1, color="darkblue", size=4, inherit.aes=FALSE) +
  guides(fill=guide_legend(title = "Timepoint")) +
  guides(shape=guide_legend(title = "Functional Pathway")) +
  theme_classic() +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=env_coord, size=1,
               alpha=0.6, color = "grey30") +
  geom_text(data=env_coord, aes(x=NMDS1, y=NMDS2), color="grey30",
            fontface = "bold", label=rownames(env_coord)) +
  geom_text_repel(data=species_scores, aes(x=NMDS1, y=NMDS2), color="grey30",
                  label=species_scores$assay)
nice_species_plot

abund_nmds_plant <-
  ggplot(data=plot_df, aes(x=NMDS1, y=NMDS2, color=treatment)) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  scale_color_brewer(palette="Dark2") +
  theme(legend.title=element_blank()) + 
  stat_ellipse()
abund_nmds_plant
ggsave("nice_abund_nmds_plant.png", dpi=300)

abund_nmds_loc <-
  ggplot(data=plot_df, aes(x=NMDS1, y=NMDS2, color=location)) +
  geom_point(alpha=0.4, size=3) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(legend.title=element_blank()) +
  stat_ellipse()
abund_nmds_loc

ggsave("nice_abund_nmds_loc.png", dpi=300)

# 
# abund_nmds_gwc <-
#   ggplot(data=plot_df, aes(x=NMDS1, y=NMDS2, color=gwc)) +
#   geom_point(alpha=0.8) +
#   theme_classic() +
#   viridis::scale_color_viridis(option="B", name="GWC")
# abund_nmds_gwc



#final plot - both locations
final_nice_abund<- 
  abund_nmds_time +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2), data=env_coord, size=1,
              alpha=0.6, color = "grey30") +
  geom_text(data=env_coord, aes(x=NMDS1 +0.004, y=NMDS2), color="grey30",
            fontface = "bold", label=rownames(env_coord)) +
   geom_point(data=species_scores, aes(x=NMDS1, y=NMDS2), position="identity",
              inherit.aes=FALSE, shape=18, size=4, alpha=0.5,
              color = "royalblue")
ggsave(dpi=300, "nice_abund_nmds_envfit.png")
