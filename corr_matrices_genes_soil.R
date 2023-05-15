library(tidyverse)
library(lubridate)
library(SNCAnalysis)
library(viridis)
library(corrplot)
library(Hmisc)
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/02_Station_N_Cycling/final_data")

#####formatting data for correlograms ####
assay <- read_csv("assay_list_2.csv") %>%
  select(primer_pair, acronym)

only_16S <- read_csv("16S_data.csv")[,2:4] %>%
  mutate(acronym = ifelse(grepl("Arch_16S", acronym), "16S_Arch", acronym))

final_nice <- read_csv("final_nice_std_data.csv") %>%
  select(sample, log_conc_16S_truesoil, log_conc_truesoil, pH,
         gwc, timepoint, location, primer_F, primer_R, assay) %>%
  mutate(pH = ifelse(is.na(pH), 5.14, pH),
         primer_pair = paste(primer_F, primer_R, sep=" / ")) %>%
  select(-c(primer_F, primer_R)) %>%
  merge(assay, by= "primer_pair")

primers <- final_nice %>%
  select(sample, log_conc_truesoil, acronym) %>%
  rbind(only_16S)

meta_dat <- full_dat %>%
  select(location, timepoint, sample)

full_dat <- full_dat %>%
  mutate(DON = TN - mg_kg_NH4N - mg_kg_NO3N) %>%
  rename(PMN = NH4N_net)

long_data <- primers %>% 
  rbind(only_16S) %>%
  merge(meta_dat, by="sample")

spread_primers <- spread(primers, key=acronym, value=log_conc_truesoil) 

##### summarize gene abundance data by locxtime ####
# tab_sum <- long_data %>%
#   group_by(location, timepoint, acronym) %>%
#   summarise(mean = mean(log_conc_truesoil), stdd = sd(log_conc_truesoil),
#             stde = (stdd/sqrt(28)))
# 
# write_csv(tab_sum, "nice_gene_abunds.csv")

#### wide dataframe for correlogram: everything ####
full_spread <- spread_primers %>%
  merge(full_dat, by="sample") %>%
  select(-c(block, timepoint, location, mmol_kg_hr, date, TN, truesoil,
            log_auc_dea, auc_dea, treatment, auc_np, sample)) %>%
  mutate(pH = ifelse(is.na(pH), 5.14, pH),
         `16S` = as.numeric(`16S`))

cor_mat <- cor(full_spread, method="spearman", use="complete.obs") %>%
  round(2)

#### cor plot for soil features x gene abundances ####
genes <- full_spread[,1:22]
genes <- genes %>%
  mutate(`16S` = ifelse(is.na(`16S`), mean(`16S`), `16S`))

soil <- full_spread[, 23:32]

gene_by_soil <- cor(genes, soil, method="spearman", use="complete.obs") %>%
  round(2)

gene_by_soil_2 <- rcorr(as.matrix(genes), as.matrix(soil),
                        type="spearman")

r_vals <- gene_by_soil_2$r
r_vals <- r_vals[23:32 ,1:22]
p_vals <- gene_by_soil_2$P
p_vals <- p_vals[23:32, 1:22]

soil_gene_cor <- 
  corrplot(r_vals, tl.srt = 45, tl.col="black",
         tl.cex = 0.7, p.mat = p_vals, 
         sig.level=0.05, insig='blank')

#### cor plot for soil features ####
soil_cor_mat <- rcorr(as.matrix(soil), type="spearman")
soil_p <- soil_cor_mat$P
soil_r <- soil_cor_mat$r

corrplot(soil_r, tl.srt = 45, tl.col="black", tl.cex=0.7, p.mat=soil_p, 
         type="upper", diag = FALSE,
         sig.level=0.05, insig='blank')

#### cor plot for gene abundances ####

genes_cor_mat <- rcorr(as.matrix(genes), type="spearman")
genes_p <- genes_cor_mat$P
genes_r <- genes_cor_mat$r

gene_cor_plot <- 
  corrplot(genes_r, tl.srt = 45, tl.col="black", tl.cex=0.7, p.mat=genes_p, 
         type="upper", diag = FALSE,
         sig.level=0.05, insig='blank')

