setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/02_Station_N_Cycling/final_data")

library(RColorBrewer)
library(glmmLasso)
library(SNCAnalysis)
library(lme4)
library(MASS)
library(nlme)
library(tidyverse)
library(PerformanceAnalytics)
library(lmerTest)  # p-values for MEMs based on the Satterthwaite approximation
library(psycho)    # mainly for an "analyze()" function
library(broom)     # for tidy results
library(knitr)     # beautifying tables
library(sjPlot)    # for visualising MEMs
library(effects)   # for visualising MEMs
library(report)    # for describing models
library(emmeans)   # for post-hoc analysis

source("nxrB_amoA_analysis.R")
source("LASSOprep_nxr_amoA.R")

#set reference group to Napier
merge_df$treatment <- as.factor(merge_df$treatment)
merge_df$timepoint <- as.factor(merge_df$timepoint)
merge_df$block <- as.factor(merge_df$block)
levels(merge_df$treatment)
merge_df$treatment <-  C(merge_df$treatment, contr.treatment, base=7)

merge_df <- merge_df %>%
  dplyr::select(-c(log_NH4N, log_TON)) 

#### DEA ####

##vector for all lambdas to try
grid <- 10^seq(15, 1, length=100)

#initialize vectors
BIC <- rep(Inf, length(grid))
AIC <- rep(Inf, length(grid))
Devianz_ma <- NULL
Coeff_ma <- NULL

j<-1

for (j in 1:length(BIC)){
  print(paste("Iteration ", j, sep=""))
  
  mod <- try(
    glmmLasso(log_dea~ 
                # `Nspira_amoA_F1R2` +
                # `Nspira_nxrB_Arch_amoA_FA` +
                # `Nspira_nxrB_Arch_amoA_FB` +
                # `Nspira_nxrB_Beta_amoA` +
                # `Nspira_nxrB_comaB` +
                # `Nspira_nxrB_Gamma_amoA_F1R1` +
                # `Nspira_nxrB_Gamma_amoA_F2R1` +
                # `Nbacter_nxrB_amoA_F1R2` +
                # `Nbacter_nxrB_Arch_amoA_FA` +
                # `Nbacter_nxrB_Arch_amoA_FB` +
                # `Nbacter_nxrB_Beta_amoA` +
                # `Nbacter_nxrB_comaB` +
                # `Nbacter_nxrB_Gamma_amoA_F1R1` +
                # `Nbacter_nxrB_Gamma_amoA_F2R1` +
                # `Nspira_nxrB_hzocl` +
                # `Nbacter_nxrB_hzocl` +
              #  `Nspira_nxrB_hao` +
              #  `Nbacter_nxrB_hao` +
              #  `nirK_876_16S_bact` +
              #  `nirK_876_16S_arch` +
              #  `nirK_FlaCu_16S_bact` +
              #  `nirK_FlaCu_16S_arch` +
              #  `nirK_fungi_16S_bact` +
              #  `nirK_fungi_16S_arch` +
                #`nirK_876_hao` +
                #`nirK_FlaCu_hao` +
               # `nirK_fungi_hao` +
              #  `nosZI_16S_bact` +
                #`nosZI_16S_arch` +
               # `nosZII_16S_bact` +
                #`nosZII_16S_arch` +
               # `nosZI_hao` +
                #`nosZII_hao`+
                Nitrobacter_nxrB +
                Nitrospira_nxrB +
                nirK_Fungi +
                nosZ_1F +
                `16S_arch` +
                norB_2 +
                Gamma_amoA_F1R1 +
                nirK_FlaCu +
                `qnorB_2F-5R` + 
                nirK_876 +
                Nitrobacter_nxrB+
                Arch_amoA_FB+
                comaB+
                Annamox_hzocl+
                Arch_amoA_F+
                Gamma_amoA_F2R1+
                Gamma_amoA_F1R2+
                #`16S_bact`+
                nosZ_912F +
                nifH_IGK3 +
                Beta_amoA +
                Arch_amoA_FA +
                Proteo_hao +
                location +
                gwc, #gwc as covariate  
              rnd = list(block=~1, timepoint=~1), 
              family = gaussian(link = "identity"),
              lambda=grid[j], data = merge_df, switch.NR=TRUE, final.re=TRUE))

  # code to make it continue anyway if an error occurs
  if(class(mod)!="try-error")
  {
    
    #save BIC, AIC
    BIC[j]<-mod$bic
    AIC[j]<-mod$aic
    
    #save coefficient outputs
    Coeff_ma<-cbind(Coeff_ma,mod$coefficients)
    
    #save error (deviance) values
    #y.hat<-predict(glm1,NP_P600_EightElectrodeROI)
    # Devianz_ma[j]<-sum(family$dev.resids(NP_P600_EightElectrodeROI$Amplitude,y.hat,wt=rep(1,length(y.hat))))
    
  }
}

final_bic_lambda <- grid[which.min(BIC)] 
#without ratios or 16S bact: 135.3048
final_aic_lambda <- grid[which.min(AIC)] 
#without ratios: 13.84886

#examine final models
bic_mod <- glmmLasso(log_dea~ 
                       # `Nspira_amoA_F1R2` +
                       # `Nspira_nxrB_Arch_amoA_FA` +
                       # `Nspira_nxrB_Arch_amoA_FB` +
                       # `Nspira_nxrB_Beta_amoA` +
                       # `Nspira_nxrB_comaB` +
                       # `Nspira_nxrB_Gamma_amoA_F1R1` +
                       # `Nspira_nxrB_Gamma_amoA_F2R1` +
                       # `Nbacter_nxrB_amoA_F1R2` +
                       # `Nbacter_nxrB_Arch_amoA_FA` +
                       # `Nbacter_nxrB_Arch_amoA_FB` +
                       # `Nbacter_nxrB_Beta_amoA` +
                       # `Nbacter_nxrB_comaB` +
                       # `Nbacter_nxrB_Gamma_amoA_F1R1` +
                       # `Nbacter_nxrB_Gamma_amoA_F2R1` +
                       # `Nspira_nxrB_hzocl` +
                       # `Nbacter_nxrB_hzocl` +
                       # `Nspira_nxrB_hao` +
                       # `Nbacter_nxrB_hao` +
                       # `nirK_876_16S_bact` +
                       # `nirK_876_16S_arch` +
                       # `nirK_FlaCu_16S_bact` +
                       # `nirK_FlaCu_16S_arch` +
                       # `nirK_fungi_16S_bact` +
                       # `nirK_fungi_16S_arch` +
                       # `nirK_876_hao` +
                       # `nirK_FlaCu_hao` +
                       # `nirK_fungi_hao` +
                       # `nosZI_16S_bact` +
                       # `nosZI_16S_arch` +
                       # `nosZII_16S_bact` +
                       # `nosZII_16S_arch` +
                       # `nosZI_hao` +
                       # `nosZII_hao`+
                       Nitrobacter_nxrB +
                       Nitrospira_nxrB +
                       nirK_Fungi +
                       nosZ_1F +
                       `16S_arch` +
                       norB_2 +
                       Gamma_amoA_F1R1 +
                       nirK_FlaCu +
                       `qnorB_2F-5R` + 
                       nirK_876 +
                       Nitrobacter_nxrB+
                       Arch_amoA_FB+
                       comaB+
                       Annamox_hzocl+
                       Arch_amoA_F+
                       Gamma_amoA_F2R1+
                       Gamma_amoA_F1R2+
                      # `16S_bact`+
                       nosZ_912F +
                       nifH_IGK3 +
                       Beta_amoA +
                       Arch_amoA_FA +
                       Proteo_hao +
                       location + 
                       gwc, 
                     rnd = list(block=~1, timepoint=~1), 
                     family = gaussian(link = "identity"),
                     lambda=final_bic_lambda, data = merge_df, switch.NR=TRUE, final.re=TRUE)

summary(bic_mod) 
#build model in lmer for clean summary
dea_lasso_bic <- lmer(log_dea ~  nosZ_912F + 
                        (1|block) + (1|timepoint), data=merge_df)
summary(dea_lasso_bic)
dea_lasso_null <- lmer(log_dea ~ location + nosZ_912F + (1|block) + (1|timepoint), data=merge_df)
anova(dea_lasso_bic, dea_lasso_null) #dea_lasso preferred

sjPlot::tab_model(dea_lasso_bic)

ggplot(merge_df, aes(x=`nosZII_16S_bact`, y=log_dea), group=location)+
  geom_point(aes(color=timepoint))+
  ggpubr::stat_cor(method="spearman") +
  scale_color_brewer(palette = "PRGn") +
  geom_smooth(method="lm") +
  theme_minimal() +
 # theme(legend.position="none") +
  ylab("log(DEA)") +
  xlab("nosZII/16S") 

#### NP model ####

##vector for all lambdas to try
grid <- 10^seq(100, 1, length=100)

#initialize vectors
BIC <- rep(Inf, length(grid))
AIC <- rep(Inf, length(grid))
Devianz_ma <- NULL
Coeff_ma <- NULL

j<-1

for (j in 1:length(BIC)){
  print(paste("Iteration ", j, sep=""))
  
  mod <- try(
    glmmLasso(NP~ 
                # `Nspira_amoA_F1R2` +
                # `Nspira_nxrB_Arch_amoA_FA` +
                # `Nspira_nxrB_Arch_amoA_FB` +
                # `Nspira_nxrB_Beta_amoA` +
                # `Nspira_nxrB_comaB` +
                # `Nspira_nxrB_Gamma_amoA_F1R1` +
                # `Nspira_nxrB_Gamma_amoA_F2R1` +
                # `Nbacter_nxrB_amoA_F1R2` +
                # `Nbacter_nxrB_Arch_amoA_FA` +
                # `Nbacter_nxrB_Arch_amoA_FB` +
                # `Nbacter_nxrB_Beta_amoA` +
                # `Nbacter_nxrB_comaB` +
                # `Nbacter_nxrB_Gamma_amoA_F1R1` +
                # `Nbacter_nxrB_Gamma_amoA_F2R1` +
                # `Nspira_nxrB_hzocl` +
                # `Nbacter_nxrB_hzocl` +
                # `Nspira_nxrB_hao` +
                # `Nbacter_nxrB_hao` +
                # `nirK_876_16S_bact` +
                # `nirK_876_16S_arch` +
                # `nirK_FlaCu_16S_bact` +
                # `nirK_FlaCu_16S_arch` +
                # `nirK_fungi_16S_bact` +
                # `nirK_fungi_16S_arch` +
                # `nirK_876_hao` +
                # `nirK_FlaCu_hao` +
                # `nirK_fungi_hao` +
                # `nosZI_16S_bact` +
                # `nosZI_16S_arch` +
                # `nosZII_16S_bact` +
                # `nosZII_16S_arch` +
                # `nosZI_hao` +
                # `nosZII_hao`+
                Nitrobacter_nxrB +
                Nitrospira_nxrB +
                nirK_Fungi +
                nosZ_1F +
                `16S_arch` +
                norB_2 +
                Gamma_amoA_F1R1 +
                nirK_FlaCu +
                `qnorB_2F-5R` + 
                nirK_876 +
                Nitrobacter_nxrB+
                Arch_amoA_FB+
                comaB+
                Annamox_hzocl+
                Arch_amoA_F+
                Gamma_amoA_F2R1+
                Gamma_amoA_F1R2+
               # `16S_bact`+
                nosZ_912F +
                nifH_IGK3 +
                Beta_amoA +
                Arch_amoA_FA +
                Proteo_hao +
                location + treatment + gwc, 
              rnd = list(block=~1, timepoint=~1), 
              family = gaussian(link = "identity"),
              lambda=grid[j], data = merge_df, switch.NR=TRUE, final.re=TRUE))
  
  # code to make it continue anyway if an error occurs
  if(class(mod)!="try-error")
  {
    
    #save BIC, AIC
    BIC[j]<-mod$bic
    AIC[j]<-mod$aic
    
    #save coefficient outputs
    Coeff_ma<-cbind(Coeff_ma,mod$coefficients)
    
    #save error (deviance) values
    #y.hat<-predict(glm1,NP_P600_EightElectrodeROI)
    # Devianz_ma[j]<-sum(family$dev.resids(NP_P600_EightElectrodeROI$Amplitude,y.hat,wt=rep(1,length(y.hat))))
    
  }
}

final_bic_lambda <- grid[which.min(BIC)] #100 with no ratios
final_aic_lambda <- grid[which.min(AIC)] #10 with no ratios

#examine final models
np_bic_mod <- glmmLasso(NP~ 
                       # `Nspira_amoA_F1R2` +
                       # `Nspira_nxrB_Arch_amoA_FA` +
                       # `Nspira_nxrB_Arch_amoA_FB` +
                       # `Nspira_nxrB_Beta_amoA` +
                       # `Nspira_nxrB_comaB` +
                       # `Nspira_nxrB_Gamma_amoA_F1R1` +
                       # `Nspira_nxrB_Gamma_amoA_F2R1` +
                       # `Nbacter_nxrB_amoA_F1R2` +
                       # `Nbacter_nxrB_Arch_amoA_FA` +
                       # `Nbacter_nxrB_Arch_amoA_FB` +
                       # `Nbacter_nxrB_Beta_amoA` +
                       # `Nbacter_nxrB_comaB` +
                       # `Nbacter_nxrB_Gamma_amoA_F1R1` +
                       # `Nbacter_nxrB_Gamma_amoA_F2R1` +
                       # `Nspira_nxrB_hzocl` +
                       # `Nbacter_nxrB_hzocl` +
                       # `Nspira_nxrB_hao` +
                       # `Nbacter_nxrB_hao` +
                       # `nirK_876_16S_bact` +
                       # `nirK_876_16S_arch` +
                       # `nirK_FlaCu_16S_bact` +
                       # `nirK_FlaCu_16S_arch` +
                       # `nirK_fungi_16S_bact` +
                       # `nirK_fungi_16S_arch` +
                       # `nirK_876_hao` +
                       # `nirK_FlaCu_hao` +
                       # `nirK_fungi_hao` +
                       # `nosZI_16S_bact` +
                       # `nosZI_16S_arch` +
                       # `nosZII_16S_bact` +
                       # `nosZII_16S_arch` +
                       # `nosZI_hao` +
                       # `nosZII_hao`+
                       Nitrobacter_nxrB +
                       Nitrospira_nxrB +
                       nirK_Fungi +
                       nosZ_1F +
                       `16S_arch` +
                       norB_2 +
                       Gamma_amoA_F1R1 +
                       nirK_FlaCu +
                       `qnorB_2F-5R` + 
                       nirK_876 +
                       Nitrobacter_nxrB+
                       Arch_amoA_FB+
                       comaB+
                       Annamox_hzocl+
                       Arch_amoA_F+
                       Gamma_amoA_F2R1+
                       Gamma_amoA_F1R2+
                      # `16S_bact`+
                       nosZ_912F +
                       nifH_IGK3 +
                       Beta_amoA +
                       Arch_amoA_FA +
                       Proteo_hao +
                       location + treatment + gwc , 
                     rnd = list(block=~1, timepoint=~1), 
                     family = gaussian(link = "identity"),
                     lambda=final_bic_lambda, data = merge_df, switch.NR=TRUE, final.re=TRUE)

summary(np_bic_mod) #nxrB, nifH, treatment
levels(merge_df$treatment)
np_lasso_bic <- lmer(NP ~ `nifH_IGK3` + `Nitrobacter_nxrB` +
                       treatment -1 + (1|block) + (1|timepoint), data=merge_df)

sjPlot::tab_model(np_lasso_bic)
anova(np_lasso_bic)

#aic retains a lot, including mostly non-significant terms

ggplot(merge_df, aes(x=nifH_IGK3, y=NP), group=location)+
  geom_point(aes(color=timepoint))+
  ggpubr::stat_cor(method="spearman") +
  scale_color_brewer(palette = "PRGn") +
  geom_smooth(method="lm") +
  theme_minimal() +
  # theme(legend.position="none") +
  ylab("NP") +
  xlab("nifH_IGK3 log10(copies g-1 soil)") 

