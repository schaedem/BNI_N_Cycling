library(SNCAnalysis)
library(tidyverse)
library(emmeans)
library(ggforce)
library(ggforestplot)
library(lmerTest)
library(multcomp)
library(plotrix) #for std error
library(broom.mixed) #tidy mixed effects model results
library(dotwhisker)
source("nxrB_amoA_analysis.R")
source("LASSOprep_nxr_amoA.R")
#set reference group to Napier - redo with lasso_df
#merge_df uses scaled, centered data so units are meaningless
# merge_df$timepoint <- as.factor(merge_df$timepoint)
# merge_df$treatment <- as.factor(merge_df$treatment)
# merge_df$treatment<-factor(merge_df$treatment,
#                             levels=c("Napier grass", "Napier + Desmodium", 
#                                      "Brachiaria cv. Mulato II", "Brachiaria + Desmodium",
#                                      "Desmodium distortum", 
#                                      "maize monoculture", "Maize + Desmodium"),
#                             ordered=TRUE)
# 
# karama_df <- merge_df %>% filter(location=="Karama") 
# rubona_df <- merge_df %>% filter(location=="Rubona") 
# trtlevel <- c("Napier", "Napier + Desmodium", "Brachiaria cv. Mulato II", 
#               "Brachiaria + Desmodium",
#               "Desmodium distortum", "Maize", "Maize + Desmodium")
# levels(karama_df$treatment)

lasso_df$timepoint <- as.factor(lasso_df$timepoint)
lasso_df$treatment <- as.factor(lasso_df$treatment)
lasso_df$treatment<-factor(lasso_df$treatment,
                           levels=c("Napier grass", "Napier + Desmodium", 
                                    "Brachiaria cv. Mulato II", "Brachiaria + Desmodium",
                                    "Desmodium distortum", 
                                    "maize monoculture", "Maize + Desmodium"),
                           ordered=TRUE)

karama_df <- lasso_df %>% filter(location=="Karama") 
rubona_df <- lasso_df %>% filter(location=="Rubona") 
trtlevel <- c("Napier", "Napier + Desmodium", "Brachiaria cv. Mulato II", 
              "Brachiaria + Desmodium",
              "Desmodium distortum", "Maize", "Maize + Desmodium")
levels(karama_df$treatment)

gene_lmer_fun <- function(target_gene, input_df) {

  gene <- eval(substitute(target_gene), input_df)
  trtlevel <- c("Napier grass", "Napier + Desmodium", "Brachiaria cv. Mulato II", 
                            "Brachiaria + Desmodium",
                            "Desmodium distortum", "maize monoculture", "Maize + Desmodium")
  ##### fit mixed model ######
  model <- lmer(gene ~ treatment -1 + (1|timepoint) + (1|block), data=input_df)
  mod_df <- summary(model)
  
  mod_df1 <- as.data.frame(mod_df$coefficients) %>%
    rownames_to_column(var="treatment") %>%
    rename(estimate = Estimate,
           se = `Std. Error`,
           pval = `Pr(>|t|)`) %>% 
    mutate(treatment = str_split(treatment, "treatment", simplify=TRUE)[,2]) 

  #####create dfs for graphing#####
  napier <- mod_df1 %>% filter(treatment=="Napier grass")
  vline <- napier$estimate
  results <- mod_df1 %>% anti_join(napier) %>% as.data.frame()
  place <- cbind(trtlevel, c(7, 6, 2, 1, 3, 5, 4)) %>%
    as_tibble() %>%
    dplyr::rename(treatment = trtlevel,
           pos = `V2`) %>%
    mutate(treatment = as.character(treatment),
           pos = as.numeric(pos))

  #### extract sig contrasts ####
  t <- tidy(glht((model), linfct = mcp(treatment = "Dunnett"), alternative = "two.sided")) %>%
    as_tibble() %>%
    mutate(treatment = str_split_fixed(contrast, "-", 2)[,1],
           treatment = as.character(str_trim(treatment, side="right"))) %>%
    dplyr::select(-c(estimate)) %>%
    left_join(place, by="treatment") %>%
    left_join(results, by="treatment") %>%
    filter(as.numeric(adj.p.value) < 0.07) 
  
  vals <- c()
  trts <- c()
  est <- c()
  min <- c()
  
  for(i in 1:nrow(t)) {       # for-loop over rows
    vals[i] <- t$adj.p.value[i]
    trts[i] <- t$pos[i]
    est[i] <- t$estimate
    min[i] <- t$std.error
  }
  
  trts_pos <- trts 
    
  xpos <- est + min
  

  #####make custom forest plot#####
  title <- paste(substitute(target_gene))
 # pal2 <- c("#B00000", "#FFB4B3",  "#065F9B", "#C0E8FD", "#6A3D9A", "#006D00", "#CCF9A4")
  pal2 <- c("#C0E8FD", "#065F9B", "#6A3D9A", "#CCF9A4", "#006D00","#FFB4B3")
  
  
  forestplot <-
  ggplot(data=results, aes(y=treatment, x=estimate, xmin= estimate - se,
      xmax=estimate + se, color=treatment)) +
  geom_point(size=2) +
  geom_errorbarh(height=0, aes(color=treatment), size=2) +
  scale_color_manual(values = pal2) +
  geom_vline(xintercept=vline, color='black', linetype='dashed', alpha=.8) +
  theme_classic() +
  theme(legend.position="none", 
        axis.text.y = element_blank())+
       # legend.title = element_blank()) +
  ylab("") +
  xlab("Emmean") +
   geom_rect(data = napier, aes(xmin = estimate-se, xmax = estimate + se,
                                    ymin = -Inf, ymax = Inf),
             fill="#FFB3B3", alpha = 0.2, inherit.aes=FALSE) +
  #ggtitle(label=title) +
  theme(axis.ticks.y = element_blank()) +
  annotate("text", x=xpos, y=trts_pos+0.2, label= "*", size=7) +
  annotate("text", x=-Inf, y=Inf, label=title, size=3.5, hjust=0, vjust="top") 

  final <- list("summary" = results, "plot" = forestplot, "model" = model, "dunnetts" = t)
  return(final)
}

num <- gene_lmer_fun(Beta_amoA, karama_df)
# 
# 
# num$plot + theme(legend.position = "bottom")
# num$summary
# num$plot
# y <- num$summary
# conft <- confint(num$model)
# 
# sjPlot::tab_model(num$model)

#### testing ####
# test <- lmer(Arch_amoA_FA ~ treatment + (1|block) + (1|timepoint), data=karama_df)
# df <- summary(test)
# df <- as.data.frame(df$coefficients) %>%
#   rownames_to_column(var="treatment") %>%
#   rename(estimate = Estimate,
#          se = `Std. Error`,
#          pval = `Pr(>|t|)`) %>%
#   mutate(treatment = factor(treatment, labels=trtlevel)) %>%
#   filter(!grepl("Napier", treatment))


####ggplot graph custom forest #####
# result_nap <- df %>%
#   filter(treatment=="Napier")
# vline <- result_nap$estimate
# 
# plotest<-
#   ggplot(data=df,
#        aes(y=treatment, x=estimate, xmin= estimate - se,
#            xmax=estimate + se, color=treatment)) +
#   geom_point() +
#   geom_errorbarh(height=.1, aes(color=treatment)) +
#   #scale_y_continuous(breaks=1:nrow(df), labels=df$study) +
#   #labs(title='Effect Size by Study', x='Effect Size', y = 'Study') +
#   geom_vline(xintercept=vline, color='black', linetype='dashed', alpha=.5) +
#   theme_classic() +
#   theme(legend.position="none") +
#   ylab("") +
 #  geom_rect(data = result_nap, aes(xmin = estimate-se, xmax = estimate + se,
 #                                   ymin = -Inf, ymax = Inf), alpha = 0.1, inherit.aes=FALSE) +
 #  ggtitle("norB_2") +
 # # annotate("text", x=0.7, y=5, label="p=0.03") +
 # # annotate("text", x=0.73, y=6, label="p=0.06") +
 #  xlim(-0.7,0.75)
# 
# #   
# # emmeans(test, trt.vs.ctrl~treatment, ref=7, type="response")
# # #why does emmeans give no sig differences and lmer output does?
# ## answer: output indicates whether coefficient is sig diff from 0

