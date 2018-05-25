setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/Behavioural/Model Comparison')
data <- read.csv('RF2_trans_expData.csv',header=T)
mod1 <- read.csv('mod1_predict.csv', header=T)
mod2 <- read.csv('mod2_predict.csv', header=T)
# Exclude all rows with 'NA' 
# Recode Money_Amount as categorical
dataRaw$context <- as.factor(dataRaw$context)
# Recode Money_Amount as categorical
mod1$context <- ifelse(mod1$gainVal==0.5 & mod1$lossVal==0.5, "FF",
                       ifelse(mod1$gainVal==1.0 & mod1$lossVal==1.0, "OO",
                              ifelse(mod1$gainVal==0.5 & mod1$lossVal==1.0, "FO",
                                     ifelse(mod1$gainVal==1.0 & mod1$lossVal==0.5, "OF","NA"))))
mod1$context <- as.factor(mod1$context)
mod2$context <- ifelse(mod2$gainVal==0.5 & mod2$lossVal==0.5, "FF",
                       ifelse(mod2$gainVal==1.0 & mod2$lossVal==1.0, "OO",
                              ifelse(mod2$gainVal==0.5 & mod2$lossVal==1.0, "FO",
                                     ifelse(mod2$gainVal==1.0 & mod2$lossVal==0.5, "OF","NA"))))
mod2$context <- as.factor(mod2$context)

library("lme4")
library("effects")
library("ggplot2")
library("lmerTest")
source('../sortData.R')

# Clean and sort out data
# Data sorting function
dataClean <- sortData(data,removeNA=TRUE) 

# Only accurate trials: rt ~ context*stim
# Parse data
dataAcc <- dataClean[ which(dataClean$Accuracy_inTime==1), ]
mod1_Acc <- mod1[ which(mod1$Accuracy_inTime==1), ]
mod2_Acc <- mod2[ which(mod2$Accuracy_inTime==1), ]

# Compute glm
mod_accRT <- lmer(rt_sec ~ context*stim + (1 | subj_idx/contextNo), data=dataAcc)
anova(mod_accRT)


# Look at comp model predictions
compMod_1 <- lm(rt_sec ~ context*stim, data=mod1_Acc); anova(compMod_1)
compMod_2 <- lm(rt_sec ~ context*stim,  data=mod2_Acc); anova(compMod_2)



# Simple effect tests! 
# Recoding for simple effects (Aiken and West 1991)
dataAcc$lowContrast <- ifelse(dataAcc$context == "FF", 0, 1)
dataAcc$highContrast <- ifelse(dataAcc$context == "OO", 0,1)
dataAcc$prevContrast <- ifelse(dataAcc$context == "FO", 0,1)
dataAcc$promContrast <- ifelse(dataAcc$context == "OF", 0,1)
dataAcc$gainVnoloss <- ifelse(dataAcc$stim == "Gain  ", 1, ifelse(dataAcc$stim == "NoLoss", -1, 0))
modAcc_simple1 <- lmer(rt ~ lowContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataAcc); anova(modAcc_simple1)
modAcc_simple2 <- lmer(rt ~ highContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataAcc); anova(modAcc_simple2)
modAcc_simple3 <- lmer(rt ~ prevContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataAcc); summary(modAcc_simple3)
modAcc_simple4 <- lmer(rt ~ promContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataAcc); summary(modAcc_simple4)

# Computational model simple effects
mod2_Acc$lowContrast <- ifelse(mod2_Acc$context == "FF", 0, 1)
mod2_Acc$highContrast <- ifelse(mod2_Acc$context == "OO", 0,1)
mod2_Acc$prevContrast <- ifelse(mod2_Acc$context == "FO", 0,1)
mod2_Acc$promContrast <- ifelse(mod2_Acc$context == "OF", 0,1)
mod2_Acc$gainVnoloss <- ifelse(mod2_Acc$stim == "Gain  ", 1, ifelse(mod2_Acc$stim == "NoLoss", -1, 0))
mod2Acc_simple1 <- lm(rt ~ lowContrast*gainVnoloss, data=mod2_Acc); anova(mod2Acc_simple1)
mod2Acc_simple2 <- lm(rt ~ highContrast*gainVnoloss, data=mod2_Acc); anova(mod2Acc_simple2)
mod2Acc_simple3 <- lm(rt ~ prevContrast*gainVnoloss, data=mod2_Acc); summary(mod2Acc_simple3)
mod2Acc_simple4 <- lm(rt ~ promContrast*gainVnoloss, data=mod2_Acc); summary(mod2Acc_simple4)


# Plot effects
ef_accRT <- as.data.frame(effect('context:stim', mod_accRT,  confidence.level = 0.95))
ef_compMod_1 <- as.data.frame(effect('context:stim', compMod_1,  confidence.level = 0.95))
ef_compMod_2 <- as.data.frame(effect('context:stim', compMod_2,  confidence.level = 0.95))

ctxtLabels <- c(
  `FF`="LM",
  `OO`="HM",
  `FO`="Prev",
  `OF`="Prom")

ggplot(ef_accRT, aes(x = stim, y=fit, width=0.75)) + 
  facet_grid(~context, labeller=as_labeller(ctxtLabels), switch="both") + 
  geom_col(aes(fill=stim),color="black",size=1,width=0.6,show.legend = TRUE, position="dodge") +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8,position=position_dodge(width=0.75)) + 
  #geom_point(data=ef_compMod_1,aes(stim, fit),size=3, color="darkgreen",alpha=0.8) +
  #geom_line(data=ef_compMod_1,aes(group=context),linetype="dashed",size=0.8, color="darkgreen",alpha=0.8) +
  geom_point(data=ef_compMod_2,aes(stim, fit),size=3, color="red3",alpha=0.8) +
  geom_line(data=ef_compMod_2,aes(group=context),size=0.8, color="red3",alpha=0.8) +
  scale_fill_manual(labels = c("Control","Gain","No Loss"),values=c("gray20","gray70","gray40")) +
  ylab("Reaction Time (s)") +
  xlab("Context") + 
  labs(fill = "Cue") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = NA),
        strip.background = element_rect(colour="white", fill="white"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 14, family="Arial"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") +
  coord_cartesian(ylim = c(0.3, 0.35)) 
  



## For errors when NEUTRAL cue 
dataErr <- dataClean[ which(dataClean$Accuracy!=1), ]
# Parse data
dataErr_Neut <- dataErr[ which(dataErr$stim=="Cntrl "), ]
# Recode Response_Type so that No-loss response is 0 and gain is 1
dataErr_Neut$Response_Type <- relevel(dataErr_Neut$Response_Type, ref="NoLoss")

# Look at errors in computational models
mod1_Err <- mod1[ which(mod1$Accuracy!=1), ]
mod1_Err_Neut <- mod1_Err[ which(mod1_Err$stim=="Cntrl "), ]

mod2_Err <- mod2[ which(mod2$Accuracy!=1), ]
mod2_Err_Neut <- mod2_Err[ which(mod2_Err$stim=="Cntrl "), ]


# Run glms on errors
mod_errNeut_RT <- lmer(rt_sec ~ context*Response_Type + (1 | subj_idx/contextNo), data=dataErr_Neut)
anova(mod_errNeut_RT)


# Look at comp model predictions
compMod_1_err <- lm(rt_sec ~ context*Response_Type, data=mod1_Err_Neut); anova(compMod_1_err)
compMod_2_err <- lm(rt_sec ~ context*Response_Type, data=mod2_Err_Neut); anova(compMod_1_err)



# Simple effect tests! 
# Recoding for simple effects (Aiken and West 1991)
dataErr_Neut$lowContrast <- ifelse(dataErr_Neut$context == "FF", 0, 1)
dataErr_Neut$highContrast <- ifelse(dataErr_Neut$context == "OO", 0,1)
dataErr_Neut$prevContrast <- ifelse(dataErr_Neut$context == "FO", 0,1)
dataErr_Neut$promContrast <- ifelse(dataErr_Neut$context == "OF", 0,1)
dataErr_Neut$gainVnoloss <- ifelse(dataErr_Neut$Response_Type == "Gain  ", 1, ifelse(dataErr_Neut$stim == "Response_Type", -1, 0))
modAcc_simple3 <- lmer(rt ~ prevContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataErr_Neut); summary(modAcc_simple3)
modAcc_simple4 <- lmer(rt ~ promContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataErr_Neut); summary(modAcc_simple4)

# Plot effects
ef_errNeutRT <- as.data.frame(effect('context:Response_Type', mod_errNeut_RT,  confidence.level = 0.95))
ef_compMod_1_errNeutRT <- as.data.frame(effect('context:Response_Type', compMod_1_err, xlevels=2, confidence.level = 0.95))
ef_compMod_2_errNeutRT <- as.data.frame(effect('context:Response_Type', compMod_2_err, xlevels=2,  confidence.level = 0.95))


ggplot(ef_errNeutRT, aes(context, fit, group = Response_Type)) +
  geom_col(aes(fill=Response_Type),color="black",size=1,width=0.6,show.legend = TRUE, position="dodge") +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8,position=position_dodge(width=0.6)) + 
  geom_pointrange(data=ef_compMod_1_errNeutRT, aes(x=context, y=fit,  ymin=lower, ymax=upper), color='red', position=position_dodge(width=0.75), alpha=0.6) +
  geom_pointrange(data=ef_compMod_2_errNeutRT, aes(x=context, y=fit, ymin=lower, ymax=upper), color='dark green', position=position_dodge(width=0.75), alpha=0.6) + 
  labs(color = "Response") +
  scale_x_discrete(labels=c("LM", "Prev", "Prom", "HM")) + 
  scale_fill_grey(labels = c("Gain","No Loss"), start=0.8, end=0.4) +
  labs(fill="Response") + 
  ylab("Reaction Time (s)") +
  xlab("Context") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 14,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "top") + 
  coord_cartesian(ylim=c(0.325,0.375))





### Supplementary analyses

# Overall reaction time ~ context*stim, collapse across and accuracy
mod_genRT <- lmer(rt ~ context*stim + (1 | subj_idx/contextNo), data=dataClean)
anova(mod_genRT)
# Simple effect tests!
# Recoding for simple effects (Aiken and West 1991)
dataClean$lowContrast <- ifelse(dataClean$context == "FF", 0, 1)
dataClean$highContrast <- ifelse(dataClean$context == "OO", 0,1)
dataClean$prevContrast <- ifelse(dataClean$context == "FO", 0,1)
dataClean$promContrast <- ifelse(dataClean$context == "OF", 0,1)
dataClean$gainVnoloss <- ifelse(dataClean$stim == "Gain  ", 1, ifelse(dataClean$stim == "NoLoss", -1, 0))
mod_simple1 <- lmer(rt ~ lowContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataClean); anova(mod_simple1)
mod_simple2 <- lmer(rt ~ highContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataClean); anova(mod_simple2)
mod_simple3 <- lmer(rt ~ prevContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataClean); summary(mod_simple3)
mod_simple4 <- lmer(rt ~ promContrast*gainVnoloss + (1 | subj_idx/contextNo), data=dataClean); summary(mod_simple4)


# Plot effects (omnibus)
ef_genRT <- effect('context:stim', mod_genRT,  confidence.level = 0.95)
summary(ef_genRT)

ggplot(as.data.frame(ef_genRT), aes(context, fit, group = stim, width=0.75)) +
  geom_point(position=position_dodge(0.25)) +
  geom_pointrange(aes(ymin=lower, ymax=upper, color = stim), size = 1, position=position_dodge(0.25)) +
  #geom_jitter(data=dataClean,aes(x=context,y=rt,color=stim),position=position_dodge(0.25)) +
  scale_colour_manual(labels = c("Control","Gain","No Loss"), values = c("black", gray(0.35), gray(0.70))) +
  labs(color = "Cue") +
  scale_x_discrete(labels=c("Low Reward", "Prevention", "Promotion", "High Reward")) +
  ylab("Reaction Time (s)") +
  xlab("Context") +
  ggtitle("RT Across Hits and errors") +
  theme_bw()+ theme(axis.text.x = element_text(size = 12), legend.text = element_text(size = 14))



## For errors when GAIN cue 
# Parse data
dataErr_Gain <- dataErr[ which(dataErr$stim=="Gain  "), ]

mod_errGain_RT <- lmer(rt ~ context*Response_Type + (1 | subj_idx/contextNo), data=dataErr_Gain)
anova(mod_errGain_RT)

# Plot effects
ef_errGainRT <- effect('context:Response_Type', mod_errGain_RT,  confidence.level = 0.95)
summary(ef_errGainRT)

ggplot(as.data.frame(ef_errGainRT), aes(context, fit, group = Response_Type, width=0.25)) +
  scale_colour_manual(labels = c("Gain","No Loss"), values = c("#9EBCDA", "#4D004B")) +
  geom_crossbar(aes(ymin=lower, ymax=upper, color = Response_Type), size = 1, position=position_dodge(0.35)) + 
  geom_pointrange(aes(ymin=fit-se, ymax=fit+se, color = Response_Type), size = 1, position=position_dodge(0.35)) +
  labs(color = "Response") +
  scale_x_discrete(labels=c("Low Reward", "Prevention", "Promtion", "High")) + 
  ylab("Reaction Time (s)") +
  xlab("Context") +
  ggtitle("RT for Errors on Gain Cue Trials") +
  theme_bw()

## For errors when NO LOSS cue 
# Parse data
dataErr_NoLoss <- dataErr[ which(dataErr$stim=="NoLoss"), ]

mod_errNL_RT <- lmer(rt ~ context*Response_Type + (1 | subj_idx/contextNo), data=dataErr_NoLoss)
anova(mod_errNL_RT)

# Plot effects
ef_errNLRT <- effect('context:Response_Type', mod_errNL_RT,  confidence.level = 0.95)
summary(ef_errNLRT)

ggplot(as.data.frame(ef_errNLRT), aes(context, fit, group = Response_Type, width=0.25)) +
  scale_colour_manual(labels = c("Gain","No Loss"), values = c("#9EBCDA", "#4D004B")) +
  geom_crossbar(aes(ymin=lower, ymax=upper, color = Response_Type), size = 1, position=position_dodge(0.35)) + 
  geom_pointrange(aes(ymin=fit-se, ymax=fit+se, color = Response_Type), size = 1, position=position_dodge(0.35)) +
  labs(color = "Response") +
  scale_x_discrete(labels=c("Low Reward", "Prevention", "Promtion", "High")) + 
  ylab("Reaction Time (s)") +
  xlab("Context") +
  ggtitle("RT for Errors on No-Loss Cue Trials") +
  theme_bw()

########## Three way interaction between response type, stimulus cue, and context in RT ################ 
## For all trials (need both acc and err trials for 3-way interaction!)

mod_fullRT <- lmer(rt ~ context*stim*Response_Type + (1 | subj_idx/contextNo), data=dataClean)
anova(mod_fullRT)

# Plot effects
ef_fullRT <- effect('context:stim:Response_Type', mod_fullRT,  confidence.level = 0.95)
summary(ef_fullRT)

# Facet labels: 
contextID <- c(
  `FF` = "Low Reward",
  `FO` = "Prevention",
  `OF` = "Promotion",
  `OO` = "HighReward"
)

ggplot(as.data.frame(ef_fullRT), aes(stim, fit, group = Response_Type, width=0.25)) + 
  facet_wrap(~context,labeller = as_labeller(contextID)) + 
  scale_colour_manual(labels = c("Control","Gain","No Loss"), values = c("#9EBCDA", "#4D004B","dark blue")) +
  geom_crossbar(aes(ymin=lower, ymax=upper, color = Response_Type), size = 1, position=position_dodge(0.35)) + 
  geom_pointrange(aes(ymin=fit-se, ymax=fit+se, color = Response_Type), size = 1, position=position_dodge(0.35)) +
  labs(color = "Response") +
  scale_x_discrete(labels= c("Control","Gain","No Loss")) + 
  ylab("Reaction Time (s)") +
  xlab("Cue") +
  ggtitle("RT for Stimulus/Response Interactions") +
  theme_bw() + theme(axis.text.x = element_text(size = 14), legend.text = element_text(size = 14))
