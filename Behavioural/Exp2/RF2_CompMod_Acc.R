setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/Behavioural/Model Comparison/')
dataRaw <- read.csv('RF2_trans_expData.csv',header=T)
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
library("lmerTest")
library("effects")
library("ggplot2")
library("car")
source('../sortData.R')

# Clean and sort out data
# Data sorting function
dataClean <- sortData(dataRaw,removeNA=TRUE) 

## ACC test:
dataClean <- dataClean[ which(dataClean$rt<=400), ]
mod1 <- mod1[ which(mod1$rt<=400), ]
mod2 <- mod2[ which(mod2$rt<=400), ]
mod3 <- mod3[ which(mod3$rt<=400), ]


## All data - collapsing across context
# Looking at accuracy type (just pressing the correct shape, not necessarily in time - to orthogonalize from rt plot)
mod_accStim <- glmer(Accuracy ~ stim + (1 | subj_idx/contextNo), family=binomial, data=dataClean)
Anova(mod_accStim)


## All data - testing at each context
# Looking at accuracy type (just pressing the correct shape, not necessarily in time - to orthogonalize from rt plot)
mod_accInter <- glmer(Accuracy ~ context*stim + (1 | subj_idx/contextNo), data=dataClean, family=binomial("logit"),control=glmerControl(optimizer = "bobyqa"),nAGQ = 1)
mod_accInter_nest <- glmer(Accuracy ~ stim + (1 | subj_idx/contextNo), data=dataClean, family=binomial("logit"),control=glmerControl(optimizer = "bobyqa"),nAGQ = 1)
anova(mod_accInter,mod_accInter_nest)
Anova(mod_accInter)


# Look at comp model predictions
compMod_1 <- glm(Accuracy ~ context*stim, family=binomial, data=mod1); Anova(compMod_1)
compMod_2 <- glm(Accuracy ~ context*stim, family=binomial, data=mod2); 
summary(compMod_2)

# Plot effects
ef_mod_accInter <- as.data.frame(effect('context:stim', mod_accInter, confidence.level = 0.95))
ef_compMod_1 <- as.data.frame(effect('context:stim', compMod_1,  confidence.level = 0.95))
ef_compMod_2 <- as.data.frame(effect('context:stim', compMod_2,  confidence.level = 0.95))


ggplot(ef_mod_accInter, aes(x = stim, y=fit, width=0.75)) + 
  facet_grid(~context, labeller=as_labeller(ctxtLabels), switch="both") + 
  geom_col(aes(fill=stim),color="black",size=1,width=0.6,show.legend = TRUE, position="dodge") +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8,position=position_dodge(width=0.75)) + 
  #geom_point(data=ef_compMod_1,aes(stim, fit), size=3, color="darkgreen",alpha=0.8) +
  #geom_line(data=ef_compMod_1,aes(group=context),linetype="dashed",size=0.8, color="darkgreen",alpha=0.8) +
  geom_point(data=ef_compMod_2,aes(stim, fit),size=3, color="red3",alpha=0.8) +
  geom_line(data=ef_compMod_2,aes(group=context),size=0.8, color="red3",alpha=0.8) +
  scale_fill_manual(labels = c("Control","Gain","No Loss"),values=c("gray20","gray70","gray40")) +
  ylab("Accuracy") +
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
  coord_cartesian(ylim = c(0.57, 1)) 


# Look at errors -- when people don't respond accurately, is there a bias in what they press? 
dataErr <- dataClean[ which(dataClean$Accuracy!=1), ]
## For errors when NEUTRAL cue 
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
mod_respNeut <- glmer(Response_Type ~ context + (1 | subj_idx/contextNo), family=binomial, data=dataErr_Neut)
mod_respNeut_nest <- glmer(Response_Type ~  (1 | subj_idx/contextNo), family=binomial, data=dataErr_Neut)
anova(mod_respNeut,mod_respNeut_nest)


# Look at comp model predictions
compMod_1_err <- glm(Response_Type ~ context, family=binomial, data=mod1_Err_Neut); Anova(compMod_1_err)
compMod_2_err <- glm(Response_Type ~ context, family=binomial, data=mod2_Err_Neut); Anova(compMod_1_err)

# Plot effects
ef_mod_respNeut <- as.data.frame(effect('context', mod_respNeut, confidence.level = 0.95))
ef_compMod_1_respNeut <- as.data.frame(effect('context', compMod_1_err,  confidence.level = 0.95))
ef_compMod_2_respNeut <- as.data.frame(effect('context', compMod_2_err,  confidence.level = 0.95))


ggplot(ef_mod_respNeut, aes(x=context, y=fit)) + 
  geom_col(aes(group=context),fill="gray80",color="black",size=1,width=0.4,show.legend = FALSE) +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8) + 
  geom_pointrange(data=ef_compMod_1_respNeut, aes(x=context, y=fit, ymin=lower, ymax=upper), color='red', position=position_dodge(width=0.75), alpha=0.6) +
  geom_pointrange(data=ef_compMod_2_respNeut, aes(x=context, y=fit,ymin=lower, ymax=upper), color='dark green', position=position_dodge(width=0.75), alpha=0.6) + 
  scale_fill_grey(start=0.3, end=0.7) +
  scale_x_discrete(labels=c("LM", "Prev", "Prom", "HM")) + 
  ylab("Proportion Response to Non-loss") +
  xlab("Context") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 14,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") 
  #coord_cartesian(ylim = c(0.3, 0.6)) 


