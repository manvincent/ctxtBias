setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/Behavioural')
dataRaw <- read.csv('RFS2_MLM.csv',header=T)
# Exclude all rows with 'NA' 
# Recode Money_Amount as categorical
dataRaw$context <- as.factor(dataRaw$context)

library("lme4")
library("lmerTest")
library("effects")
library("ggplot2")
library("car")
source('sortData.R')

# Clean and sort out data
# Data sorting function
dataClean <- sortData(dataRaw,removeNA=TRUE) 


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


# Plot effects
ef_mod_accInter <- as.data.frame(effect('context:stim', mod_accInter, confidence.level = 0.95))

ggplot(ef_mod_accInter, aes(x = context, y=fit, group=stim, width=0.75)) + 
  geom_col(aes(fill=stim),color="black",size=1,width=0.6,show.legend = TRUE, position="dodge") +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8,position=position_dodge(width=0.75)) + 
  scale_fill_manual(labels = c("Control","Gain","No Loss"),values=c("navyblue","lightskyblue1","steelblue")) +
  ylab("Accuracy") +
  xlab("Context") + 
  labs(fill = "Cue") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") + 
  coord_cartesian(ylim = c(0.65, 0.95)) 


# Look at errors -- when people don't respond accurately, is there a bias in what they press? 
dataErr <- dataClean[ which(dataClean$Accuracy!=1), ]
## For errors when NEUTRAL cue 
# Parse data
dataErr_Neut <- dataErr[ which(dataErr$stim=="Cntrl "), ]
# Recode Response_Type so that No-loss response is 0 and gain is 1
dataErr_Neut$Response_Type <- relevel(dataErr_Neut$Response_Type, ref="NoLoss")


# Run glms on errors
mod_respNeut <- glmer(Response_Type ~ context + (1 | subj_idx/contextNo), family=binomial, data=dataErr_Neut)
mod_respNeut_nest <- glmer(Response_Type ~  (1 | subj_idx/contextNo), family=binomial, data=dataErr_Neut)
anova(mod_respNeut,mod_respNeut_nest)


# Plot effects
ef_mod_respNeut <- as.data.frame(effect('context', mod_respNeut, confidence.level = 0.95))

ggplot(ef_mod_respNeut, aes(x=context, y=fit)) + 
  geom_col(aes(group=context),fill="navyblue",color="black",size=1,width=0.4,show.legend = FALSE, alpha=0.8) +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8) + 
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
        legend.position = "bottom") +
  coord_cartesian(ylim = c(0.3, 0.55)) 


