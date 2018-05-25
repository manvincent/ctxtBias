setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex1/Results')
data <- read.csv('RFlex1_data.csv',header=T)
compModData <- read.csv('./Model Comparison/mod1_predict.csv', header=T)
# Exclude all rows with 'NA' 
data <- data[ which(data$rt!=0), ]
# Recode Money_Amount as categorical
data$Money_Amount <- as.factor(data$Money_Amount)
compModData$Money_Amount <- as.factor(compModData$Money_Amount)

library("lme4")
library("effects")
library("ggplot2")
library("lmerTest")

# Only accurate trials: rt ~ context*stim
# Parse data
dataAcc <- data[ which(data$Accuracy_inTime==1), ]
compDataAcc <- compModData[ which(compModData$Accuracy_inTime==1), ]

# Run lmer on empirical data
mod_accRT <- lmer(rt_sec ~ stim*Money_Amount + (1 | subj_idx), data=dataAcc)
anova(mod_accRT)

# Run lm on model data
compMod_accRT <- lm(rt_sec ~ stim*Money_Amount , data=compDataAcc)
anova(compMod_accRT)

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


# Plot effects
ef_accRT <- effect('stim:Money_Amount', mod_accRT,  confidence.level = 0.95)
# Plot comp modl effects
ef_compMod_accRT <- effect('stim:Money_Amount', compMod_accRT, xlevels=2, confidence.level = 0.95)

ggplot(as.data.frame(ef_accRT), aes(Money_Amount, fit, group = stim)) + 
  geom_col(aes(fill=stim),color="black",size=1,width=0.4,show.legend = TRUE, position="dodge") +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8,position=position_dodge(width=0.4)) + 
  geom_pointrange(data= as.data.frame(ef_compMod_accRT), aes(x=Money_Amount, y=fit, group=stim, ymin=lower, ymax=upper), color='red', position=position_dodge(width=0.4),alpha=0.6) +
  labs(fill = "Cue") +
  scale_fill_manual(labels = c("Control","Gain","No Loss"),values=c("gray20","gray70","gray40")) +
  scale_x_discrete(labels = c("LM","HM"))+
  ylab("Reaction Time (s)") +
  xlab("Outcome Condition") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  coord_cartesian(ylim = c(0.28, 0.35)) +
  theme(axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 14,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") 




# Look at errors -- when people don't respond accuracy, is there a bias in what they press? 
dataErr <- data[ which(data$Accuracy!=1), ]
## For errors when NEUTRAL cue 
# Parse data
dataErr_Neut <- dataErr[ which(dataErr$stim=="Cntrl "), ]
# Recode Response_Type so that No-loss response is 0 and gain is 1
dataErr_Neut$Response_Type <- relevel(dataErr_Neut$Response_Type, ref="NoLoss")

mod_respNeut <- lmer(rt_sec ~ Response_Type*Money_Amount + (1 | subj_idx),na.action=na.exclude, data=dataErr_Neut)
anova(mod_respNeut)

# Plot effects
ef_mod_respNeut <- effect('Response_Type:Money_Amount', mod_respNeut, confidence.level = 0.95)
summary(ef_mod_respNeut)

ggplot(as.data.frame(ef_mod_respNeut), aes(x=Money_Amount, y=fit, group=Response_Type)) + 
  geom_col(aes(fill=Response_Type),color="black",size=1,width=0.3,show.legend = TRUE,position="dodge") +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8,position=position_dodge(width=0.3)) + 
  scale_fill_grey(start=0.8, end=0.4) +
  scale_x_discrete(labels=c("LM","HM")) + 
  labs(fill="Response") + 
  ylab("Reaction Time (s)") +
  xlab("Outcome Condition") +
  coord_cartesian(ylim = c(0.3, 0.38)) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 14,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "top") 



