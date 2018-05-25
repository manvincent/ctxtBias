setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/EEG/TF/TF_Interac')

library("reshape2")
library("lme4")
library("car")
library("effects")
library("ggplot2")
library("lmerTest")
source('sortData.R')

### Import data for fixation and stim parts of trial
dataTrial <- read.csv('allSub_modEstimates_Trial.csv', header=T)

# Clean 'stim' data
# Re-convert stim column
dataTrial$stim <- as.factor(dataTrial$stim)
# Create contextNo column
dataTrial <- sortData(dataTrial,removeNA=TRUE) 
# Only look at odd trials
dataTrial_split <- dataTrial[ which(dataTrial$Data_Split == 2), ]

##### Trial TF analysis #######
# Model 1:
mod_EEG_rt <- lmer(rt ~  Ave_TF*context*stim + (1 | subj_idx), data=dataTrial)
anova(mod_EEG_rt)

# Plot effects
ef_EEG_rt <- as.data.frame(effect('Ave_TF:context:stim', mod_EEG_rt, xlevels=100, confidence.level = 0.95))
# Facet labels: 
contextID <- c(
  `FF` = "LM",
  `FO` = "Prev",
  `OF` = "Prom",
  `OO` = "HM"
)
ggplot(ef_EEG_rt, aes(x=Ave_TF,  fit)) + 
  facet_grid(~context,labeller = as_labeller(contextID))+
  geom_line(aes(color=stim),size=1, show.legend = FALSE) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=stim), alpha=0.7) + 
  scale_fill_manual(labels = c("Control","Gain","No Loss"), values=c("gray10","gray75","gray40")) +
  scale_color_manual(labels = c("Control","Gain","No Loss"), values=c("gray10","gray75","gray40")) +
  labs(fill = "Cue") +
  ylab("Reaction Time (s)") +
  xlab(expression('Low '~gamma~'Power')) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.text = element_text(size = 12,family="Arial"), 
        axis.text.x = element_text(size = 12,family="Arial"), 
        axis.text.y = element_text(size = 10,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") 



# Model 1:
mod_EEG_compMod1 <- lmer(Ave_TF ~ meanAcc + rt + (1 | subj_idx), data=dataTrial)
anova(mod_EEG_compMod1)

# Model 2:
mod_EEG_compMod2 <- lmer(Ave_TF ~  probAcc_norm + rt + (1 | subj_idx), data=dataTrial)
anova(mod_EEG_compMod2)
summary(mod_EEG_compMod2)

# Plot effects
ef_EEG_compMod2 <- as.data.frame(effect('probAcc_norm', mod_EEG_compMod2, xlevels=100, confidence.level = 0.95))

ggplot(ef_EEG_compMod2, aes(x=probAcc_norm,  fit)) + 
  geom_line(color="dark red",size=2) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="dark red", alpha= 0.5,size=1, show.legend=TRUE) + 
  xlab("Choice Accumulator Activation") +
  ylab(expression('Low '~gamma~'Power'))+
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.text = element_text(size = 12,family="Arial"), 
        axis.text.x = element_text(size = 12,family="Arial"), 
        axis.text.y = element_text(size = 10,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") 
 
