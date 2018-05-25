setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex1/Results')
data <- read.csv('RF1_trans_expData.csv',header=T)
compModData <- read.csv('./Model Comparison/mod1_predict.csv', header=T)
# Exclude all rows with 'NA' 
data <- data[ which(data$rt!=0), ]
# Recode Money_Amount as categorical
data$Money_Amount <- as.factor(data$Money_Amount)
compModData$Money_Amount <- as.factor(compModData$Money_Amount)

library("lme4")
library("lmerTest")
library("effects")
library("ggplot2")
library("car")


# Clean and sort out data

## All data - collapsing across context
# Looking at accuracy type (just pressing the correct shape, not necessarily in time - to orthogonalize from rt plot)
mod_accStim <- glmer(Accuracy ~ stim * Money_Amount + (1 | subj_idx), family=binomial, data=data)
Anova(mod_accStim)

# Look at comp model predictions
compMod_accStim <- glm(Accuracy ~ stim * Money_Amount, family=binomial, data=compModData)
Anova(compMod_accStim)
0
# Plot effects ## Plot not in paper
ef_mod_accStim <- effect('stim:Money_Amount', mod_accStim, xlevels=2, confidence.level = 0.95)
ef_compMod_accStim <- effect('stim:Money_Amount', compMod_accStim, xlevels=2, confidence.level = 0.95)

ggplot(as.data.frame(ef_mod_accStim), aes(x=Money_Amount, y=fit, group=stim)) +
  geom_col(aes(fill=stim),color="black",size=1,width=0.4,show.legend = TRUE, position="dodge") +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8,position=position_dodge(width=0.4)) +
  geom_pointrange(data= as.data.frame(ef_compMod_accStim), aes(x=Money_Amount, y=fit, group=stim, ymin=lower, ymax=upper), color='red', position=position_dodge(width=0.4), alpha=0.6) + 
  scale_fill_manual(labels = c("Control","Gain","No Loss"),values=c("gray20","gray70","gray40")) +
  scale_x_discrete(labels=c("LM","HM"))+
  labs(fill="Cue") + 
  ylab("Accuracy") +
  xlab("Outcome Condition") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  coord_cartesian(ylim = c(0.55, .95)) +
  theme(axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 14,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") 

### Raw data
library(plyr)
r1<-ddply(data, .(Money_Amount,stim), summarize, mean=mean(Accuracy))

se <- function(x) sqrt(var(x)/length(x))

r2<-ddply(data, .(Money_Amount,stim), summarize, se=se(Accuracy))
r1$se <- r2$se

r1$seMin <- r1$mean-r1$se
r1$seMax <- r1$mean+r1$se






# Look at errors -- when people don't respond accuracy, is there a bias in what they press? 
dataErr <- data[ which(data$Accuracy!=1), ]
## For errors when NEUTRAL cue 
# Parse data
dataErr_Neut <- dataErr[ which(dataErr$stim=="Cntrl "), ]
# Recode Response_Type so that No-loss response is 0 and gain is 1
dataErr_Neut$Response_Type <- relevel(dataErr_Neut$Response_Type, ref="NoLoss")


mod_respNeut <- glmer(Response_Type ~ Money_Amount + (1 | subj_idx), family=binomial, data=dataErr_Neut)
Anova(mod_respNeut)

# Plot effects
ef_mod_respNeut <- effect('Money_Amount', mod_respNeut, confidence.level = 0.95)
summary(ef_mod_respNeut)

ggplot(as.data.frame(ef_mod_respNeut), aes(x=Money_Amount, y=fit)) + 
  geom_col(aes(group=Money_Amount),fill="gray80",color="black",size=1,width=0.2,show.legend = FALSE) +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8) + 
  scale_fill_grey(start=0.3, end=0.7) +
  scale_x_discrete(labels=c("LM", "HM")) + 
  ylab("Proportion Response to Non-loss") +
  xlab("Outcome Condition") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 14,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") +
  coord_cartesian(ylim = c(0.3, 0.6)) 



