setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex1/Results/')
data <- read.csv('RF1_trans_expData.csv',header=T)


# Load in different model data
mod1 <- read.csv('./Model Comparison/mod1_predict.csv', header=T)
mod2 <- read.csv('./Model Comparison/mod2_predict.csv', header=T)
mod3 <- read.csv('./Model Comparison/mod3_predict.csv', header=T)
mod4 <- read.csv('./Model Comparison/mod4_predict.csv', header=T)
mod5 <- read.csv('./Model Comparison/mod5_predict.csv', header=T)
mod6 <- read.csv('./Model Comparison/mod6_predict.csv', header=T)


# Exclude all rows with 'NA' 
data <- data[ which(data$rt!=0), ]
# Recode Money_Amount as categorical
data$Money_Amount <- as.factor(data$Money_Amount)
mod1$Money_Amount <- as.factor(mod1$Money_Amount)
mod2$Money_Amount <- as.factor(mod2$Money_Amount)
mod3$Money_Amount <- as.factor(mod3$Money_Amount)
mod4$Money_Amount <- as.factor(mod4$Money_Amount)
mod5$Money_Amount <- as.factor(mod5$Money_Amount)
mod6$Money_Amount <- as.factor(mod6$Money_Amount)


library("lme4")
library("effects")
library("ggplot2")
library("lmerTest")


## RT test:
# Only accurate trials: rt ~ context*stim
# Parse data
dataAcc <- data[ which(data$Accuracy_inTime==1), ]
mod1Acc <- mod1[ which(mod1$Accuracy_inTime==1), ]
mod2Acc <- mod2[ which(mod2$Accuracy_inTime==1), ]
mod3Acc <- mod3[ which(mod3$Accuracy_inTime==1), ]
mod4Acc <- mod4[ which(mod4$Accuracy_inTime==1), ]
mod5Acc <- mod5[ which(mod5$Accuracy_inTime==1), ]
mod6Acc <- mod6[ which(mod6$Accuracy_inTime==1), ]


# Run lmer on empirical data
dataRT <- lmer(rt_sec ~ stim*Money_Amount + (1 | subj_idx), data=dataAcc)
anova(dataRT)

# Run lm on model data 
mod1_RT <- lm(rt_sec ~ stim*Money_Amount , data=mod1Acc)
mod2_RT <- lm(rt_sec ~ stim*Money_Amount , data=mod2Acc)
mod3_RT <- lm(rt_sec ~ stim*Money_Amount , data=mod3Acc)
mod4_RT <- lm(rt_sec ~ stim*Money_Amount , data=mod4Acc)
mod5_RT <- lm(rt_sec ~ stim*Money_Amount , data=mod5Acc)
mod6_RT <- lm(rt_sec ~ stim*Money_Amount , data=mod6Acc)




# Plot effects
ef_dataRT <- as.data.frame(effect('stim:Money_Amount', dataRT,  confidence.level = 0.95))
# Plot comp modl effects
ef_mod1_RT <- as.data.frame(effect('stim:Money_Amount', mod1_RT, xlevels=2, confidence.level = 0.95))
ef_mod2_RT <- as.data.frame(effect('stim:Money_Amount', mod2_RT, xlevels=2, confidence.level = 0.95))
ef_mod3_RT <- as.data.frame(effect('stim:Money_Amount', mod3_RT, xlevels=2, confidence.level = 0.95))
ef_mod4_RT <- as.data.frame(effect('stim:Money_Amount', mod4_RT, xlevels=2, confidence.level = 0.95))
ef_mod5_RT <- as.data.frame(effect('stim:Money_Amount', mod5_RT, xlevels=2, confidence.level = 0.95))
ef_mod6_RT <- as.data.frame(effect('stim:Money_Amount', mod6_RT, xlevels=2, confidence.level = 0.95))

magLabels <- c(
  `0.5`="LM",
  `1`="HM")

ggplot(ef_dataRT, aes(stim, fit)) + 
  facet_grid(~Money_Amount, labeller=as_labeller(magLabels)) + 
  geom_pointrange(aes(ymin=lower, ymax=upper), color="black", size=0.8,show.legend=FALSE) +
  geom_line(aes(x=stim,group=Money_Amount), color="black", size=1) +
  geom_point(data=ef_mod2_RT,aes(stim, fit),size=2, color="orange",alpha=0.6) +
  geom_line(data=ef_mod2_RT,aes(group=Money_Amount),size=0.5, color="orange",alpha=0.6) +
  geom_point(data=ef_mod3_RT,aes(stim, fit),size=2, color="gold",alpha=0.6) +
  geom_line(data=ef_mod3_RT,aes(group=Money_Amount),size=0.5, color="gold",alpha=0.6) +
  geom_point(data=ef_mod4_RT,aes(stim, fit),size=2, color="dark green",alpha=0.6) +
  geom_line(data=ef_mod4_RT,aes(group=Money_Amount),size=0.5, color="dark green",alpha=0.6) +
  geom_point(data=ef_mod5_RT,aes(stim, fit),size=2, color="blue",alpha=0.6) +
  geom_line(data=ef_mod5_RT,aes(group=Money_Amount),size=0.5, color="blue",alpha=0.6) +
  geom_point(data=ef_mod6_RT,aes(stim, fit),size=2, color="purple",alpha=0.6) +
  geom_line(data=ef_mod6_RT,aes(group=Money_Amount),size=0.5, color="purple",alpha=0.6) +
  geom_point(data=ef_mod1_RT,aes(stim, fit),size=2,color="red",alpha=0.6) +
  geom_line(data=ef_mod1_RT,aes(group=Money_Amount),size=0.5, color="red",alpha=0.6) +
  scale_x_discrete(labels = c("Control","Gain","No Loss"))+
  ylab("Reaction Time (s)") +
  xlab("Cue") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 18, family="Arial"),
        axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 18,family="Arial"), 
        legend.text = element_text(size = 18,family="Arial"),
        axis.title = element_text(size = 20,family="Arial"))


## ACC test:
# Run lmer on empirical data
dataAccuracy <- glmer(Accuracy ~ stim * Money_Amount + (1 | subj_idx), family=binomial, data=data)

# Run lm on model data 
mod1_Accuracy <- glm(Accuracy ~ stim * Money_Amount, family=binomial, data=mod1)
mod2_Accuracy <- glm(Accuracy ~ stim * Money_Amount, family=binomial, data=mod2)
mod3_Accuracy <- glm(Accuracy ~ stim * Money_Amount, family=binomial, data=mod3)
mod4_Accuracy <- glm(Accuracy ~ stim * Money_Amount, family=binomial, data=mod4)
mod5_Accuracy <- glm(Accuracy ~ stim * Money_Amount, family=binomial, data=mod5)
mod6_Accuracy <- glm(Accuracy ~ stim * Money_Amount, family=binomial, data=mod6)

# Plot effects
ef_dataAcc <- as.data.frame(effect('stim:Money_Amount', dataAccuracy,  confidence.level = 0.95))
# Plot comp model effects
ef_mod1_Acc <- as.data.frame(effect('stim:Money_Amount', mod1_Accuracy, xlevels=2, confidence.level = 0.95))
ef_mod2_Acc <- as.data.frame(effect('stim:Money_Amount', mod2_Accuracy, xlevels=2, confidence.level = 0.95))
ef_mod3_Acc <- as.data.frame(effect('stim:Money_Amount', mod3_Accuracy, xlevels=2, confidence.level = 0.95))
ef_mod4_Acc <- as.data.frame(effect('stim:Money_Amount', mod4_Accuracy, xlevels=2, confidence.level = 0.95))
ef_mod5_Acc <- as.data.frame(effect('stim:Money_Amount', mod5_Accuracy, xlevels=2, confidence.level = 0.95))
ef_mod6_Acc <- as.data.frame(effect('stim:Money_Amount', mod6_Accuracy, xlevels=2, confidence.level = 0.95))

magLabels <- c(
  `0.5`="LM",
  `1`="HM")

ggplot(ef_dataAcc, aes(stim, fit)) + 
  facet_grid(~Money_Amount, labeller=as_labeller(magLabels)) + 
  geom_pointrange(aes(ymin=lower, ymax=upper), color="black", size=0.8,show.legend=FALSE) +
  geom_line(aes(x=stim,group=Money_Amount), color="black", size=1) +
  geom_point(data=ef_mod2_Acc,aes(stim, fit),size=2, color="orange",alpha=0.6) +
  geom_line(data=ef_mod2_Acc,aes(group=Money_Amount),size=0.5, color="orange",alpha=0.6) +
  geom_point(data=ef_mod3_Acc,aes(stim, fit),size=2, color="gold",alpha=0.6) +
  geom_line(data=ef_mod3_Acc,aes(group=Money_Amount),size=0.5, color="gold",alpha=0.6) +
  geom_point(data=ef_mod4_Acc,aes(stim, fit),size=2, color="dark green",alpha=0.6) +
  geom_line(data=ef_mod4_Acc,aes(group=Money_Amount),size=0.5, color="dark green",alpha=0.6) +
  geom_point(data=ef_mod5_Acc,aes(stim, fit),size=2, color="blue",alpha=0.6) +
  geom_line(data=ef_mod5_Acc,aes(group=Money_Amount),size=0.5, color="blue",alpha=0.6) +
  geom_point(data=ef_mod6_Acc,aes(stim, fit),size=2, color="purple",alpha=0.6) +
  geom_line(data=ef_mod6_Acc,aes(group=Money_Amount),size=0.5, color="purple",alpha=0.6) +
  geom_point(data=ef_mod1_Acc,aes(stim, fit),size=2,color="red",alpha=0.6) +
  geom_line(data=ef_mod1_Acc,aes(group=Money_Amount),size=0.5, color="red",alpha=0.6) +
  scale_x_discrete(labels = c("Control","Gain","No Loss"))+
  ylab("Accuracy") +
  xlab("Cue") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 18, family="Arial"),
        axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 18,family="Arial"), 
        legend.text = element_text(size = 18,family="Arial"),
        axis.title = element_text(size = 20,family="Arial"))


