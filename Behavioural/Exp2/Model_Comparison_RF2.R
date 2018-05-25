setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/Behavioural/Model Comparison/')
data <- read.csv('RF2_trans_expData.csv',header=T)


# Load in different model data
mod1 <- read.csv('mod1_predict.csv', header=T)
mod2 <- read.csv('mod2_predict.csv', header=T)
mod3 <- read.csv('mod3_predict.csv', header=T)
# mod4 <- read.csv('./Model Comparison/mod4_predict.csv', header=T)
# mod5 <- read.csv('./Model Comparison/mod5_predict.csv', header=T)
# mod6 <- read.csv('./Model Comparison/mod6_predict.csv', header=T)


# Exclude all rows with 'NA' 
data <- data[ which(data$rt!=0), ]
# Recode Money_Amount as categorical
mod1$context <- ifelse(mod1$gainVal==0.5 & mod1$lossVal==0.5, "FF",
                       ifelse(mod1$gainVal==1.0 & mod1$lossVal==1.0, "OO",
                       ifelse(mod1$gainVal==0.5 & mod1$lossVal==1.0, "FO",
                       ifelse(mod1$gainVal==1.0 & mod1$lossVal==0.5, "OF","NA"))))
mod2$context <- ifelse(mod2$gainVal==0.5 & mod2$lossVal==0.5, "FF",
                            ifelse(mod2$gainVal==1.0 & mod2$lossVal==1.0, "OO",
                                   ifelse(mod2$gainVal==0.5 & mod2$lossVal==1.0, "FO",
                                          ifelse(mod2$gainVal==1.0 & mod2$lossVal==0.5, "OF","NA"))))
mod3$context <- ifelse(mod3$gainVal==0.5 & mod3$lossVal==0.5, "FF",
                       ifelse(mod3$gainVal==1.0 & mod3$lossVal==1.0, "OO",
                              ifelse(mod3$gainVal==0.5 & mod3$lossVal==1.0, "FO",
                                     ifelse(mod3$gainVal==1.0 & mod3$lossVal==0.5, "OF","NA"))))

mod1$context <- as.factor(mod1$context)
mod2$context <- as.factor(mod2$context)
mod3$context <- as.factor(mod3$context)

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
# mod4Acc <- mod4[ which(mod4$Accuracy_inTime==1), ]
# mod5Acc <- mod5[ which(mod5$Accuracy_inTime==1), ]
# mod6Acc <- mod6[ which(mod6$Accuracy_inTime==1), ]


# Run lmer on empirical data
dataRT <- lmer(rt_sec ~ stim*context + (1 | subj_idx), data=dataAcc)
anova(dataRT)

# Run lm on model data 
mod1_RT <- lm(rt_sec ~ stim*context, data=mod1Acc); anova(mod1_RT)
mod2_RT <- lm(rt_sec ~ stim*context , data=mod2Acc);anova(mod2_RT)
mod3_RT <- lm(rt_sec ~ stim*context , data=mod3Acc)
# mod4_RT <- lm(rt_sec ~ stim*context , data=mod4Acc)
# mod5_RT <- lm(rt_sec ~ stim*context , data=mod5Acc)
# mod6_RT <- lm(rt_sec ~ stim*context , data=mod6Acc)




# Plot effects
ef_dataRT <- as.data.frame(effect('stim:context', dataRT,  confidence.level = 0.95))
# Plot comp modl effects
ef_mod1_RT <- as.data.frame(effect('stim:context', mod1_RT, xlevels=2, confidence.level = 0.95))
ef_mod2_RT <- as.data.frame(effect('stim:context', mod2_RT, xlevels=2, confidence.level = 0.95))
ef_mod3_RT <- as.data.frame(effect('stim:context', mod3_RT, xlevels=2, confidence.level = 0.95))
# ef_mod4_RT <- as.data.frame(effect('stim:context', mod4_RT, xlevels=2, confidence.level = 0.95))
# ef_mod5_RT <- as.data.frame(effect('stim:context', mod5_RT, xlevels=2, confidence.level = 0.95))
# ef_mod6_RT <- as.data.frame(effect('stim:context', mod6_RT, xlevels=2, confidence.level = 0.95))

ctxtLabels <- c(
  `FF`="LM",
  `OO`="HM",
  `FO`="Prev",
  `OF`="Prom")

ggplot(ef_dataRT, aes(stim, fit)) + 
  facet_grid(~context, labeller=as_labeller(ctxtLabels)) + 
  geom_pointrange(aes(ymin=lower, ymax=upper), color="black", size=0.8,show.legend=FALSE) +
  geom_line(aes(x=stim,group=context), color="black", size=1) +
  geom_point(data=ef_mod1_RT,aes(stim, fit),size=2, color="dark green",alpha=0.6) +
  geom_line(data=ef_mod1_RT,aes(group=context),size=0.5, color="dark green",alpha=0.6) +
  geom_point(data=ef_mod3_RT,aes(stim, fit),size=2, color="dark blue",alpha=0.6) +
  geom_line(data=ef_mod3_RT,aes(group=context),size=0.5, color="dark blue",alpha=0.6) +
  geom_point(data=ef_mod2_RT,aes(stim, fit),size=2,color="red",alpha=0.6) +
  geom_line(data=ef_mod2_RT,aes(group=context),size=0.5, color="red",alpha=0.6) +
  scale_x_discrete(labels = c("Control","Gain","No Loss"))+
  ylab("Reaction Time (s)") +
  xlab("Cue") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        #strip.text.x = element_text(size = 18, family="Arial"),
        strip.text.x=element_blank(),
        axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 18,family="Arial"), 
        legend.text = element_text(size = 18,family="Arial"),
        axis.title = element_text(size = 20,family="Arial"))
  


## ACC test:
data <- data[ which(data$rt<=400), ]
mod1 <- mod1[ which(mod1$rt<=400), ]
mod2 <- mod2[ which(mod2$rt<=400), ]
mod3 <- mod3[ which(mod3$rt<=400), ]

# Run lmer on empirical data
#dataAccuracy <- glmer(Accuracy ~ stim * context + (1 | subj_idx), family=binomial, data=data)

# Run lm on model data 
mod1_Accuracy <- glm(Accuracy ~ stim * context, family=binomial, data=mod1)
mod2_Accuracy <- glm(Accuracy ~ stim * context, family=binomial, data=mod2)
mod3_Accuracy <- glm(Accuracy ~ stim * context, family=binomial, data=mod3)
# mod4_Accuracy <- glm(Accuracy ~ stim * context, family=binomial, data=mod4)
# mod5_Accuracy <- glm(Accuracy ~ stim * context, family=binomial, data=mod5)
# mod6_Accuracy <- glm(Accuracy ~ stim * context, family=binomial, data=mod6)

# Plot effects
ef_dataAcc <- as.data.frame(effect('stim:context', dataAccuracy,  confidence.level = 0.95))
# Plot comp model effects
ef_mod1_Acc <- as.data.frame(effect('stim:context', mod1_Accuracy, xlevels=2, confidence.level = 0.95))
ef_mod2_Acc <- as.data.frame(effect('stim:context', mod2_Accuracy, xlevels=2, confidence.level = 0.95))
ef_mod3_Acc <- as.data.frame(effect('stim:context', mod3_Accuracy, xlevels=2, confidence.level = 0.95))
# ef_mod4_Acc <- as.data.frame(effect('stim:context', mod4_Accuracy, xlevels=2, confidence.level = 0.95))
# ef_mod5_Acc <- as.data.frame(effect('stim:context', mod5_Accuracy, xlevels=2, confidence.level = 0.95))
# ef_mod6_Acc <- as.data.frame(effect('stim:context', mod6_Accuracy, xlevels=2, confidence.level = 0.95))


ctxtLabels <- c(
  `FF`="LM",
  `OO`="HM",
  `FO`="Prev",
  `OF`="Prom")

ggplot(ef_dataAcc, aes(stim, fit)) + 
  facet_grid(~context, labeller=as_labeller(ctxtLabels)) + 
  geom_pointrange(aes(ymin=lower, ymax=upper), color="black", size=0.8,show.legend=FALSE) +
  geom_line(aes(x=stim,group=context), color="black", size=1) +
  geom_point(data=ef_mod1_Acc,aes(stim, fit),size=2, color="dark green",alpha=0.6) +
  geom_line(data=ef_mod1_Acc,aes(group=context),size=0.5, color="dark green",alpha=0.6) +
  geom_point(data=ef_mod3_Acc,aes(stim, fit),size=2, color="dark blue",alpha=0.6) +
  geom_line(data=ef_mod3_Acc,aes(group=context),size=0.5, color="dark blue",alpha=0.6) +
  geom_point(data=ef_mod2_Acc,aes(stim, fit),size=2,color="red",alpha=0.6) +
  geom_line(data=ef_mod2_Acc,aes(group=context),size=0.5, color="red",alpha=0.6) +
  scale_x_discrete(labels = c("Control","Gain","No Loss"))+
  ylab("Accuracy") +
  xlab("Cue") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 18, family="Arial"),
        #axis.text.x = element_text(size = 14,family="Arial"), 
        axis.text.y = element_text(size = 18,family="Arial"), 
        legend.text = element_text(size = 18,family="Arial"),
        axis.title = element_text(size = 20,family="Arial"))+ 
  coord_cartesian(ylim = c(0.55, 0.97)) 

