setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/EEG/FFT')

library("reshape2")
library("lme4")
library("effects")
library("ggplot2")
library("lmerTest")
source('sortData.R')

### Import data for fixation and stim parts of trial
dataTrialEEG_Stim <- read.csv('Group_flowGamma_FFTChan.csv', header=T)

# Clean 'stim' data
# Re-convert stim column
dataTrialEEG_Stim$stim <- ifelse(dataTrialEEG_Stim$stim == 1, "Gain", ifelse(dataTrialEEG_Stim$stim == -1, "NoLoss", "Cntrl"))
dataTrialEEG_Stim$stim <- as.factor(dataTrialEEG_Stim$stim)
# Re-convert response column
dataTrialEEG_Stim$response <- ifelse(dataTrialEEG_Stim$response == 1, "Gain", ifelse(dataTrialEEG_Stim$response == -1, "NoLoss", "Cntrl"))
dataTrialEEG_Stim$response <- as.factor(dataTrialEEG_Stim$response)
dataTrialEEG_Stim <- droplevels(dataTrialEEG_Stim)
# Create contextNo column
dataTrialEEG_Stim <- sortData(dataTrialEEG_Stim,removeNA=TRUE) 


##### Whole-trial analysis #######

# Step 1: Look at average EEG within block, after removing all non-loss trials (the stim of interest)
dataParsedEEG_Stim <- dataTrialEEG_Stim[ which(dataTrialEEG_Stim$stim != "NoLoss"),]
dataParsedEEG_Stim <- droplevels(dataParsedEEG_Stim)

# Get average EEG signal on stimulus for each block
dataParsedEEG_Stim <- data.frame(tapply(dataTrialEEG_Stim$TF_Stim,list(dataTrialEEG_Stim$subj_idx,dataTrialEEG_Stim$context),mean))
dataParsedEEG_Stim$subj_idx <- rownames(dataParsedEEG_Stim)
dataParsedEEG_Stim <- melt(dataParsedEEG_Stim,id.vars="subj_idx")
names(dataParsedEEG_Stim) <- c("subj_idx","context","BlockAve_TF_Stim")
  # Group-mean block EEG signal
dataParsedEEG_Stim$subjMean <-tapply(dataParsedEEG_Stim$BlockAve_TF,list(dataParsedEEG_Stim$subj_idx),mean)
dataParsedEEG_Stim$BlockAve_TF_mc <- dataParsedEEG_Stim$BlockAve_TF - dataParsedEEG_Stim$subjMean

# Get average EEG signal on fixation periods for each block
dataParsedEEG_Fix <- data.frame(tapply(dataTrialEEG_Stim$TF_Fix,list(dataTrialEEG_Stim$subj_idx,dataTrialEEG_Stim$context),mean))
dataParsedEEG_Fix$subj_idx <- rownames(dataParsedEEG_Fix)
dataParsedEEG_Fix <- melt(dataParsedEEG_Fix,id.vars="subj_idx")
names(dataParsedEEG_Fix) <- c("subj_idx","context","BlockAve_TF_Fix")
  # Group-mean block EEG signal
dataParsedEEG_Fix$subjMean <-tapply(dataParsedEEG_Fix$BlockAve_TF,list(dataParsedEEG_Fix$subj_idx),mean)
dataParsedEEG_Fix$BlockAve_TF_mc <- dataParsedEEG_Fix$BlockAve_TF - dataParsedEEG_Fix$subjMean

# Merge the fixation and stimulus data
dataParsedEEG_Stim$BlockAve_TF_Fix_mc <- dataParsedEEG_Fix$BlockAve_TF_mc
# Integrate this with behavioural data
dataEEG_Stim <- merge(dataTrialEEG_Stim, dataParsedEEG_Stim,, by=c("subj_idx","context"))
# Get only the accurate trials
dataAccEEG_Stim <- dataEEG_Stim[ which(dataEEG_Stim$Accuracy_inTime==1), ]


# Get mean RT per context x condition, per subject 
dataMeanCtxtRT_Stim <- aggregate(dataAccEEG_Stim$rt, by = list(dataAccEEG_Stim$subj_idx,dataAccEEG_Stim$context,dataAccEEG_Stim$stim), FUN=mean )
names(dataMeanCtxtRT_Stim) <- c("subj_idx","context","stim","meanRT")
dataStim <- merge(dataMeanCtxtRT_Stim, dataAccEEG_Stim)

# Include only non-loss trials
dataStim <- dataStim[ which(dataStim$stim=="NoLoss"), ]



#dataStim$BlockAve_TF_mc[dataStim$BlockAve_TF_mc > mean(dataStim$BlockAve_TF_mc) + 3*sd(dataStim$BlockAve_TF_mc)] = NA
#dataStim$BlockAve_TF_mc[dataStim$BlockAve_TF_mc > mean(dataStim$BlockAve_TF_mc) - 3*sd(dataStim$BlockAve_TF_mc)] = NA


# Model 
mod_RT_blockEEG_Stim <- lmer(meanRT ~ context*BlockAve_TF_mc*BlockAve_TF_Fix_mc + (1 | subj_idx), data=dataStim)
anova(mod_RT_blockEEG_Stim)

# Simple effects
dataStim$noLoss <- ifelse(dataStim$stim=="NoLoss",0,1)
dataStim$LMCtxt <- ifelse(dataStim$context=="FF",0,1)
dataStim$prevCtxt <- ifelse(dataStim$context=="FO",0,1)
dataStim$promCtxt <- ifelse(dataStim$context=="OF",0,1)
dataStim$HMCtxt <- ifelse(dataStim$context=="OO",0,1)

modAcc_simple_FFTtrial2a <- lmer(meanRT ~ prevCtxt*BlockAve_TF_mc + (1 | subj_idx), data=dataStim); summary(modAcc_simple_FFTtrial2a)
# effect size
tVal <- 8.597
tdf <-  7.981e+03
d <- (2*tVal/sqrt(tdf)); d

FVal <- tVal^2          
dfN <- 1
dfD <- tdf
eta <- ((dfN/dfD) * FVal)/(1+((dfN/dfD) * FVal)); eta

modAcc_simple_FFTtrial2b <- lmer(meanRT ~ promCtxt*BlockAve_TF_mc + (1 | subj_idx), data=dataStim); summary(modAcc_simple_FFTtrial2b)
modAcc_simple_FFTtrial2c <- lmer(meanRT ~ LMCtxt*BlockAve_TF_mc + (1 | subj_idx), data=dataStim); summary(modAcc_simple_FFTtrial2c)
modAcc_simple_FFTtrial2d <- lmer(meanRT ~ HMCtxt*BlockAve_TF_mc + (1 | subj_idx), data=dataStim); summary(modAcc_simple_FFTtrial2d)

# Plot effects
ef_RT_trialEEG_Stim <- effect('context:BlockAve_TF_mc ', xlevels=3, mod_RT_blockEEG_Stim, confidence.level = 0.95)
# Facet labels: 
contextID <- c(
  `FF` = "Low Magnitude",
  `FO` = "Prevention",
  `OF` = "Promotion",
  `OO` = "High Magnitude"
)

ggplot(as.data.frame(ef_RT_trialEEG_Stim), aes(x=BlockAve_TF_mc , fit)) + 
  geom_line(aes(color=context),show.legend=TRUE,size=2) +
  geom_linerange(aes(ymin=lower, ymax=upper,colour = context), alpha= 0.5,size=1, show.legend=TRUE) +
  scale_color_manual(labels = c("LM","Prev","Prom","HM"), values = c("dodgerblue4","darkgreen","dodgerblue1","dodgerblue3")) +
  #geom_point(data=dataStim, aes(x=BlockAve_TF_mc, y = meanRT, color=context) ) + 
  #scale_linetype_manual(breaks = c("LM","Prev","Prom","HM"), values = c(2, 1, 4, 5)) + 
  labs(color = "Context") +
  ylab("Mean Reaction Time (s)") +
  xlab(expression('Mean Block'~Delta~'Power')) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.text = element_text(size = 12,family="Arial"), 
        axis.text.x = element_text(size = 12,family="Arial"), 
        axis.text.y = element_text(size = 10,family="Arial"), 
        legend.text = element_text(size = 14,family="Arial"),
        axis.title = element_text(size = 18,family="Arial"),
        legend.title = element_text(size = 18,family="Arial"),
        legend.position = "bottom") 




 



 
# ####### Fixation analysis #####
# # Get average EEG signal for each block
# dataBlockMean_Fix <- data.frame(tapply(dataTrialEEG_Fix$Ave_TF,list(dataTrialEEG_Fix$subj_idx,dataTrialEEG_Fix$context),mean))
# dataBlockMean_Fix$subj_idx <- rownames(dataBlockMean_Fix)
# library("reshape2")
# dataBlockMean_Fix <- melt(dataBlockMean_Fix,id.vars="subj_idx")
# names(dataBlockMean_Fix) <- c("subj_idx","context","BlockAve_TF")
# # Integrate block-average EEG dataframe with trialwise dataframe
# dataEEG_Fix <- merge(dataTrialEEG_Fix, dataBlockMean_Fix)
# # Get only the accurate trials
# dataAccEEG_Fix <- dataEEG_Fix[ which(dataEEG_Fix$Accuracy_inTime==1), ]
# 
# # Get mean RT per context x condition, per subject 
# dataMeanCtxtRT_Fix <- aggregate(dataAccEEG_Fix$rt, by = list(dataAccEEG_Fix$subj_idx,dataAccEEG_Fix$context,dataAccEEG_Fix$stim), FUN=mean )
# names(dataMeanCtxtRT_Fix) <- c("subj_idx","context","stim","meanRT")
# dataFix <- merge(dataMeanCtxtRT_Fix, dataAccEEG_Fix)
# 
# 
# # Model 
# mod_RT_blockEEG_Fix <- lmer(meanRT ~ context*stim*BlockAve_TF + (1 | subj_idx/contextNo), data=dataFix)
# anova(mod_RT_blockEEG_Fix)
# 
# # Simple effects
# dataFix$noLoss <- ifelse(dataFix$stim=="NoLoss",0,1)
# dataFix$LMCtxt <- ifelse(dataFix$context=="FF",0,1)
# dataFix$prevCtxt <- ifelse(dataFix$context=="FO",0,1)
# dataFix$promCtxt <- ifelse(dataFix$context=="OF",0,1)
# dataFix$HMCtxt <- ifelse(dataFix$context=="OO",0,1)
# 
# modAcc_simple_FFTtrial1a <- lmer(meanRT ~ prevCtxt*noLoss*BlockAve_TF + (1 | subj_idx/contextNo), data=dataFix); anova(modAcc_simple_FFTtrial1a)
# modAcc_simple_FFTtrial1b <- lmer(meanRT ~ promCtxt*noLoss*BlockAve_TF + (1 | subj_idx/contextNo), data=dataFix); anova(modAcc_simple_FFTtrial1b)
# modAcc_simple_FFTtrial1c <- lmer(meanRT ~ LMCtxt*noLoss*BlockAve_TF + (1 | subj_idx/contextNo), data=dataFix); anova(modAcc_simple_FFTtrial1c)
# modAcc_simple_FFTtrial1d <- lmer(meanRT ~ HMCtxt*noLoss*BlockAve_TF + (1 | subj_idx/contextNo), data=dataFix); anova(modAcc_simple_FFTtrial1d)
# 
# 
# 
# # Plot effects
# ef_RT_trialEEG_fix <- effect('context:stim:BlockAve_TF', xlevels=2, mod_RT_blockEEG_Fix, confidence.level = 0.95)
# # Facet labels: 
# contextID <- c(
#   `FF` = "Low Magnitude",
#   `FO` = "Prevention",
#   `OF` = "Promotion",
#   `OO` = "High Magnitude"
# )
# 
# ggplot(as.data.frame(ef_RT_trialEEG_fix), aes(x=BlockAve_TF, fit, color = stim)) + 
#   facet_grid(~context,labeller=as_labeller(contextID)) + 
#   geom_line(aes(BlockAve_TF,fit)) +
#   geom_pointrange(aes(ymin=lower, ymax=upper,colour = stim), alpha= 0.5, show.legend=TRUE) +
#   #scale_color_manual(labels = c("Control","Gain","No Loss"), values = c("Black","Gray60","darkgreen")) +
#   labs(color = "Cue") +
#   ylab("Reaction Time (s)") +
#   xlab("Mean Gamma Power at Fixation") +
#   theme_bw() +
#   theme(strip.text = element_text(size = 12,family="Arial"), 
#         axis.text.x = element_text(size = 8,family="Arial"), 
#         axis.text.y = element_text(size = 14,family="Arial"), 
#         legend.text = element_text(size = 14,family="Arial"),
#         axis.title = element_text(size = 18,family="Arial"),
#         legend.title = element_text(size = 18,family="Arial"),
#         legend.position = "bottom") 
