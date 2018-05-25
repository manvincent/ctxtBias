setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/EEG/FFT')
data <- read.csv('RFS2_MLM.csv',header=T)

library( nlme )
library(RColorBrewer)
brewer.pal(9,"Greens")
colorRampPalette(brewer.pal(5,"BuPu"))(100)
library("lme4")
library("lmerTest")
library("effects")
library("ggplot2")
source('sortData.R')


freq_list = c('delta','theta','lowGamma')

for (currBand in freq_list) {  
  dataFFT <- read.csv(paste('ContextMeans_ContrastLossVal_f',currBand,'.csv',sep=""),header=T)
  print(currBand)

  
  # Create dataframe of only channel values 
  dataChan = dataFFT[c(-1:-3,-(length(dataFFT)-2):-length(dataFFT))]
  dataContrast = dataFFT[c(length(dataFFT)-2):length(dataFFT)]
  dataSub = dataFFT[1]
  
  for (ctrstNo in c(2)) {
    currCtrst = dataContrast[ctrstNo]
    currDataFrame = cbind(dataSub, dataChan, currCtrst)
    # Pull out contrast-based sum of contexts for each channel (for each subject) 
    cntrstSum <- aggregate(dataChan, by=list(dataSub[,], currCtrst[,]), FUN=sum)
    
    # Pull out contrast-based differences between context codes for each channel (for each subjet)             
    cntrstDiff <- aggregate(cntrstSum, by=list(cntrstSum[,1]), FUN=diff)
    
    # Average across channels
    cntrstAve =  as.data.frame(matrix(NA, ncol = 1, nrow = dim(cntrstDiff)[1]))
    cntrstAve$V1 = cntrstDiff$Group.1
    colnames(cntrstAve) <- "subj_idx"
    
    cntrstAve$aveAllChan <- rowMeans( cntrstDiff[c(-1:-3)])
  }
  
  # Using all subjects
  # Remove outliers +/- 3 sd in mean channel power differences
  cntrstAve$aveAllChan[cntrstAve$aveAllChan > 3*sd(cntrstAve$aveAllChan,na.rm=TRUE) | cntrstAve$aveAllChan < -3*sd(cntrstAve$aveAllChan,na.rm=TRUE)] <- NA
  
  # Mean-centre the data
  cntrstAve$aveAllChan <- scale(cntrstAve$aveAllChan,scale=F)
  dataRaw <-merge(data,cntrstAve,by.x = "subj_idx")
  
  # Clean and sort out data
  # Data sorting function
  dataClean <- sortData(dataRaw,removeNA=TRUE) 
  
  ########## Look at if RT is moderated by EEG signal ################ 
  # Only accurate trials: rt ~ context*stim
  # Parse data
  dataAcc <- dataClean[ which(dataClean$Accuracy_inTime==1), ]
  
  
  ## Testing as continuous regressor - does mean gamma signal predict mean RT across contexts
  # Get average rts per context condition
  dataMeanRT_context <- aggregate(dataAcc$rt, by = list(dataAcc$subj_idx,dataAcc$context,dataAcc$stim), FUN=mean )
  names(dataMeanRT_context) <- c("subj_idx","context","stim","meanRT")
  dataContTest <- merge(dataMeanRT_context, dataAcc)
  
  # Run lme model
  mod_accRT_EEG_cont <- lmer(meanRT ~ context*aveAllChan + (1 | subj_idx/contextNo), data=dataContTest)
  print(anova(mod_accRT_EEG_cont))
  
  # Compute high and low s.d. splits
  lowMeanChans <- mean(dataContTest$aveAllChan,na.rm=TRUE)-sd(dataContTest$aveAllChan,na.rm=TRUE)
  highMeanChans <- mean(dataContTest$aveAllChan,na.rm=TRUE)+sd(dataContTest$aveAllChan,na.rm=TRUE)
  
  # Plot effect
  ef_accRT_EEG_cont <- effect('context:aveAllChan', xlevels=2, mod_accRT_EEG_cont, confidence.level = 0.95)
  summary(ef_accRT_EEG_cont)
  
  ggplot(as.data.frame(ef_accRT_EEG_cont), aes(x=aveAllChan, fit)) + 
    geom_line(aes(color=context),size=2,show.legend=FALSE) +
    geom_linerange(aes(ymin=lower, ymax=upper,colour = context), alpha= 0.8,size=1, show.legend=TRUE) +
    geom_point(data=dataContTest,aes(x=aveAllChan,y=meanRT,color=context)) + 
    scale_color_manual(labels = c("LM","Prev","Prom","HM"), values = c("Black","darkgreen","gray75","gray50")) +
    scale_linetype_manual(breaks = c("LM","Prev","Prom","HM"), values = c(2, 1, 4, 5)) + 
    labs(color = "Context") +
    ylab("Mean Reaction Time (s)") +
    xlab(expression('Individual Mean'~Delta~'Difference')) + 
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    theme(strip.text = element_text(size = 12,family="Arial"), 
          axis.text.x = element_text(size = 12,family="Arial"), 
          axis.text.y = element_text(size = 14,family="Arial"), 
          legend.text = element_text(size = 14,family="Arial"),
          axis.title = element_text(size = 18,family="Arial"),
          legend.title = element_text(size = 18,family="Arial"),
          legend.position = "bottom") 
  
  
  
  
  # Simple effects 1: Main effect of stim (gain v loss) in the high loss (prev, high mag) contexts but not low loss (prom, low mag)
  dataContTest$prevCtxt <- ifelse(dataContTest$context=="FO",0,1)
  dataContTest$highCtxt <- ifelse(dataContTest$context=="OO",0,1)
  dataContTest$promCtxt <- ifelse(dataContTest$context=="OF",0,1)
  dataContTest$lowCtxt <- ifelse(dataContTest$context=="FF",0,1)
  
  dataContTest$highVlowFFT <- ifelse(dataContTest$aveAllChan > highMeanChans, 1, ifelse(dataContTest$aveAllChan < lowMeanChans, -1, 0))
  # Run simple effect tests
  modAcc_simple1a <- lmer(meanRT ~ highVlowFFT*prevCtxt + (1 | subj_idx/contextNo), data=dataContTest); summary(modAcc_simple1a)
  modAcc_simple1b <- lmer(meanRT ~ highVlowFFT*highCtxt + (1 | subj_idx/contextNo), data=dataContTest); summary(modAcc_simple1b)
  modAcc_simple1c <- lmer(meanRT ~ highVlowFFT*promCtxt + (1 | subj_idx/contextNo), data=dataContTest); summary(modAcc_simple1c)
  modAcc_simple1d <- lmer(meanRT ~ highVlowFFT*lowCtxt + (1 | subj_idx/contextNo), data=dataContTest); summary(modAcc_simple1d)
  
  # Template for computing mixed-effects model effect sizes
  tVal <- # model t avlue
  tdf <- # model t DoF
  d <- (2*tVal/sqrt(tdf)); d
  
  FVal <- tVal^2     
  dfN <- # model F DoF numerator
  dfD <- tdf
  eta <- ((dfN/dfD) * FVal)/(1+((dfN/dfD) * FVal)); eta
  
  

  
  
  

###### Testing as categorical (split 2 groups)  #######
dataContTest$powDiffCode = ifelse(dataContTest$aveAllChan > 0,1,-1)

mod_accRT_EEG <- lmer(meanRT ~ context*powDiffCode + (1 | subj_idx/contextNo), data=dataContTest)
anova(mod_accRT_EEG)


# Plot effects
ef_accRT_EEG <- effect('context:powDiffCode', mod_accRT_EEG, xlevels=2, confidence.level = 0.95)
summary(ef_accRT_EEG)
# Facet labels: 
fftID <- c(
  `-1` = "Gamma Desynchrony to Contextual Loss",
  `1` = "Gamma Synchrony to Contextual Loss"
)
ggplot(as.data.frame(ef_accRT_EEG), aes(context, fit, width = 0.75)) +
  facet_grid( ~ powDiffCode, labeller = as_labeller(fftID)) +
  geom_col(aes(fill=context),color="black",size=1,width=0.4,show.legend = FALSE, position="dodge") +
  geom_linerange(aes(ymin=lower, ymax=upper),size=0.8,position=position_dodge(width=0.4)) + 
  scale_fill_manual(values = c( "gray20","darkgreen","gray75","gray50")) +
  scale_x_discrete(labels = c("LM", "Prev", "Prom", "HM"))+
  ylab("Mean Reaction Time (s)") +
  xlab("Context") +
  theme_bw() +
  theme(strip.text = element_text(size = 12, family = "Arial"),
        axis.text.x = element_text(size = 14, family = "Arial",angle=90),
        axis.text.y = element_text(size = 14, family = "Arial"),
        legend.text = element_text(size = 14, family = "Arial"),
        axis.title = element_text(size = 18, family = "Arial"),
        legend.title = element_text(size = 18, family = "Arial"),
        legend.position = "none" ) + 
  coord_cartesian(ylim=c(0.305, 0.328))

# Simple effects: Main effect of stim (gain v loss) in the high loss (prev, high mag) contexts but not low loss (prom, low mag)
dataContTest$prevVprom <- ifelse(dataContTest$context == "FO", 1, ifelse(dataContTest$context == "OF", -1, 0))


dataContTest$gammaDesync <- ifelse(dataContTest$powDiffCode==-1,0,1)
dataContTest$gammaSync <- ifelse(dataContTest$powDiffCode==1,0,1)


modAcc_simple2a <- lmer(meanRT ~ gammaDesync*prevVprom + (1 | subj_idx/contextNo),data=dataContTest); summary(modAcc_simple2a)
modAcc_simple2b <- lmer(meanRT ~ gammaSync*prevVprom + (1 | subj_idx/contextNo), data=dataContTest); summary(modAcc_simple2b)

}
