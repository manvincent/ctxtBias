setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/EEG/FFT')
library( nlme )
library(RColorBrewer)
display.brewer.all()
brewer.pal(9,"Greens")
#colorRampPalette(brewer.pal(n,palette))(n of spectrum)
colorRampPalette(brewer.pal(5,"BuPu"))(100)
library("lme4")
library("lmerTest")


freq_list = c('delta','theta','lowGamma')

for (currBand in freq_list) {  
  dataFFT <- read.csv(paste('ContextMeans_ContrastLossVal_f',currBand,'.csv',sep=""),header=T)
  print(currBand)
  attach(dataFFT)
  # Average across channels
  dataFFT$Channels <- rowMeans( dataFFT[c(-1:-3,-(length(dataFFT)-2):-length(dataFFT))])
 
  # Simple effects test! 
  # Recoding for simple effects (Aiken and West 1991)
  dataFFT$GainContrast <- ifelse(dataFFT$context == "OO" |  dataFFT$context == "OF", 1, ifelse(dataFFT$context == "FO" | dataFFT$context == "FF", -1, 0))
  dataFFT$LossContrast <- ifelse(dataFFT$context == "OO" |  dataFFT$context == "FO", 1, ifelse(dataFFT$context == "OF" | dataFFT$context == "FF", -1, 0))
  dataFFT$Informativeness <- ifelse(dataFFT$context == "OF" |  dataFFT$context == "FO", 1, ifelse(dataFFT$context == "OO" | dataFFT$context == "FF", -1, 0))
  
  
  mod_simple1 <- lmer(Channels ~ GainContrast + LossContrast + Informativeness + (1 | subj_idx), data=dataFFT)
  print(anova(mod_simple1))
  detach(dataFFT)
  }  

  