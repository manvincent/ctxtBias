setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex2/Results/SubmitAnalyses/Demographics')

library("vcd")
library("plyr")
library("lsr")
library("effects")
library("ggplot2")

demographics_data <- read.csv("Intro_Form_Study2_FinalSample_Scored.csv", header = T)
# Cross tab the nominal data for Chi-squared tests
# Compute counts
nominal_counts <- count(demographics_data, "Sex"); nominal_counts

attach(demographics_data)
# get AGE descriptives
mean(demographics_data$Age)
sd(demographics_data$Age)
min(demographics_data$Age)
max(demographics_data$Age)

# get years edu descriptives 
mean(demographics_data$Years.of.education..beginning.from.Grade.1, na.rm=T)
sd(demographics_data$Years.of.education..beginning.from.Grade.1, na.rm=T)
min(demographics_data$Years.of.education..beginning.from.Grade.1, na.rm=T)
max(demographics_data$Years.of.education..beginning.from.Grade.1, na.rm=T)
