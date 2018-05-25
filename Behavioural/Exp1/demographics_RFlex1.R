setwd('/Users/vincentman/Dropbox/Graduate_School/current_projects/RFlex/Studies/RFlex1/Results')
demographics_data <- read.csv('Rflex1_Final_Demographics.csv',header=T)

library("vcd")
library("plyr")
library("lsr")
library("effects")
library("ggplot2")
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
mean(demographics_data$Education, na.rm=T)
sd(demographics_data$Education, na.rm=T)
min(demographics_data$Education, na.rm=T)
max(demographics_data$Education, na.rm=T)

