library(tidyverse) 
library(ggplot) #for plots
library(nlme)#for mixed effect models to test effects of treatments

# Import csv file, call it data
setwd("~/Desktop/PhD Research/Manuscripts/In prep/LTM")
data<-read.csv("fungi.csv",header=T)
str(data)
levels(data$ppt_trt)#check levels of precipitation treatment factor
levels(data$nut_trt)#check levels of nutrient treatment factor
levels(data$fungi)#check levels of fungi


