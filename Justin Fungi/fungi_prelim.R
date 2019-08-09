library(tidyverse) 
library(ggplot2) #for plots
library(nlme)#for mixed effect models to test effects of treatments

# Import csv file, call it data
setwd("C:/Users/Owner/Desktop")
data<-read.csv("compost.fungi.csv",header=T)
str(data)
levels(data$ppt_trt)#check levels of precipitation treatment factor
levels(data$nut_trt)#check levels of nutrient treatment factor
levels(data$fungi)#check levels of fungi
data$block <- as.factor(data$block)
data$root <- as.factor(data$root)
data$rep <- as.factor(data$rep)

#colonization of amf by ppt, nut,root, and block
colonization <- data %>% group_by(block, ppt_trt, nut_trt, root, fungi) %>% filter(count != "NA") %>%
  summarize(percent=sum(count)/length(count))

#mean and standard deviation
col.plot.1 <- colonization %>% group_by(ppt_trt, nut_trt, fungi) %>%
  summarize(mean=mean(percent), stdev=sd(percent))

#make a graph!
ggplot(col.plot.1, aes(x=ppt_trt, y=mean, fill=fungi)) +
  facet_wrap(~nut_trt) +
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev), position=position_dodge(0.9))
