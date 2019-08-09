library(tidyverse) 
library(ggplot2) #for plots
library(nlme)#for mixed effect models to test effects of treatments
library(lsmeans)#post hoc test for significance between treatments

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
  summarize(mean=mean(percent), stdev= sd(percent), se=sd(percent)/length(percent))

#make a graph!
ggplot(col.plot.1, aes(x=ppt_trt, y=mean, fill=fungi)) +
  facet_wrap(~nut_trt) +
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))

#Another graph!
ggplot(col.plot.1, aes(x=ppt_trt, y=mean, fill=nut_trt)) +
  facet_wrap(~fungi) +
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))

#statistical significance test
colonization <- as.data.frame(colonization)
m1 <- lme(subset(colonization, fungi=="amf"), percent ~ nut_trt*ppt_trt, random=~1|root/block, na.action=na.exclude)

col.amf <- colonization %>% filter(fungi=="amf")
col.amf <- as.data.frame(col.amf)
m1 <- lme(percent ~ nut_trt*ppt_trt, random=~1|block,col.amf, na.action=na.exclude)
summary(m1)
anova(m1)

m1.lsm <- lsmeans(m1, ~nut_trt*ppt_trt)
contrast(m1.lsm, "pairwise")

m1.ppt <- lsmeans(m1, ~ppt_trt)
contrast(m1.ppt, "pairwise")

#new graph!
ggplot(subset(colonization, fungi=="amf"), aes(x=nut_trt, y=percent, fill=ppt_trt))+
  geom_boxplot()
