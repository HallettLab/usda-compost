library(tidyverse) 
library(ggplot2) #for plots
library(nlme)#for mixed effect models to test effects of treatments
library(lsmeans)#post hoc test for significance between treatments

# Import csv file, call it data. Import soil moisture data, call it moisture.data
setwd("C:/Users/Owner/Desktop")
data<-read.csv("compost.fungi.csv",header=T)
str(data)
levels(data$ppt_trt)#check levels of precipitation treatment factor
levels(data$nut_trt)#check levels of nutrient treatment factor
levels(data$fungi)#check levels of fungi
data$block <- as.factor(data$block)
data$root <- as.factor(data$root)
data$rep <- as.factor(data$rep)

moisture.data <- read.csv("moisture.csv", header=T)
str(moisture.data)
moisture.data$block <- as.factor(moisture.data$block)
levels(moisture.data$block)
levels(moisture.data$ppt_trt)
levels(moisture.data$nut_trt)

#colonization of amf by ppt, nut,root, and block
colonization <- data %>% group_by(block, ppt_trt, nut_trt, root, fungi) %>% filter(count != "NA") %>%
  summarize(percent=sum(count)/length(count))

#mean and standard deviation
col.plot.1 <- colonization %>% group_by(ppt_trt, nut_trt, fungi) %>%
  summarize(mean=mean(percent), stdev= sd(percent), se=sd(percent)/sqrt(length(percent)))

#make a graph! 
#AS: description - bar plot with error bars of mean colonization by nut and ppt treatment, colored by fungi
ggplot(col.plot.1, aes(x=ppt_trt, y=mean, fill=fungi)) +
  facet_wrap(~nut_trt) +
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))

#Another graph!
#AS: description - bar plot with error bars of mean colonization by ppt and fungi, colored by nut trt
ggplot(col.plot.1, aes(x=ppt_trt, y=mean, fill=nut_trt)) +
  facet_wrap(~fungi) +
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))

#AS: we use a mixed model to test if nutrient and ppt treatments affect AMF colonization
#AS: in this model, nutrient and ppt treatment are fixed effects, block is a random effect
col.amf <- colonization %>% filter(fungi=="amf")
col.amf <- as.data.frame(col.amf)
m1 <- lme(percent ~ nut_trt*ppt_trt, random=~1|block,col.amf, na.action=na.exclude)
summary(m1)
anova(m1)

#AS: this is a post-hoc test to show differences between treatments (interaction term)
m1.lsm <- lsmeans(m1, ~nut_trt*ppt_trt)
contrast(m1.lsm, "pairwise")

#AS: this is a post-hoc test to show differences in ppt treatment only
m1.ppt <- lsmeans(m1, ~ppt_trt)
contrast(m1.ppt, "pairwise")

#new graph! 
#AS: description - boxplot of only AMF colonization by nut and ppt treatments
ggplot(subset(colonization, fungi=="amf"), aes(x=nut_trt, y=percent, fill=ppt_trt))+
  geom_boxplot()

#formating moisture.data. Calculating soil moisture
#AS: I changed the formula to calculate % water out of DRY soil 
moisture.data$dry_wt <- moisture.data$dry_soil_tin - moisture.data$tin_wt
moisture.data$water_wt <- moisture.data$wet_soil - moisture.data$dry_wt
moisture.data$percent_moisture <- (moisture.data$water_wt / moisture.data$dry_soil) * 100 #changed to dry soil 

#mean, sd, and se of soil moisture data
#AS: fixed error in se calculation (needed square root of n, my mistake on thursday)
moisture.stat <- moisture.data %>% group_by(ppt_trt, nut_trt) %>%
  summarize(mean=mean(percent_moisture), se=sd(percent_moisture)/sqrt(length(percent_moisture)))

#add soil moisture to colonization data
#AS: nice job!
col.moist.plot <- full_join(colonization, moisture.data)

#plotting moisture and AMF colonization. whoops! didn't work, or don't really know how to code?
#AS: I changed this to a scatterplot (geom_point) and added a linear regression by nutrient treatment (geom_smooth)
#AS: I also added some code to make a more customized, prettier plot - use this as a template for making pretty plots for the poster
ggplot(subset(col.moist.plot,fungi=="amf"), aes(y=percent,x=percent_moisture, color=nut_trt))+
  geom_point()+ #plots points for scatterplot
  geom_smooth(method="lm", se=F)+ #adds linear regression to the plot
  ylab("AMF (% colonization)")+ #change y-axis label
  xlab("Soil Moisture (% g/g)")+ #change x-axis label
  theme_classic() + #a nicer theme without gray background
  scale_color_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                     guide = guide_legend(title = "Amendment"), #change legend title
                     labels=c("Compost", "Fertilizer", "None")) #change labels in the legend

ggplot(moisture.stat,aes(x=nut_trt, y=mean, fill=ppt_trt))+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))

ggplot(moisture.stat,aes(x=ppt_trt, y=mean, fill=nut_trt))+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))
