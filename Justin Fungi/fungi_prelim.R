library(tidyverse) 
library(ggplot2) #for plots
library(nlme)#for mixed effect models to test effects of treatments
library(lsmeans)#post hoc test for significance between treatments

# Import csv file, call it data. Import soil moisture data, call it moisture.data
setwd("C:/Users/Owner/Desktop")
data<-read.csv("compost.fungi.csv",header=T) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(d="d", xc="xc", w="w"))) %>% #orders factors
  mutate(nut_trt=ordered(nut_trt, levels = c(c="c", f="f", n="n"))) #orders factors
  
str(data)
levels(data$ppt_trt)#check levels of precipitation treatment factor
levels(data$nut_trt)#check levels of nutrient treatment factor
levels(data$fungi)#check levels of fungi
data$block <- as.factor(data$block)
data$root <- as.factor(data$root)
data$rep <- as.factor(data$rep)

#import soil moisture data
moisture.data <- read.csv("moisture.csv", header=T) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(d="d", xc="xc", w="w"))) %>% #orders factors
  mutate(nut_trt=ordered(nut_trt, levels = c(c="c", f="f", n="n"))) %>%
  
str(moisture.data)
moisture.data$block <- as.factor(moisture.data$block)
levels(moisture.data$block)
levels(moisture.data$ppt_trt)
levels(moisture.data$nut_trt)

#import root biomass data (belowground net primary productivity, BNPP)
BNPP <- read.csv("BNPP.csv", header=T) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(d="d", xc="xc", w="w"))) %>% #orders factors
  mutate(nut_trt=ordered(nut_trt, levels = c(c="c", f="f", n="n"))) 
  
str(BNPP)
BNPP$block <- as.factor(BNPP$block)
levels(BNPP$block)
levels(BNPP$ppt_trt)
levels(BNPP$nut_trt)

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

#JD: description - bar plot of only AMF colonization by nut and ppt reatments.
ggplot(subset(col.plot.1, fungi=="amf"), aes(x=nut_trt, y=mean, fill=ppt_trt)) +
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9)) +
  ylab("AMF % root colonization")+
  xlab("soil amendment")+
  ggtitle("Percent AMF colonization Across Soil and Precipitation Treatments")+
  scale_x_discrete(labels=c("Compost", "Fertilizer","No Amendment")) +
  scale_fill_manual(values = c("indianred1",  "lightgoldenrod2", "skyblue2"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                    guide = guide_legend(title = "Amendment"), #change legend title
                    labels=c("Drought", "Ambient", "High")) #change labels in the legend

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
#AS: nice job joining these!! but I think I would use the average values (all 5 roots averaged per block)
#AS: I added col.plot.2 to average colonization, leaving block in 
#AS: Then I joined moisture to the averaged colonization data in col.moist.plot2
col.moist.plot <- full_join(colonization, moisture.data)
col.plot.2 <- colonization %>% group_by(block, ppt_trt, nut_trt, fungi) %>%
  summarize(mean=mean(percent), stdev= sd(percent), se=sd(percent)/sqrt(length(percent)))
col.moist.plot2 <- full_join(col.plot.2, moisture.data)

#add BNPP data to colonization and moisture data
col.moist.plot2<-merge(col.moist.plot2, BNPP)

#plotting moisture and AMF colonization. whoops! didn't work, or don't really know how to code?
#AS: I changed this to a scatterplot (geom_point) and added a linear regression by nutrient treatment (geom_smooth)
#AS: I also added some code to make a more customized, prettier plot 
#AS: feel free to use this as a template for adjusting other plots/making pretty plots for the poster
ggplot(subset(col.moist.plot2,fungi=="amf"), aes(y=mean,x=percent_moisture, color=nut_trt))+
  geom_point()+ #plots points for scatterplot
  geom_smooth(method="lm", se=F)+ #adds linear regression to the plot
  facet_wrap(~nut_trt)+ #creates 3 panels
  ylab("AMF (% colonization)")+ #change y-axis label
  xlab("Soil Moisture (% g/g)")+ #change x-axis label
  theme_classic() + #a nicer theme without gray background
  scale_color_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                     guide = guide_legend(title = "Amendment"), #change legend title
                     labels=c("Compost", "Fertilizer", "None")) #change labels in the legend
#colors() will give you a list of colors, or google "r colors"

ggplot(moisture.stat,aes(x=nut_trt, y=mean, fill=ppt_trt))+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))+
  ylab("Soil Moisture (%)")+ #change y-axis label
  xlab("Amendment Treatment") #change x-axis label

ggplot(moisture.stat,aes(x=ppt_trt, y=mean, fill=nut_trt))+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))+
  ylab("Soil Moisture (%)")+ #change y-axis label
  xlab("Amendment Treatment") #change x-axis label

##AS: I THINK THIS ONE FOR POSTER!!
##Plot AMF colonization vs. root biomass
ggplot(subset(col.moist.plot2,fungi=="amf"), aes(y=mean,x=BNPP, color=nut_trt))+
  geom_point()+ #plots points for scatterplot
  facet_wrap(~nut_trt)+
  geom_smooth(method="lm", se=F)+ #adds linear regression to the plot
  ylab("AMF (% colonization)")+ #change y-axis label
  xlab("Root Biomass (g))")+ #change x-axis label
  theme_classic() + #a nicer theme without gray background
  scale_color_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                     guide = guide_legend(title = "Amendment"), #change legend title
                     labels=c("Compost", "Fertilizer", "None")) #change labels in the legend
#colors() will give you a list of colors, or google "r colors"

##Plot root biomass vs. soil moisture
ggplot(subset(col.moist.plot2,fungi=="amf"), aes(y=BNPP,x=percent_moisture, color=nut_trt))+
  geom_point()+ #plots points for scatterplot
  geom_smooth(method="lm", se=F)+ #adds linear regression to the plot
  facet_wrap(~nut_trt)+
  ylab("Root Biomass (g)")+ #change y-axis label
  xlab("Soil Moisture(%))")+ #change x-axis label
  theme_classic() + #a nicer theme without gray background
  scale_color_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                     guide = guide_legend(title = "Amendment"), #change legend title
                     labels=c("Compost", "Fertilizer", "None")) #change labels in the legend

#boxplot to show root biomass by treatments (root biomass is highest in compost/wet)
ggplot(col.moist.plot2,aes(x=ppt_trt, y=BNPP, fill=nut_trt))+
  geom_boxplot() +
  ylab("Root biomass (g)")+ #change y-axis label
  xlab("Precipitation Treatment") #change x-axis label

#JD - Boxplot to show root biomass by treatment (precipitation fill and soil amendment x)
ggplot(col.moist.plot2,aes(x=nut_trt, y=BNPP, fill=ppt_trt))+
  geom_boxplot() +
  ylab("Root biomass (g)")+ #change y-axis label
  xlab("Soil Amendment")+ #change x-axis label
  ggtitle("Root Biomass (g) and Precipitation Treatments")+
  scale_x_discrete(labels=c("Compost", "Fertilizer","No Amendment")) +
  scale_fill_manual(values = c("indianred1",  "lightgoldenrod2", "skyblue2"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                    guide = guide_legend(title = "Amendment"), #change legend title
                    labels=c("Drought", "Ambient", "High")) #change labels in the legend


#ANOVA for nut_trt*percent_moisture on percent colonization
#I'm not entirely sure that I did this analysis correctly
#AS: This is correct for a linear model, no significant effects though :(
amf.moist <- col.moist.plot2 %>% filter(fungi=="amf")
options(contrasts = c("contr.treatment", "contr.poly"))
m2 = lm ( mean ~ nut_trt + percent_moisture + nut_trt:percent_moisture,
              data = amf.moist)
summary(m2)
anova(m2)

##AS: test if nut*root biomass affect AMF colonization in the COMPOST treatment
m3 <- lm(mean ~  BNPP, 
         data = subset(amf.moist, nut_trt=="c"))
summary(m3)
anova(m3) #yes, this is significant!
#we can conclude that AMF significantly declines with increasing root biomass
#under the compost treatment only (see below for FERT and NO amendments, which are
#not significant. 
#linear equation is AMF(%) = 0.96 - 0.75(BNPP)
#R2=0.33 (translates as 33% of variation in AMF colonization is explained by BNPP under compost)

##AS: test if nut*root biomass affect AMF colonization in the FERT treatment
m4 <- lm(mean ~  BNPP, 
         data = subset(amf.moist, nut_trt=="f"))
summary(m4)
anova(m4) #not significant!

##AS: test if nut*root biomass affect AMF colonization in the NO AMEND treatment
m5 <- lm(mean ~  BNPP, 
         data = subset(amf.moist, nut_trt=="n"))
summary(m5)
anova(m5) #not significant!
