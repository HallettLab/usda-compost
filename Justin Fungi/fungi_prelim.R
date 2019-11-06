library(tidyverse) 
library(ggplot2) #for plots
library(nlme)#for mixed effect models to test effects of treatments
library(lsmeans)#post hoc test for significance between treatments
library(vegan)

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
  mutate(nut_trt=ordered(nut_trt, levels = c(c="c", f="f", n="n"))) 
  
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

#import plant composition data
plant.data <- read.csv("Compost_Cover_LongClean.csv", header=T)
levels(plant.data$ppt_trt) <- c("D"="d","W"="w","XC"="xc")#Change factors to lower case
levels(plant.data$nut_trt) <- c("C"="c", "F"="f", "N"="n")

str(plant.data)
plant.data$block <- as.factor(plant.data$block)
levels(plant.data$block)
levels(plant.data$ppt_trt)
levels(plant.data$nut_trt)
levels(plant.data$fxnl_grp)
levels(plant.data$Duration)
levels(plant.data$nativity)
levels(plant.data$date)

#percent grass/forb
plant1 <- plant.data%>%
  dplyr::select(block, nut_trt, ppt_trt, pct_grass, pct_forb, pct_bare, pct_litter, litter_depth_cm)%>%
  group_by(block, ppt_trt, nut_trt)%>%
  summarise(pct.grass=max(pct_grass), pct.forb = max(pct_forb), pct.bare = max(pct_bare), pct.litter=max(pct_litter),litter.depth.cm=max(litter_depth_cm))

plant2 <- merge(amf.moist, plant1)

#species data/ diversity
plant3 <- plant.data%>%
  dplyr::select(block, ppt_trt, nut_trt, species, pct_cover, date)%>%
  filter(date!="2019-04-19", date!="2019-04-20")%>%
  spread(species, pct_cover)

cover <- plant3%>%
  dplyr::select(5:56)

cover[is.na(cover)] <- 0

plant2$diversity <- diversity(cover)

#functional group
plant4 <- plant.data%>%
  dplyr::select(block, ppt_trt, nut_trt, fxnl_grp, pct_cover, date)%>%
  filter(date!="2019-04-19", date!="2019-04-20")

plant4 <- plant4%>%  
  group_by(block, nut_trt, ppt_trt, fxnl_grp)%>%
  summarise(percent=sum(pct_cover))%>%
  spread(fxnl_grp, percent)

plant4 <- merge(plant4, plant2)

plant4 <- plant4%>%
  select(-pct.grass, -pct.forb)

#Add plant composition data to colonization by plot
all.data <- merge(plant.data, amf.moist)

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
col.amf$z<-scale(col.amf$percent) #try converting amf percent to z scores
m1 <- lme(percent ~ nut_trt*ppt_trt, random=~1|block,col.amf, na.action=na.exclude)
summary(m1)
anova(m1)
shapiro.test(residuals(m1)) #shapiro test for normality, normal if p>0.05, these data are not normally distributed

#histogram of nut_trt 
ggplot(data=col.amf, aes(x=percent, group=nut_trt, color=nut_trt))+
  geom_density()
  #geom_vline(aes(xintercept = grp.mean, group=nut_trt), linetype = "dashed", size = 0.6)

#histogram of ppt_trt 
ggplot(data=col.amf, aes(x=percent, group=ppt_trt, color=ppt_trt))+
  geom_density()

#histogram of all data, play around with transformations
#note closest to normal was asin transformation (also tried log, log(x+1), sqrt, asin(x^0.5))
ggplot(data=col.amf, aes(x=percent))+
  geom_density()

#AS: this is a post-hoc test to show differences between treatments (interaction term)
m1.lsm <- lsmeans(m1, ~nut_trt*ppt_trt)
contrast(m1.lsm, "pairwise")

#AS: this is a post-hoc test to show differences in ppt treatment only
m1.ppt <- lsmeans(m1, ~ppt_trt)
contrast(m1.ppt, "pairwise")

#new graph! 
#AS: description - boxplot of only AMF colonization by nut and ppt treatments
ggplot(subset(colonization, fungi=="amf"), aes(x=ppt_trt, y=percent, fill=nut_trt))+
  geom_boxplot()

#scatterplot of moisture vs AMF
ggplot(plant4, aes(x=percent_moisture, y=mean, color=ppt_trt, shape=nut_trt))+
  geom_point()+
  geom_smooth(method="lm", aes(group=1))

#scatterplot of moisture vs root biomass
ggplot(plant4, aes(x=percent_moisture, y=BNPP, color=ppt_trt, shape=nut_trt))+
  geom_point()+
  geom_smooth(method="lm", aes(group=1))

#JD: description - bar plot of only AMF colonization by nut and ppt reatments.
ggplot(subset(col.plot.1, fungi=="amf"), aes(x=nut_trt, y=mean, fill=ppt_trt)) +
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9)) +
  ylab("AMF (% root colonization)")+
  xlab("")+
  ylim(0, 1.1)+
  theme(legend.position=c(0.85,0.9),legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"))+
  ggtitle("AMF colonization across nutrient and precipitation treatments")+
  scale_x_discrete(labels=c("Compost", "Fertilizer","No Amendment")) +
  scale_fill_manual(values = c("indianred1",  "lightgoldenrod2", "skyblue2"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                    guide = guide_legend(title = "Precipitation"), #change legend title
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

#JD BNPP mean, sd, se
BNPP.stat <- BNPP %>% group_by(nut_trt, ppt_trt)%>%
  summarize(mean=mean(BNPP), stdev= sd(BNPP), se=sd(BNPP)/sqrt(length(BNPP)))

#add BNPP data to colonization and moisture data
col.moist.plot2<-merge(col.moist.plot2, BNPP)

#plotting moisture and AMF colonization. whoops! didn't work, or don't really know how to code?
#AS: I changed this to a scatterplot (geom_point) and added a linear regression by nutrient treatment (geom_smooth)
#AS: I also added some code to make a more customized, prettier plot 
#AS: feel free to use this as a template for adjusting other plots/making pretty plots for the poster
col.moist.plot2<-col.moist.plot2 %>% mutate(nut_trt=ifelse(nut_trt=="c", "Compost", 
                                                            ifelse(nut_trt=="f", "Fertilizer", 
                                                                   ifelse(nut_trt=="n", "No Amendment", nut_trt))))
ggplot(subset(col.moist.plot2,fungi=="amf"), aes(y=mean,x=percent_moisture, color=nut_trt))+
  geom_point()+ #plots points for scatterplot
  geom_smooth(method="lm", se=F)+ #adds linear regression to the plot
  facet_wrap(~nut_trt)+ #creates 3 panels
  ylab("AMF (% root colonization)")+ #change y-axis label
  xlab("Soil Moisture (% g/g)")+ #change x-axis label
  theme_classic() +#a nicer theme without gray background
  theme(legend.position="none", axis.text=element_text(size=12))+
  ggtitle("AMF colonization by soil moisture")+
  scale_color_manual(values = c("green3", "orange", "skyblue3"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                     guide = guide_legend(title = "Amendment"), #change legend title
                     labels=c("Compost", "Fertilizer", "None")) #change labels in the legend
#colors() will give you a list of colors, or google "r colors"
moisture.stat<-moisture.stat%>%mutate(nut_trt=ifelse(nut_trt=="c", "Compost", 
                                                     ifelse(nut_trt=="f", "Fertilizer", 
                                                            ifelse(nut_trt=="n", "No Amendment", nut_trt))))
ggplot(moisture.stat,aes(x=nut_trt, y=mean, fill=ppt_trt))+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))+
  ylab("Soil Moisture (%)")+ #change y-axis label
  xlab("")+ #change x-axis label
  theme(legend.position=c(0.1,0.9),legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"))+
  ggtitle("Soil moisture across nutrient and precipitation treatments")+
  scale_x_discrete(labels=c("Compost", "Fertilizer","No Amendment")) +
  scale_fill_manual(values = c("indianred1",  "lightgoldenrod2", "skyblue2"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                    guide = guide_legend(title = "Precipitation"), #change legend title
                    labels=c("Drought", "Ambient", "High")) #change labels in the legend


##JD - barplot showing soil moisture against all treatments
ggplot(moisture.stat,aes(x=ppt_trt, y=mean, fill=nut_trt))+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))+
  ylab("Soil Moisture (%)")+ #change y-axis label
  xlab("Precipitation Treatment")+#change x-axis label
  scale_x_discrete(labels=c("Drought", "Ambient","High")) +
  theme(legend.position=c(0.2,0.8), legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"))+
  ggtitle ("Soil moisture across nutrient and precipitation treatments")+
  scale_fill_manual(values = c("green3", "orange", "skyblue3"),
                     guide = guide_legend(title = "Nutrient Treatment"), #change legend title
                     labels=c("Compost", "Fertilizer", "None")) 

##AS: I THINK THIS ONE FOR POSTER!!
##Plot AMF colonization vs. root biomass
ggplot(subset(col.moist.plot2,fungi=="amf"), aes(y=mean,x=BNPP, color=nut_trt))+
  geom_point()+ #plots points for scatterplot
  facet_wrap(~nut_trt)+
  geom_smooth(method="lm", se=F)+ #adds linear regression to the plot
  ylab("AMF (% colonization)")+ #change y-axis label
  xlab("Root Biomass (g)")+ #change x-axis label
  ggtitle("Regression of AMF colonization on root biomass")+
  theme_classic() + #a nicer theme without gray background
  theme(legend.position="none", axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"), strip.text.x = element_text(size = 16))+
  scale_color_manual(values = c("green3", "orange", "skyblue3"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
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

#JD - Barplot to show root biomass by treatment (precipitation fill and soil amendment x)
ggplot(BNPP.stat,aes(x=nut_trt, y=mean, fill=ppt_trt))+
  geom_bar(stat="identity", position="dodge") +
  ylab("Root biomass (g)")+ #change y-axis label
  xlab("")+ #change x-axis label
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9))+
  ggtitle("Root biomass across nutrient and precipitation treatments")+
  scale_x_discrete(labels=c("Compost", "Fertilizer","No Amendment")) +
  theme(legend.position=c(0.8,0.8), legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"))+
  scale_fill_manual(values = c("indianred1",  "lightgoldenrod2", "skyblue2"), #changes colors of points (use scale_fill_manual for boxplots and bar plots)
                    guide = guide_legend(title = "Precipitation Treatment"), #change legend title
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

#AMF colonization against:
#functional group (fxnl_grp)
#percent cover of forbs/ grass (pct_forb, pct grass)
#percent cover of each species
#plot composition?
#nativity