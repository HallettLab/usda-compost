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

#colonization of amf by ppt, nut,root, and block
colonization <- data %>% group_by(block, ppt_trt, nut_trt, root, fungi) %>% filter(count != "NA") %>%
  summarize(percent=sum(count)/length(count))

#mean and standard deviation
col.plot.1 <- colonization %>% group_by(ppt_trt, nut_trt, fungi) %>%
  summarize(mean=mean(percent), stdev= sd(percent), se=sd(percent)/sqrt(length(percent)))


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

#ANOVA for nut_trt*percent_moisture on percent colonization
#I'm not entirely sure that I did this analysis correctly
#AS: This is correct for a linear model, no significant effects though :(
amf.moist <- col.moist.plot2 %>% filter(fungi=="amf")
options(contrasts = c("contr.treatment", "contr.poly"))
m2 = lm ( mean ~ nut_trt + percent_moisture + nut_trt:percent_moisture,
          data = amf.moist)
summary(m2)
anova(m2)

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

plant2 <- full_join(amf.moist, plant1)

#species data/ diversity
plant3 <- plant.data%>%
  dplyr::select(block, ppt_trt, nut_trt, species, pct_cover, date)%>%
  filter(date!="2019-04-19", date!="2019-04-20")%>%
  spread(species, pct_cover)

cover <- plant3%>%
  dplyr::select(5:56)

cover[is.na(cover)] <- 0

plant2$diversity <- diversity(cover)

#richness
plant2$richness <- specnumber(cover)

#Evenness diversity
#
#Needs Debugging
plant2$evenness <- plant2$diversity/log(specnumber(cover))

#functional group
plant4 <- plant.data%>%
  dplyr::select(block, ppt_trt, nut_trt, fxnl_grp, pct_cover, date)%>%
  filter(date!="2019-04-19", date!="2019-04-20")%>%
  mutate(ppt_trt=ordered(ppt_trt, levels=c("d","xc","w")))%>%
  mutate(nut_trt=ordered(nut_trt, levels=c("n","f","c")))
  
levels(plant4$ppt_trt)
levels(plant4$nut_trt)

plant4 <- plant4%>%  
  group_by(block, nut_trt, ppt_trt, fxnl_grp)%>%
  summarise(percent=sum(pct_cover))%>%
  spread(fxnl_grp, percent)

plant4 <- merge(plant4, plant2)

plant4 <- plant4%>%
  select(-pct.grass, -pct.forb)

str(plant4)
colnames(plant4)[colnames(plant4) == "N-fixer"] <- "nfixer"


#PLANT COMPOSITION STATS
#
#

#ANOVA  for AMF and diversity
#significant relationship between AMF colonization and diversity (AMF decline with increasing diversity)
p1 = lme ( mean ~ diversity, random=~1|block, plant4, na.action=na.exclude)
summary(p1)
anova(p1)

ggplot(plant4, aes(x=diversity, y=mean))+
  geom_point()+
  geom_smooth(method="lm")

#ANOVA  for AMF and richness
#significant relationship between AMF colonization and diversity (AMF decline with increasing diversity)
p1a = lme ( mean ~ richness, random=~1|block, plant4, na.action=na.exclude)
summary(p1a)
anova(p1a)

ggplot(plant4, aes(x=richness, y=mean))+
  geom_point()+
  geom_smooth(method="lm")

#ANOVA for AMF and Forb
#no significance
p2 = lme ( mean ~ Forb, random=~1|block, plant4, na.action=na.exclude)
summary(p2)
anova(p2)

#ANOVA for AMF and Grass
#no significance
p3 = lme ( mean ~ Grass, random=~1|block, plant4, na.action=na.exclude)
summary(p3)
anova(p3)

#ANOVA for AMF and N-fixer
#significant effects of nfixers on AMF, where AMF declines with increasing Nfixer cover
p4 = lme ( mean ~ nfixer, random=~1|block, plant4, na.action=na.exclude)
summary(p4)
anova(p4)

ggplot(plant4, aes(x=nfixer, y=mean))+
  geom_point()+
  geom_smooth(method="lm")

#ANOVA for AMF and evenness
#ANOVA for AMF richness

#ANOVA for forb and treatment
#significance for diversity X nut_trt, but not combined treatments
m1 = lm (diversity ~ ppt_trt + nut_trt + ppt_trt:nut_trt,
         data = plant4)
summary(m1)
anova(m1)

#ANOVA for nfixer and treatment
#no significance
m2 = lm (nfixer ~ ppt_trt + nut_trt + ppt_trt:nut_trt,
         data = plant4)
summary(m2)
anova(m2)

#ANOVA for forb and treatment
#no significance
m3 = lm (Forb ~ ppt_trt + nut_trt + ppt_trt:nut_trt,
         data = plant4)
summary(m3)
anova(m3)

#ANOVA for grass and treatment
#no significance
m4 = lm (Grass ~ ppt_trt + nut_trt + ppt_trt:nut_trt,
         data = plant4)
summary(m4)
anova(m4)

#across treatments
q1 = lme ( mean ~ diversity*nut_trt*ppt_trt, random=~1|block, plant4, na.action=na.exclude)
summary(q1)
anova(q1)

#Richness and AMF
#Significant intercept with ppt_trt
q2 = lme ( richness ~ mean*nut_trt*ppt_trt, random=~1|block, plant4, na.action=na.exclude)
summary(q2)
anova(q2)

#difference in richness with colonization
lq2 <- lsmeans(q2, ~mean*ppt_trt)
contrast(lq2, "pairwise")

#evenness
q3 = lme ( evenness ~ mean*nut_trt*ppt_trt, random=~1|block, plant4, na.action=na.exclude)
summary(q3)
anova(q3)


#FiGURES
#
#

#new data set for plots specifically
plot_data <- plant4

plot_data<- plot_data %>% mutate(nut_trt=ifelse(nut_trt=="c", "Compost", 
                                                           ifelse(nut_trt=="f", "Fertilizer", 
                                                                  ifelse(nut_trt=="n", "No Amendment", nut_trt))))
plot_data<- plot_data %>% mutate(ppt_trt=ifelse(ppt_trt=="d", "Drought", 
                                                ifelse(ppt_trt=="xc", "Ambient", 
                                                       ifelse(ppt_trt=="w", "Wet", ppt_trt))))
plot_data <- plot_data%>%
  mutate(ppt_trt=ordered(ppt_trt, levels=c("Drought","Ambient","Wet")))%>%
  mutate(nut_trt=ordered(nut_trt, levels=c("No Amendment","Fertilizer","Compost")))
  
str(plot_data)
levels(plot_data$ppt_trt)
levels(plot_data$nut_trt)
#diversity*amf
ggplot(subset(plot_data,fungi=="amf"), aes(y=diversity,x=mean))+
  geom_point()+ 
  geom_smooth(method="lm", se=F)+ 
  facet_wrap(~nut_trt)+
  ylab("diversity")+ 
  xlab("AMF colonization (% root)")+ 
  ggtitle("AMF vs. diversity")+
  theme_classic() + 
  theme(legend.position="none", axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"), strip.text.x = element_text(size = 16))

#nfixer*amf
ggplot(subset(plot_data,fungi=="amf"), aes(y=mean,x=nfixer))+
  geom_point()+ 
  geom_smooth(method="lm", se=F)+ 
  ylab("AMF colonization")+ 
  xlab("nitrogen fixers")+ 
  ggtitle("AMF vs. nfixer")+
  theme_classic() + 
  theme(legend.position="none", axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"), strip.text.x = element_text(size = 16))

#richness*amf
ggplot(subset(plot_data,fungi=="amf"), aes(y=richness,x=mean, color=ppt_trt))+
  geom_point()+ 
  geom_smooth(method="lm", se=F)+ 
  facet_wrap(~ppt_trt)+
  xlab("AMF colonization")+ 
  ylab("richness")+ 
  ggtitle("Regression of Plot Richness with AMF Colonization")+
  theme_classic()+
  theme(legend.position="none", axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"), strip.text.x = element_text(size = 16))+
  scale_color_manual(values = c( "indianred1","lightgoldenrod2","skyblue2" ), 
                     guide = guide_legend(title = "Precipitation Treatment"),
                     labels=c("Drought", "Ambient", "High"))


#diversity*nutrients
ggplot(plot_data,aes(x=nut_trt, y=diversity))+
  geom_bar(stat="identity", position="dodge") +
  ylab("diversity")+
  xlab("")+
  ggtitle("")+
  scale_x_discrete(labels=c("Compost", "Fertilizer","No Amendment")) +
  theme(legend.position=c(0.8,0.8), legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"))+


#evenness*amf
ggplot(subset(plot_data,fungi=="amf"), aes(y=evenness,x=mean))+
  geom_point()+ 
  geom_smooth(method="lm", se=F)+ 
  ylab("evenness")+ 
  xlab("AMF colonization (% root)")+ 
  ggtitle("AMF vs. evenness")+
  theme_classic() + 
  theme(legend.position="none", axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(size = 18, face = "bold"), strip.text.x = element_text(size = 16))

#Boxplot of amf colonization across nut and ppt treatments.
ggplot(plot_data, aes(x=nut_trt, y=mean, fill=ppt_trt))+
  geom_boxplot()+
  scale_fill_manual(values = c( "indianred1","lightgoldenrod2","skyblue2" ), 
                     guide = guide_legend(title = "Precipitation Treatment"),
                     labels=c("Drought","Ambient", "Wet"))

#Boxplot of BNPP across nut and ppt treatments.
ggplot(plot_data, aes(x=nut_trt, y=BNPP, fill=ppt_trt))+
  geom_boxplot()+
  scale_fill_manual(values = c( "indianred1","lightgoldenrod2","skyblue2" ), 
                    guide = guide_legend(title = "Precipitation Treatment"),
                    labels=c( "Drought", "Ambient","Wet"))
