soil <- read_csv("~/Downloads/Soil CN Data USDA (2).csv")
soil <- soil %>%
  mutate(site = "Moderate",
         site = ifelse(block == "B1" | block == "B2", "Heavy", site))

soil$ppt_trt_f = factor(soil$subtrtm, levels=c('D','CN','W'))
soil$nut_trt = factor(soil$treatment, levels=c('C','F','CN'))

ggplot(subset(soil, year!="2018"&depth_category!="30-40cm"&depth_category!="40-50cm"), aes(x=ppt_trt_f, y = N_percent, fill=treatment)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~depth_category)+
  labs(x="", y="Nitrogen (%)") +
  scale_fill_manual(values = c("indianred4",   "darkgoldenrod", "dodgerblue1"), guide = guide_legend(title = "Amendment"))#, labels=c("Compost", "Fertilizer", "None"))

m1<-lme(sqrt(N_percent) ~treatment, random=~1|depth_category/year, subset(soil, year!="2018"&depth_category!="30-40cm"&depth_category!="40-50cm"), na.action=na.exclude)
summary(m1)
anova(m1)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))

LS1<-lsmeans(m1, ~treatment)
contrast(LS1, "pairwise")
#significant effect of nut_trt


soil2<- subset(soil, depth_category!="30-40cm"&depth_category!="40-50cm") %>% group_by(treatment, subtrtm)%>%
  summarize(N=mean(N_percent), SE=sd(N_percent)/sqrt(length(N_percent)))
soil2$ppt_trt_f = factor(soil2$subtrtm, levels=c('W','CN','D'))

ggplot(subset(soil2), aes(x=treatment, y = N, color=treatment)) + 
  geom_point(size=4)+ # facet_wrap(~ppt_trt)+
  geom_errorbar(aes(ymax = N+SE, ymin = N-SE), width=.25)+
  #theme_bw()+
  facet_wrap(~ppt_trt_f)+
  scale_x_discrete(labels=c("C" = "Compost", "CN" = "None", "F" = "Fertilizer"))+
  labs(x="", y="Nitrogen (%)") +
  scale_color_manual(values = c("indianred4",   "darkgoldenrod", "dodgerblue1"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

soil3<- subset(soil, depth_category!="30-40cm"&depth_category!="40-50cm") %>% group_by(treatment, subtrtm,depth_category)%>%
  summarize(N=mean(N_percent), SE=sd(N_percent)/sqrt(length(N_percent)))
soil3$ppt_trt_f = factor(soil3$subtrtm, levels=c('D','CN','W'))

ggplot(subset(soil3), aes(x=depth_category, y = N, color=treatment)) + 
  geom_point()+ 
  geom_errorbar(aes(ymax = N+SE, ymin = N-SE), width=.25)+
  #scale_x_discrete(labels=c("C" = "Compost", "CN" = "None", "F" = "Fertilizer"))+
  #theme_bw()+
  facet_wrap(~ppt_trt_f)+
  labs(x="", y="Nitrogen (%)") +
  scale_color_manual(values = c("indianred4",   "darkgoldenrod", "dodgerblue1"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

soil4<- subset(soil, depth_category!="30-40cm"&depth_category!="40-50cm") %>% group_by(treatment)%>%
  summarize(N=mean(N_percent), SE=sd(N_percent)/sqrt(length(N_percent)))

ggplot(subset(soil4), aes(x=treatment, y = N, color=treatment)) + 
  geom_point(size=4)+ # facet_wrap(~ppt_trt)+
  geom_errorbar(aes(ymax = N+SE, ymin = N-SE), width=.25)+
  #theme_bw()+
  #facet_wrap(~ppt_trt_f)+
  labs(x="", y="Nitrogen (%)") +
  scale_x_discrete(labels=c("C" = "Compost", "CN" = "None", "F" = "Fertilizer"))+
  scale_color_manual(values = c("indianred4",   "darkgoldenrod", "dodgerblue1"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(subset(soil, year!="2018"&depth_category!="30-40cm"&depth_category!="40-50cm"), aes(x=treatment, y = C_percent, fill=treatment)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~depth_category)+
  labs(x="", y="Carbon (%)") +
  scale_fill_manual(values = c("indianred4",   "darkgoldenrod", "dodgerblue1"), guide = guide_legend(title = "Amendment"))#, labels=c("Compost", "Fertilizer", "None"))

m2<-lme(sqrt(C_percent) ~treatment*subtrtm, random=~1|depth_category/year, subset(soil, year!="2018"&depth_category!="30-40cm"&depth_category!="40-50cm"), na.action=na.exclude)
summary(m2)
anova(m2)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))

LS1<-lsmeans(m2, ~treatment)
contrast(LS1, "pairwise")
#significant effect of nut_trt


soilC<- subset(soil, depth_category!="30-40cm"&depth_category!="40-50cm") %>% group_by(treatment, subtrtm)%>%
  summarize(C=mean(C_percent), SE=sd(C_percent)/sqrt(length(C_percent)))
soilC$ppt_trt_f = factor(soilC$subtrtm, levels=c('D','CN','W'))
soilC$nut_trt_f = factor(soilC$treatment, levels=c('C','F','CN'))

f2f<-ggplot(subset(soilC), aes(x=ppt_trt_f, y = C, color=treatment)) + 
  geom_point(size=4, position=position_dodge(0.9))+  #facet_wrap(~ppt_trt)+
  geom_errorbar(aes(ymax = C+SE, ymin = C-SE), width=.25, position=position_dodge(0.9))+
  #theme_bw()+
  #facet_wrap(~ppt_trt_f)+
  #scale_x_discrete(labels=c("C" = "Compost", "CN" = "None", "F" = "Fertilizer"))+
  labs(x="", y="") +
  scale_x_discrete(labels=c("D" = "Dry", "CN" = "Ambient", "W" = "Wet"))+
  scale_color_manual(values = c("indianred4", "darkgoldenrod","dodgerblue1"), guide = guide_legend(title = "Amendment"), labels=c("Compost",  "None", "Fertilizer"))
f2f

soilC3<- subset(soil, depth_category!="30-40cm"&depth_category!="40-50cm") %>% group_by(treatment,depth_category)%>%
  summarize(C=mean(C_percent), SE=sd(C_percent)/sqrt(length(C_percent)))
soilC3$ppt_trt_f = factor(soilC3$subtrtm, levels=c('D','CN','W'))

f2e<-ggplot(subset(soilC3), aes(x=depth_category, y = C, color=treatment)) + 
  geom_point(size=4,position=position_dodge(0.9))+ 
  geom_errorbar(aes(ymax = C+SE, ymin = C-SE), width=.25, position=position_dodge(0.9))+
  #scale_x_discrete(labels=c("C" = "Compost", "CN" = "None", "F" = "Fertilizer"))+
  #theme_bw()+
  #facet_wrap(~ppt_trt_f)+
  labs(x="", y="") +
  #scale_x_discrete(labels=c("D" = "Dry", "CN" = "Ambient", "W" = "Wet"))+
  scale_color_manual(values = c("indianred4",   "darkgoldenrod", "dodgerblue1"), guide = guide_legend(title = "Amendment"), labels=c("Compost",  "None", "Fertilizer"))
f2e

soilC4<- subset(soil, depth_category!="30-40cm"&depth_category!="40-50cm") %>% group_by(treatment)%>%
  summarize(C=mean(C_percent), SE=sd(C_percent)/sqrt(length(C_percent)))

f2d<-ggplot(subset(soilC4), aes(x=treatment, y = C, color=treatment)) + 
  geom_point(size=4)+ # facet_wrap(~ppt_trt)+
  geom_errorbar(aes(ymax = C+SE, ymin = C-SE), width=.25)+
  #theme_bw()+
  #facet_wrap(~ppt_trt_f)+
  labs(x="", y="Carbon (%)") +
  scale_x_discrete(labels=c("C" = "Compost", "CN" = "None", "F" = "Fertilizer"))+
  scale_color_manual(values = c("indianred4",   "darkgoldenrod", "dodgerblue1"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "None", "Fertilizer"))
f2d

p_col2<-plot_grid(
  f2e + theme(legend.position="none"), f2f + theme(legend.position="none"),
  labels = c("E","F"), ncol = 1
)


f2D<-plot_grid(
  f2d + theme(legend.position="none"), p_col2,
  labels = c("D", "", ""), ncol = 3)
f2D

f2_all<-plot_grid(
  f2, f2D, ncol = 1)
f2_all

f3b<-ggplot(subset(soil, depth_category=="0-10cm"), aes(x=treatment, y = C_percent, fill=treatment)) + geom_boxplot() + facet_grid(~site)+
  #theme_bw()+
  labs(x="",y="Carbon (%)")+
  ylim(1,4)+
  scale_x_discrete(labels=c("C" = "Compost", "CN" = "None", "F" = "Fertilizer"))+
  scale_fill_manual(values = c("indianred4",   "darkgoldenrod","dodgerblue1"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "None", "Fertilizer"))
f3b

p_site<-plot_grid(
  f3a + theme(legend.position="none"), f3b + theme(legend.position="none"),
  labels = c("A","B"), ncol = 1, align = "v"
)
p_site
