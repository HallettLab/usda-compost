library(lubridate)
library(tidyverse)

setwd("~/Dropbox/USDA-compost/Data")
dat<-read.csv("SoilMoisture/SoilMoisture_CleanedData/SoilMosture_all_clean.csv")

#show vwc by amendment treatment across precip treatment for low blocks
smdat2 <- dat %>%
  #filter(comp_trt=="control")%>%
  filter(vwc>0)%>%
  filter(block!=3)%>%
  group_by(nut_trt, ppt_trt, date, block) %>%
  summarize(sm = mean(vwc))

ggplot(smdat2, aes(x=date, y=sm, 
                 group = interaction(nut_trt, ppt_trt),
                 color = nut_trt)) + facet_wrap(~ppt_trt) + geom_line()

#show low drought only
smdat3 <- smdat2 %>%
  tbl_df() %>%
  filter(ppt_trt=="D")%>%
  group_by(nut_trt, block, date) %>%
  summarize(sm = mean(sm, na.rm=T)) 

ggplot(smdat3, aes(x=date, y=sm, color = nut_trt))  + geom_line(lwd = 1.5) +
  #facet_wrap(~ppt_trt)+
  theme_bw() + 
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  labs(y = expression(paste("Soil moisture (m"^3, " m"^-3,")")), x = "Day of 2018-2019 growing season") +
  scale_x_discrete(breaks = levels(smdat3$date)[c(T, rep(F, 30))])+
  scale_color_manual(name="Amendment",
                     breaks=c("C", "F", "N"),
                     labels=c("Compost", "Fertilizer", "None"),
                     values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"))

smdat4 <- dat %>%
  filter(comp_trt=="control")%>%
  filter(vwc>0)%>%
  filter(nut_trt!="C")%>%
  tbl_df() %>%
  group_by(ppt_trt, block) %>%
  summarize(sm = mean(vwc, na.rm=T)) 
  
