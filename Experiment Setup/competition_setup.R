library(tidyverse)
#Shelter layout
shelterlayout<-read_csv("~/Dropbox/USDA-compost/site-map/treatment_plots.csv")
#Identify and randomly fill remaining subplots
mytrick<-as.data.frame(cbind(plot=rep(seq(1,36,1),each=6), subplot=rep(seq(1,6,1),6)))
tog<-merge(shelterlayout, mytrick, all=T)
spptreatment<-as.vector(replicate(36,
                                  sample(c("ERBO", "AVBA", "HOMU", "LOMU",  "TRHI", "TACA", "Control"))))
tog$background<-spptreatment
View(tog)
#Identify and randomly fill remaining subplots
mytrick<-as.data.frame(cbind(plot=rep(seq(1,36,1),each=6), subplot=rep(seq(1,7,1),6)))
#Identify and randomly fill remaining subplots
mytrick<-as.data.frame(cbind(plot=rep(seq(1,36,1),each=6), subplot=rep(seq(1,7,1),7)))
#Identify and randomly fill remaining subplots
mytrick<-as.data.frame(cbind(plot=rep(seq(1,36,1),each=7), subplot=rep(seq(1,7,1),7)))
#Identify and randomly fill remaining subplots
mytrick<-as.data.frame(cbind(plot=rep(seq(1,36,1),each=7), subplot=rep(seq(1,7,1),7)))
#Identify and randomly fill remaining subplots
mytrick<-as.data.frame(cbind(plot=rep(seq(1,36,1),each=7), subplot=rep(seq(1,7,1),7)))
#Identify and randomly fill remaining subplots
mytrick<-as.data.frame(cbind(plot=rep(seq(1,36,1),each=7), subplot=rep(seq(1,6,1),7)))
View(mytrick)
#Identify and randomly fill remaining subplots
mytrick<-as.data.frame(cbind(plot=rep(seq(1,36,1),each=7), subplot=rep(seq(1,7,1))))
View(mytrick)
tog<-merge(shelterlayout, mytrick, all=T)
spptreatment<-as.vector(replicate(36,
                                  sample(c("ERBO", "AVBA", "HOMU", "LOMU",  "TRHI", "TACA", "Control"))))
tog$background<-spptreatment
View(tog)
mytrick2<- as.data.frame(cbind(plot=rep(seq(1,36,1),each=36), subplot=rep(seq(1,7,1),each=6),phytonum=rep(seq(1,6,1), 6)))
tog2<-merge(tog, mytrick2, all=T)
View(tog2)
mytrick2<- as.data.frame(cbind(plot=rep(seq(1,36,1),each=36), subplot=rep(seq(1,7,1),each=6),phytonum=rep(seq(1,6,1), 7)))
View(mytrick2)
mytrick2<- as.data.frame(cbind(plot=rep(seq(1,36,1),each=36), subplot=rep(seq(1,7,1),each=6),phytonum=rep(seq(1,6,1))))
mytrick2<- as.data.frame(cbind(plot=rep(seq(1,36,1),each=42), subplot=rep(seq(1,7,1),each=6),phytonum=rep(seq(1,6,1))))
View(mytrick2)
tog2<-merge(tog, mytrick2, all=T)
View(tog2)
spptreatment2<-as.vector(replicate(252,
                                   sample(c("ERBO", "AVBA", "HOMU", "LOMU",  "TRHI", "TACA"))))
tog2$phyto<-spptreatment2
View(tog2)
#reorder dataframe by column name
layout <- tog2[c("block", "nut_trt", "ppt_trt", "plot", "subplot","background","phytonum","phyto")]
View(layout)
layout <-layout %>% arrange(plot, subplot, phytonum)
View(layout)
#write.csv(layout, "~/Dropbox/USDA-compost/site-map/plotlayout.csv")
