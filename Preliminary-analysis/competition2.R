# read in data
dat <- read_csv("~/Desktop/competition_seeds_summer2021_phyto.csv")
tri <- read_csv("~/Desktop/trifolium_seeds_competition_summer2021.csv")
dat <- dat %>% mutate(ps=as.numeric(seeds)/as.numeric(stems))
tri <- tri %>% group()
ggplot(dat, aes(x=phyto, y = ps, fill=nut_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  theme_bw()+
  #facet_wrap(~background)+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
ggplot(dat, aes(x=phyto, y = ps, fill=nut_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  theme_bw()+
  facet_wrap(~ppt_trt)+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
ggplot(dat, aes(x=background, y = ps, fill=nut_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  theme_bw()+
  facet_wrap(~phyto)+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
ggplot(dat, aes(x=background, y = ps, fill=phyto)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  theme_bw()+
  facet_wrap(~ppt_trt)#+
ggplot(dat, aes(x=background, y = ps, fill=nut_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  theme_bw()+
  facet_wrap(~phyto)+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
