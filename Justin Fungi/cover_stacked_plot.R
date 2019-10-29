##run AMF reorg.R file first
library(OneR)
plant4$bin<-bin(plant4$mean, nbins = 2, labels = c("low", "high"), na.omit = FALSE)
plant4$richbin<-bin(plant4$richness, nbins = 3, labels = c("low","med", "high"), na.omit = FALSE)

bin<-plant4 %>% select(ppt_trt, nut_trt, block, bin, richbin)

cover_long<-right_join(bin,plant.data, by=c("ppt_trt","nut_trt","block"))

sp_graph <- cover_long %>%
  group_by (ppt_trt, richbin, genus, fxnl_grp) %>%
  summarise(cover=mean(pct_cover), secover=sd(pct_cover)/sqrt(length(pct_cover)))

sp_graph <- sp_graph %>% arrange(fxnl_grp) %>% filter(cover>3)%>%filter(richbin!="NA")%>%filter(genus!="NA")%>% ungroup()%>%
  mutate(ppt_trt2=ordered(ppt_trt, levels = c( d="d",xc="xc", w="w"))) 
levels(sp_graph$ppt_trt2) <- c( "d"="Drought","xc"="Ambient", "w"="Wet")
sp_graph$genus<- factor(sp_graph$genus, c("Aphanes", "Cerastium", "Chamomilla", "Erodium","Navarretia","Triphysaria", "Triteleia","Aira", "Avena", "Briza","Bromus","Hordeum", "Juncus","Lolium", "Stipa", "Taeniatherum", "Vulpia","Trifolium", "Vicia"))
ggplot(sp_graph, aes(fill=genus, colour=fxnl_grp,  y=cover, x=richbin)) +
  facet_wrap(~ppt_trt2)+
  theme_bw()+
  scale_fill_manual(values = c("yellow", "orange", "orangered", "firebrick","indianred4", "indianred3", "indianred1", "darkolivegreen1","seagreen1","palegreen3",  "springgreen4","yellowgreen", "green2", "green4","darkolivegreen","darkgreen", "darkseagreen4", "skyblue2","royalblue3" ))+
  geom_bar( stat="identity", position='stack')+
  xlab("Richness Bin")+ 
  ylab("Mean cover (%)")+ 
  ggtitle("Abundant plant species by richness")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(plant4, aes(x=Grass, y=mean, color=ppt_trt, shape=nut_trt))+
  geom_point()

plant3[is.na(plant3)] <- 0
plant4<-merge(plant3,plant4)

ggplot(plant4, aes(x=`Avena barbata` , y=mean, color=ppt_trt, shape=ppt_trt))+
  geom_point()+
  geom_smooth(method="lm", se=F)+ 
  facet_wrap(~ppt_trt) 

ggplot(plant4, aes(x=Grass , y=mean, color=ppt_trt, shape=ppt_trt))+
  geom_point()+
  geom_smooth(method="lm", se=F)+ 
  facet_wrap(~ppt_trt)

ggplot(plant4, aes(x=(`Trifolium hirtum`+ `Trifolium subterraneum`), y=mean, color=ppt_trt, shape=nut_trt))+
  geom_point()+
  geom_smooth(method="lm", se=F)
  

ggplot(plant4, aes(x=(Forb + Grass + nfixer) , y=mean, color=ppt_trt))+
  geom_point()+
  geom_smooth(method="lm", se=F)+ 
  facet_wrap(~ppt_trt)

ggplot(plant4, aes(x=`Hordeum murinum` , y=mean, color=ppt_trt))+
  geom_point()+
  geom_smooth(method="lm", se=F)+ 
  facet_wrap(~ppt_trt)

ggplot(plant4, aes(x=`Triphysaria eriantha` , y=mean, color=ppt_trt))+
  geom_point()+
  geom_smooth(method="lm", se=F)+ 
  facet_wrap(~ppt_trt*nut_trt)


