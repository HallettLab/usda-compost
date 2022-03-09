# clean and compile compost phenology dataset
# author(s): as (e.ashley.shaw@gmail.com), ctw (caitlin.t.white@colorado.edu)
# created: may 2019 (will be modified over project lifetime)

# script purpose:
# 1) read in cleaned main phenology dataset, look at trends over time within season and across years (as gather more data)
# 2) read in cleaned winter phenology dataset, look at trends in veg and litter height by all treatments and blocks
# > composition plots not raked in fall 2019 so standing and accumulated litter still present jan 2020
# > ctw also curious how grass vs. forb gs 2019 dominance influences winter 2020 phenology and litter



# -- SETUP ----
rm(list = ls()) # clean environment
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())

# specify dropbox pathway (varies by user)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## LMH and AS
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# read in main phenology datatset
pheno<-read.csv(paste0(datpath, "Phenology/Phenology_CleanedData/Compost_Phenology_Clean.csv"))
# read in winter 2020 phenology dataset
pheno_jan20 <- read.csv(paste0(datpath,"Phenology/Phenology_CleanedData/Compost_PhenologyJan2020_Clean.csv"))

# specify treatment colors for plotting
nut_cols <- c(C = "indianred4", `F` = "dodgerblue1", N = "darkgoldenrod")

# read in comp dat to ID dominant species and fxnl grp in each plot
compdat <- read.csv(paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_LongClean.csv"))



# -- PREP DATA -----
glimpse(pheno) # convert date to date
glimpse(pheno_jan20) # convert date to date
glimpse(compdat) # convert date todate

pheno$date <- as.Date(pheno$date, format = "%Y-%m-%d")
pheno_jan20$date <- as.Date(pheno_jan20$date)
compdat$date <- as.Date(compdat$date)

# summarize composition for plotting with pheno
# describe dom spp, dom fxnl grp, and keep pct gram and pct forb
domfxnl <- group_by(compdat, plot, block, ppt_trt, nut_trt) %>%
  summarise(mean_gram = mean(pct_grass),
         mean_forb = mean(pct_forb)) %>%
  ungroup() %>%
  mutate(dom = ifelse(mean_gram > mean_forb, "grass", "forb"))

domspp <- group_by(compdat, plot) %>%
  filter(pct_cover == max(pct_cover)) %>%
  ungroup() %>%
  dplyr::select(plot, code4:ncol(.)) %>%
  distinct()

# summarize winter pheno (avg litter and veg heights)
jan20_means <- dplyr::select(pheno_jan20,plot:notes) %>%
  gather(met, val, vht_1:lht_4) %>%
  mutate(rep = parse_number(met),
         met = ifelse(grepl("^v", met), "veg", "litter")) %>%
  rename_all(casefold) %>%
  grouped_df(names(.)[grep("pl|bl|nu|ppt|pct|not|met", names(.))]) %>%
  summarise(mean = mean(val),
            se = sd(val)/length(sqrt(val))) %>%
  ungroup()



# -- VISUALIZE MAIN (GROWING SEASON) PHENOLOGY DATASET -----
#show phenology overtime by nut_trt, ppt_trt and elevation
pheno_gf<-pheno %>% filter(plot!="NA")%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site),
         date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) %>%
  group_by(date, site, nut_trt, ppt_trt, yr) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

ggplot(data=pheno_gf, aes(x=date, y=meanPG, color=nut_trt))+
  facet_wrap(~site*ppt_trt*yr)+
  geom_point(cex=1.5)+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))

#show phenology for 2019 by nut_trt, ppt_trt and elevation
pheno_19<-pheno %>% filter(plot!="NA"&yr==2019)%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site),
         date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) %>%
  group_by(date, site, nut_trt, ppt_trt) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

ggplot(data=pheno_19, aes(x=date, y=meanPG, color=nut_trt))+
  facet_wrap(~site*ppt_trt)+
  geom_point(cex=1.5)+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))

#show phenology for 2019 by nut_trt, ppt_trt and elevation
pheno_20<-pheno %>% filter(plot!="NA"&yr==2020)%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site),
         date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) %>%
  group_by(date, site, nut_trt, ppt_trt) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

ggplot(data=pheno_20, aes(x=date, y=meanPG, color=nut_trt))+
  facet_wrap(~site*ppt_trt)+
  geom_point(cex=1.5)+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))

#show phenology for 2019 by nut_trt, ppt_trt and elevation
pheno_21<-pheno %>% filter(plot!="NA"&yr==2021)%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site),
         date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) %>%
  group_by(date, site, nut_trt, ppt_trt) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

ggplot(data=pheno_21, aes(x=date, y=meanPG, color=nut_trt))+
  facet_wrap(~site*ppt_trt)+
  geom_point(cex=1.5)+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))


#pheno plot for low dry plots only
pheno_gf2<-pheno %>% filter(plot!="NA")%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site),
         date = as.Date(date, format = "%Y-%m-%d")) %>%
  filter(site=="low")%>%
  filter(ppt_trt=="D")%>%
  group_by(date, nut_trt) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

ggplot(data=pheno_gf2, aes(x=date, y=meanPG, color=nut_trt))+
  geom_point(cex=1.5)+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=20))+
  theme_bw()+
  scale_color_manual(name = "Treatment", values = nut_cols,labels=c("Compost", "Fertilizer", "None"))



# -- WINTER 2020 PHENOLOGY ----
ggplot(jan20_means, aes(ppt_trt, mean, col = met)) +
  geom_errorbar(aes(ymax = mean + se, min = mean - se), position = position_dodge(width = 0.3), width = 0.2) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(y = "Mean height (cm)") +
  scale_color_manual(values = c(veg = "green", litter = "brown")) +
  theme_bw() +
  facet_grid(block ~ nut_trt)

# plot veg and litter heights by dominant fxnl grp instead of block
left_join(jan20_means, domfxnl) %>%
  ggplot(aes(ppt_trt, mean, col = met)) +
  geom_errorbar(aes(ymax = mean + se, min = mean - se), position = position_dodge(width = 0.3), width = 0.2) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(y = "Mean height (cm)") +
  scale_color_manual(values = c(veg = "chartreuse4", litter = "burlywood3")) +
  theme_bw() +
  facet_grid(dom ~ nut_trt)

# plot by dominant 2019 spp 
left_join(jan20_means, domspp) %>%
  ggplot(aes(ppt_trt, mean, col = met, group = block)) +
  geom_errorbar(aes(ymax = mean + se, min = mean - se), position = position_dodge(width = 0.5), width = 0) +
  geom_point(position = position_dodge(width = 0.5)) +
  labs(y = "Mean height (cm)") +
  scale_color_manual(values = c(veg = "chartreuse4", litter = "burlywood3")) +
  theme_bw() +
  facet_grid(fxnl_grp ~ nut_trt) #code4  .. not as interesting. Tall veg in FW TRSU/N-fix-dominant is Stipha pulchra

# plot pct green, brown, etc (usual pheno vars)
distinct(dplyr::select(pheno_jan20, plot:pct_bare, block:ppt_trt)) %>%
  gather(met, val, pct_litter:pct_bare) %>%
  ggplot( aes(ppt_trt, val, col = met)) +
  #geom_errorbar(aes(ymax = mean + se, min = mean - se), position = position_dodge(width = 0.3), width = 0.2) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(y = "Percent") +
  scale_color_manual(values = c(pct_litter = "burlywood3", pct_green = "chartreuse4", pct_bare = "brown", pct_brown = "pink")) +
  theme_bw() +
  facet_grid(block ~ nut_trt)

# plot by dominant species fxnl grp
distinct(dplyr::select(pheno_jan20, plot:pct_bare, block:ppt_trt)) %>%
  gather(met, val, pct_litter:pct_bare) %>%
  left_join(domspp) %>%
  ggplot( aes(ppt_trt, val, col = met)) +
  #geom_errorbar(aes(ymax = mean + se, min = mean - se), position = position_dodge(width = 0.3), width = 0.2) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(y = "Percent") +
  scale_color_manual(values = c(pct_litter = "burlywood3", pct_green = "chartreuse4", pct_bare = "brown", pct_brown = "pink")) +
  theme_bw() +
  facet_grid(fxnl_grp ~ nut_trt)


# plot means by dom fxnl grp
distinct(dplyr::select(pheno_jan20, plot:pct_bare, block:ppt_trt)) %>%
  gather(met, val, pct_litter:pct_bare) %>%
  left_join(domspp) %>%
  group_by(fxnl_grp, nut_trt, ppt_trt, met) %>%
  summarise(meancov = mean(val),
          secov = sd(val)/sqrt(length(val)),
          nobs = length(val)) %>%
  ungroup() %>%
  ggplot(aes(ppt_trt, meancov, col = met)) +
  geom_errorbar(aes(ymax = meancov+secov, ymin = meancov - secov), position = position_dodge(width = 0.2), width = 0.1) +
  geom_point( position = position_dodge(width = 0.2)) +
  labs(y = "Percent cover") +
  scale_color_manual(name = NULL, values = c(pct_litter = "burlywood3", pct_green = "chartreuse4", pct_bare = "brown", pct_brown = "pink")) +
  theme_bw() +
  facet_grid(fxnl_grp ~ nut_trt)

