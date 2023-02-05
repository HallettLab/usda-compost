# native multivariate analyses
# author(s): ctw
# questions: caitlin.t.white@colorado.edu


# -- SETUP -----
# load needed libraries
library(tidyverse)
library(vegan)





# -- MULTIVARIATE TO ASSESS NEIGHBORHOOD RESPONSE TO SEEDING ----- 

# prep species abundance dataset with native seeded removed
# need companion environmental dataset with treatments per plot (ppt, nut, herbicided or not, seeded or not, hillslope position)

# assess tot cov before running NMDSes to see how removing natives seeded impacts totcov
# or can just run on relativized and unrelativized to compare if any major diffs
allbutnats <- subset(natlong, !code4 %in% nats, select = c(plot:mon, code4,pct_cover)) %>%
  # choose max cov within season
  group_by(plot, herbicide, seedtrt, code4) %>%
  summarise(pct_cover = max(pct_cover), .groups = "drop_last")

# prep cover -- choose yr 2021 only, and max cov by species (not sure how I can adjust size bc cov still relative)
sppcomp2021 <- subset(coverlong, yr ==2021 & plot != 33, select = c(plot:sample_event, code4, pct_cover)) %>%
  group_by(plot, code4) %>%
  summarise(pct_cover = max(pct_cover), .groups = "drop_last") %>%
  mutate(herbicide = "Non-herbicided", seedtrt = "Unseeded")

# join to all but nats
allbutnats <- rbind(allbutnats, sppcomp2021) %>%
  distinct() %>%
  # remove rare (present in fewer than 5% of plots)
  group_by(code4) %>%
  mutate(freq = length(pct_cover[pct_cover > 0])) %>%
  ungroup() %>%
  mutate(plots = length(unique(paste(plot, herbicide, seedtrt))),
         prop_pres = freq/plots) %>%
  subset(prop_pres >=0.05, select = -c(freq, plots, prop_pres)) %>%
  spread(code4, pct_cover, fill = 0) %>%
  #create unique ID code for each subplot
  mutate(natplotid = 1:nrow(.))

# make enviro dat
allbutnats_env <- dplyr::select(allbutnats, natplotid, plot:seedtrt) %>%
  # join main experiment trtment info
  left_join(trtkey) %>%
  # make hillpos col again
  mutate(hillpos = ifelse(block < 3, "low", "high"),
         hillpos = factor(hillpos, levels = c("low", "high")))

# prep dat for nmds and permanova
summary(rownames(allbutnats) == allbutnats$natplotid) # rowname is already what natplotid is, so don't need to keep col
allbutnats_cov <- dplyr::select(allbutnats, select = -c(plot:seedtrt, natplotid))
# check plotcov
allbutnats_env$totcov <- apply(allbutnats_cov,1, sum)
# add native seeded cover
allbutnats_env2 <- left_join(allbutnats_env, subset(neighborhood, select =c(plot, herbicide, seedtrt, BRCA:NEMA))) %>%
  # fill NAs for seeded spp with 0? or -1 for not there (no chance)
  replace_na(list(FEMI = -1, BRCA = -1, NEMA = -1, ESCA = -1))
# make PA dat too to try to draw hulls
allbutnats_envPA <- mutate_at(allbutnats_env2, .vars = c("ESCA", "FEMI", "BRCA", "NEMA"), function(x) ifelse(x>0, "P","A"))

allbutnats_bc <- vegdist(allbutnats_cov, method = "bray")
allbutnats_nmds <- metaMDS(allbutnats_cov, k = 3, trymax = 100)
# try relativized also
allbutnats_nmds_rel <- metaMDS(decostand(allbutnats_cov, method = "total"), trymax = 100)
plot(allbutnats_nmds)
plot(allbutnats_nmds_rel)
stressplot(allbutnats_nmds)
stressplot(allbutnats_nmds_rel)


ordiplot(allbutnats_nmds, type = "n")
orditorp(allbutnats_nmds, display="species",col="red", air=0.01)
orditorp(allbutnats_nmds, display="sites",col="grey50", air=0.01)
ordisurf(allbutnats_nmds,allbutnats_env2$FEMI,col="forestgreen")
ordisurf(allbutnats_nmds,allbutnats_env2$NEMA,col="chocolate3")
ordisurf(allbutnats_nmds,allbutnats_env2$ESCA,col="orchid")
ordisurf(allbutnats_nmds,allbutnats_env2$BRCA, display = "sites", col="blue") #brca is linear! on axes 1 and 2
ordihull(allbutnats_nmds, paste(allbutnats_env2$seedtrt, allbutnats_env2$herbicide), label = T, col = 1:4, lwd=3)

# axes 1 and 2
ordiplot(allbutnats_nmds, type = "n")
orditorp(allbutnats_nmds, display="species",col="red", air=0.01)
orditorp(allbutnats_nmds, display="sites",col="grey50", air=0.01)
ordiellipse(allbutnats_nmds, paste(allbutnats_env2$seedtrt, allbutnats_env2$herbicide), label = T, col = 1:4, lwd=3)
ordiellipse(allbutnats_nmds, allbutnats_envPA$BRCA, label = T, lwd=3, col = 1:2, choices = c(1,3), kind = "se")

# axes 1 and 3
ordiplot(allbutnats_nmds, type = "n", choices = c(1,3))
orditorp(allbutnats_nmds, display="species",col="red", air=0.01, choices = c(1,3))
orditorp(allbutnats_nmds, display="sites",col="grey50", air=0.01, choices = c(1,3))
ordiellipse(allbutnats_nmds, paste(allbutnats_env2$seedtrt, allbutnats_env2$herbicide), label = T, col = 1:4, lwd=3, choices = c(1,3))
ordisurf(allbutnats_nmds,allbutnats_env2$FEMI,col="forestgreen", choices = c(1,3)) # a little more linear
ordisurf(allbutnats_nmds,allbutnats_env2$ESCA,col="orchid", choices = c(1,3)) # a little more linear
ordisurf(allbutnats_nmds,allbutnats_env2$NEMA,col="chocolate3", choices = c(1,3))
ordiellipse(allbutnats_nmds, allbutnats_envPA$FEMI, label = T, lwd=3, col = 1:2, choices = c(1,3))
ordiellipse(allbutnats_nmds, allbutnats_envPA$ESCA, label = T, lwd=3, col = 5:6, choices = c(1,3), kind = "se")
ordiellipse(allbutnats_nmds, allbutnats_envPA$ESCA, label = T, lwd=3, col = 5:6, choices = c(1,3))
ordiellipse(allbutnats_nmds, allbutnats_envPA$herbicide, label = T, lwd=3, col = 3:4, choices = c(1,3), kind = "se")

# axes 2 and 3
ordiplot(allbutnats_nmds, type = "n", choices = c(2,3))
orditorp(allbutnats_nmds, display="species",col="red", air=0.01, choices = c(2,3))
orditorp(allbutnats_nmds, display="sites",col="grey50", air=0.01, choices = c(2,3))
ordiellipse(allbutnats_nmds, paste(allbutnats_env2$seedtrt, allbutnats_env2$herbicide), label = T, col = 1:4, lwd=3, choices = c(2,3))
ordisurf(allbutnats_nmds,allbutnats_env2$FEMI,col="forestgreen", choices = c(2,3)) # funny (or not) femi and esca capture in same space bc respond similar in models too
ordisurf(allbutnats_nmds,allbutnats_env2$ESCA,col="orchid", choices = c(2,3))
ordisurf(allbutnats_nmds,allbutnats_env2$NEMA,col="chocolate3", choices = c(2,3)) #linear!
ordiellipse(allbutnats_nmds, allbutnats_envPA$NEMA, label = T, lwd=3, choices = c(2,3), kind = "sd")
ordihull(allbutnats_nmds, allbutnats_envPA$seedtrt, label = T, col = 1:2,lwd=3, choices = c(2,3))
allbutnats_envPA$fulltrt <- with(allbutnats_envPA, paste(seedtrt, herbicide, ppt_trt, sep = "_"))
ordiellipse(allbutnats_nmds, allbutnats_envPA$fulltrt, label = T, col = 1:10,lwd=3, choices = c(2,3), kind = "se")

# relativized plot
ordiplot(allbutnats_nmds_rel, type = "n")
orditorp(allbutnats_nmds_rel, display="species",col="red", air=0.01)
orditorp(allbutnats_nmds_rel, display="sites",col="grey50", air=0.01)
ordiellipse(allbutnats_nmds_rel, allbutnats_envPA$NEMA, label = T, lwd=2)
ordiellipse(allbutnats_nmds_rel, allbutnats_envPA$FEMI, label = T, lwd=3, col = 2:3, kind = "se")
ordiellipse(allbutnats_nmds_rel, allbutnats_envPA$ESCA, label = T, lwd=3, col = 4:5, kind = "ehull")
ordisurf(allbutnats_nmds_rel,allbutnats_env2$NEMA,col="chocolate3") #linear!
ordisurf(allbutnats_nmds_rel,allbutnats_env2$FEMI,col="forestgreen")
ordisurf(allbutnats_nmds_rel,allbutnats_env2$ESCA,col="orchid")
ordisurf(allbutnats_nmds_rel,allbutnats_env2$BRCA,col="blue")
ordiellipse(allbutnats_nmds_rel, allbutnats_envPA$fulltrt, label = T, col = 1:10,lwd=1, kind = "se")
ordihull(allbutnats_nmds_rel, allbutnats_envPA$seedtrt, label = T, col = 1:10,lwd=1)

# permanova to see if native seeded, or herbicide explain comm diss at all
adonis2(allbutnats_bc ~ seedtrt + herbicide +ppt_trt + nut_trt, data = allbutnats_envPA, permutations = how(blocks = allbutnats_env2$block))
# for now I'm saying herbicide, ppt_trt, nut_trt, and FEMI explain community diss the most.. seeding treatment does to
adonis2(allbutnats_bc ~ seedtrt + herbicide + FEMI + ppt_trt + nut_trt, data = allbutnats_envPA, permutations = how(blocks = allbutnats_env2$block))
adonis2(allbutnats_bc ~ seedtrt * nut_trt + herbicide + FEMI + ppt_trt, data = allbutnats_envPA, permutations = how(blocks = allbutnats_env2$block))

# plot all but nats in ggplot so can see plot trtments better
nmds_allbutnats_gg <- cbind(allbutnats_env2, allbutnats_nmds$points)

ggplot(nmds_allbutnats_gg, aes(MDS1, MDS2, col = herbicide)) +
  geom_label(aes(label = block)) +
  facet_grid(ppt_trt~nut_trt)

# one more nmds using fxnl grps -- compare experimental subplots only 
mds_nats_df <- subset(natlong, seedtrt != "Unseeded", select = c(plot:mon, code4,pct_cover)) %>%
  # choose max cov within season
  group_by(plot, herbicide, seedtrt, code4) %>%
  summarise(pct_cover = max(pct_cover), .groups = "drop_last") %>%
  distinct() %>%
  # remove rare (present in fewer than 5% of plots)
  group_by(code4) %>%
  mutate(freq = length(pct_cover[pct_cover > 0])) %>%
  ungroup() %>%
  mutate(plots = length(unique(paste(plot, herbicide, seedtrt))),
         prop_pres = freq/plots) %>%
  subset(prop_pres >=0.05, select = -c(freq, plots, prop_pres)) %>%
  spread(code4, pct_cover, fill = 0) %>%
  #create unique ID code for each subplot
  mutate(natplotid = 1:nrow(.))

mds_nats_nmds <- metaMDS(subset(mds_nats_df, select = -c(plot:seedtrt, natplotid)), k = 2, trymax = 100)
plot(mds_nats_nmds)
mds_nats_results <- left_join(mds_nats_df[c("plot", "herbicide", "seedtrt")], trtkey) %>%
  cbind(mds_nats_nmds$points)

mds_nats_spp <- data.frame(mds_nats_nmds$species) %>%
  mutate(code4 = rownames(.)) %>%
  left_join(distinct(subset(spplist, select = -c(species, compost_synonym)))) %>%
  mutate(nativity = gsub("Unknown", "Exotic", nativity))

ggplot(mds_nats_results, aes(MDS1, MDS2)) +
  geom_label(aes(label = block)) + # nmds separates on axis 1 more by hillpos than anything else
  facet_grid(ppt_trt~nut_trt)

ggplot(mds_nats_results, aes(MDS1, MDS2)) +
  geom_label(aes(label = block)) + # nmds separates on axis 1 more by hillpos than anything else
  facet_grid(ppt_trt~nut_trt)

ggplot(data = mds_nats_results, aes(MDS1, MDS2)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  geom_point(aes(shape = paste(block >2, herbicide)), alpha = 0.9, size = 3) +
  geom_text(data = subset(mds_nats_spp, !code4 %in% nats), aes(MDS1, MDS2, color = paste(nativity,fxnl_grp), label = code4), size = 2.5) +
  geom_text(data = subset(mds_nats_spp, code4 %in% nats), aes(MDS1, MDS2, color = paste(nativity,fxnl_grp), label = code4), size = 3) +
  scale_shape_manual(name = "hillpos + herbicide", values = c(16,1, 17, 2), labels = c("lo, herbicide", "lo, control", "hi, herbicide", "hi, control")) +
  scale_color_manual(name = NULL, values = c('Exotic Forb' = "red", 'Native Forb' = "purple", 'Exotic Grass' = "forestgreen", 'Native Grass' = "dodgerblue", "Exotic N-fixer" = "orange")) +
  facet_grid(ppt_trt~nut_trt)
ggsave("/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/compost_natdiv_spp_mds.pdf")

# check on abundance of each nat per treatment
left_join(seededspp, trtkey) %>%
  gather(code4, cover, BRCA:TRCI) %>%
  group_by(herbicide, block > 2) %>%
  ggplot(aes(code4, cover, group = paste(herbicide, block >2), col = block > 2)) +
  geom_jitter(aes(shape = herbicide), width = 0.2,height = 0, alpha = 0.5, size = 2) +
  stat_summary(aes(shape = herbicide), alpha = 0.7) +
  facet_grid(ppt_trt~nut_trt) +
  coord_flip()

# plot by species to see cover change across treatments relative for each spp
left_join(seededspp, trtkey) %>%
  gather(code4, cover, BRCA:TRCI) %>%
  subset(code4 != "TRCI") %>%
  mutate(hillpos = ifelse(block > 2, "Upslope", "Downslope")) %>%
  #group_by(herbicide, block > 2) %>%
  ggplot(aes(fulltrt, cover, group = paste(fulltrt, herbicide), col = herbicide)) +
  #geom_jitter(aes(shape = herbicide), width = 0.2,height = 0, alpha = 0.5, size = 2) +
  stat_summary(aes(shape = nut_trt), position = position_dodge(width = 0.4), alpha = 0.9) +
  scale_shape_manual(name = "Amendment", values = c(C = 16, F = 15, N = 4), labels = c("Compost", "Fert", "None")) +
  #scale_shape_manual(name = "Subplot", values = c(21,22)) +
  #scale_color_brewer(name = "Amendment", palette = "Set2", labels = c("Compost", "Fert", "None")) +
  scale_color_brewer(name = "Subplot", palette = "Set2") +
  #scale_fill_discrete(name = "Precip trt") +
  facet_grid(code4~hillpos, scales = "free")




