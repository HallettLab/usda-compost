# native multivariate analyses
# author(s): ctw
# questions: caitlin.t.white@colorado.edu


# -- SETUP -----
# load needed libraries
library(tidyverse)
library(vegan)

# source prepped data
# > this should set default settings for datpath, plot theme, no factors by default
source("Native-analysis/native_prep_cover.R", )



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
mds_nats_nmds
plot(mds_nats_nmds)

# create bc matrix of all spp, only seededplots
allspp_bc <- vegdist(subset(mds_nats_df, select = -c(plot:seedtrt, natplotid)), method = "bray")
seeded_env <- subset(mds_nats_df, select = c(plot:seedtrt, natplotid)) %>%
  left_join(distinct(natlong[natlong$seedtrt != "Unseeded",], plot, block, ppt_trt, nut_trt, herbicide))
# check signif of treatments
adonis2(allspp_bc ~ herbicide * ppt_trt *nut_trt, data = seeded_env, permutations = how(blocks = seeded_env$block))

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  seeded_env$block 
# Permutation: free
# Number of permutations: 199

adonis2(formula = allspp_bc ~ herbicide * nut_trt * ppt_trt, data = seeded_env, permutations = how(blocks = seeded_env$block))
adonis2(formula = allspp_bc ~nut_trt * ppt_trt *herbicide, data = seeded_env, permutations = how(blocks = seeded_env$block))
# > doesn't matter if change order of terms
# Df SumOfSqs      R2      F Pr(>F)   
# herbicide                  1   0.8831 0.04454 3.2606  0.005 **
# nut_trt                    2   1.4821 0.07475 2.7361  0.005 **
# ppt_trt                    2   0.6851 0.03455 1.2648  0.035 * 
# herbicide:nut_trt          2   0.6313 0.03184 1.1654  0.040 * 
# herbicide:ppt_trt          2   0.5236 0.02641 0.9667  0.125   
# nut_trt:ppt_trt            4   0.9785 0.04935 0.9032  0.170   
# herbicide:nut_trt:ppt_trt  4   0.5601 0.02825 0.5170  0.870   
# Residual                  52  14.0833 0.71031                 
# Total                     69  19.8270 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(formula = allspp_bc ~ herbicide + nut_trt * ppt_trt, data = seeded_env, permutations = how(blocks = seeded_env$block))

# Df SumOfSqs      R2      F Pr(>F)   
# herbicide          1   0.8831 0.04454 3.3910  0.005 **
# nut_trt            2   1.4821 0.07475 2.8456  0.005 **
# ppt_trt            2   0.6851 0.03455 1.3154  0.025 * 
# herbicide:nut_trt  2   0.6313 0.03184 1.2121  0.035 * 
# Residual          62  16.1455 0.81432                 
# Total             69  19.8270 1.0000

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

# make hulls (quickplot.. spider will be easier on eyes to show maybe)
hulldat <- data.frame()
for(h in unique(mds_nats_results$herbicide)){
  for(f in unique(mds_nats_results$fulltrt)){
    grphull <- mds_nats_results[mds_nats_results$fulltrt == f & mds_nats_results$herbicide == h, ][chull(mds_nats_results[mds_nats_results$fulltrt == f & mds_nats_results$herbicide == h, 
                                                                                                                          c("MDS1", "MDS2")]), ]
    # add to df
    hulldat <- rbind(hulldat, grphull)  
  }
}

# kns says show species in one panel, treatments in other

# make pretty names for spp nmds
mds_nats_spp <- mutate(mds_nats_spp, 
                       pretty_name = ifelse(code4 %in% nats, paste(paste0(substr(genus,1,1), "."), epithet),paste(nativity, fxnl_grp)),
                       pretty_name = gsub("Exotic", "Non-native", pretty_name),
                       pretty_name = gsub("Grass", "grasses", pretty_name),
                       pretty_name = gsub("Forb", "forbs", pretty_name),
                       pretty_name = gsub("N-fixer", "N-fixers", pretty_name),
                       pretty_name = gsub("Native fo", "Background native fo", pretty_name),
                       pretty_name = factor(pretty_name, 
                                            levels = c("Non-native forbs", "Non-native grasses", "Non-native N-fixers",
                                                       "Background native forbs", "E. californica", "N. maculata", 
                                                       "B. carinatus", "F. microstachys"))
                       )
unique(mds_nats_spp$pretty_name)
#mds_nats_spp$pretty_name <- with(mds_)

mds_spp <- ggplot(data = mutate(mds_nats_results, 
                                ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W"))
                                #ppt_trt = relevel(ppt_trt, "D")
                                ), aes(MDS1, MDS2)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  #geom_polygon(data=hulldat, aes(fill=ppt_trt,group=ppt_trt, col = ppt_trt), alpha=0.35) +
  #geom_point(aes(shape = block<3), col = "grey30", alpha = 0.7, size = 3) +
  #geom_point(data = subset(mds_nats_spp, !code4 %in% nats), aes(MDS1, MDS2, color = pretty_name)) +
    geom_point(data = mds_nats_spp, aes(MDS1, MDS2, fill = pretty_name), pch = 21, col = "grey30", alpha = 0.85, size = 3) +
  #geom_text(data = subset(mds_nats_spp, !code4 %in% nats), aes(MDS1, MDS2, color = paste(nativity,fxnl_grp), label = code4), size = 3) +
  ggrepel::geom_text_repel(data = subset(mds_nats_spp, code4 %in% nats), 
                           aes(MDS1, MDS2, color = pretty_name, label = pretty_name), show.legend = FALSE, size = 4) +
  scale_color_manual(name = NULL, values = plant_cols2) +
    scale_fill_manual(name = NULL, values = plant_cols2) +
  #theme( legend.text = element_text(size = 10))
  #scale_color_manual(name = NULL, values = c("orchid", "seagreen2", "red", "purple4", "seagreen4")) +
  theme(legend.position = c(0.01,0.99),
        legend.text = element_text(size = 10),
        legend.justification = c("left", "top"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.1, "pt"))
                       #c('Exotic Forb' = "red", 'Native Forb' = "purple", 'Exotic Grass' = "forestgreen", 'Native Grass' = "dodgerblue", "Exotic N-fixer" = "orange"))

# plot hulls
mds_plots <- ggplot(data = mutate(mds_nats_results, 
                                  herbicide = factor(herbicide, levels = c("Non-herbicided", "Herbicided")),
              ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W")),
              #ppt_trt = relevel(ppt_trt, "D"),
              nut_trt = factor(nut_trt, levels = c("C", "F", "N"), labels = c("Compost", "Fertilizer", "No amendment"))),
       aes(MDS1, MDS2, col = ppt_trt, fill = ppt_trt)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  geom_polygon(data=mutate(hulldat, 
                           herbicide = factor(herbicide, levels = c("Non-herbicided", "Herbicided")),
                           ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W")),
                           nut_trt = factor(nut_trt, levels = c("C", "F", "N"), labels = c("Compost", "Fertilizer", "No amendment"))
                           ), 
               aes(fill=ppt_trt,group=ppt_trt, col = ppt_trt), alpha=0.35) +
  geom_point(aes(shape = block<3), col = "grey30", alpha = 0.7, size = 3) + #col = "grey30",
  geom_point(aes(shape = block<3), alpha = 0.7, size = 3) + #col = "grey30",
  #geom_point(alpha = 0.7, size = 3) +
  scale_shape_manual(name = "Hillslope", values = c(24, 25), labels = c("Up", "Down"), guides(color = "none")) +
  #scale_color_manual(name = NULL, values = c('Exotic Forb' = "red", 'Native Forb' = "purple", 'Exotic Grass' = "forestgreen", 'Native Grass' = "dodgerblue", "Exotic N-fixer" = "orange")) +
  scale_color_manual(name = "Precipitation\ntreatment", 
                     labels = c("D" = "Drought", "XC" = "Ambient", "W" = "Wet"),
                     values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  scale_fill_manual(name = "Precipitation\ntreatment", 
                    labels = c("D" = "Drought", "XC" = "Ambient", "W" = "Wet"), 
                    values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  #guides( = "none") +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 10)) +
  facet_grid(herbicide~nut_trt)

cowplot::plot_grid(mds_plots, mds_spp,
                   nrow = 1,
                   rel_widths = c(1.1,0.9))

ggsave(paste0(figpath,"natfxnl_seededspp_nmds2.pdf"), units = "in", width = 10.2, height = 5)


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




