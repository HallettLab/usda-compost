# ------------
colnames(traits.means.comp)
envfit(nmds.comp.trait ~ Height..cm._mean + 
         Fresh.leaf.mass..g._mean + 
         Dry.leaf.mass..g._mean +
         LDMC_mean +
         Leaf.Area..cm2._mean +
         SLA..cm2.g._mean +
         Shoot.dry.biomass..g._mean +
         Root.dry.biomass..g._mean +
         Total.biomass..g._mean +
         RMF_mean +
         Root.volume..cm3._mean +
         Root.density..g.cm3._mean +
         Coarse.root.diameter..mm._mean +
         Length..mm._mean +
         Fine.root.length..mm._mean +
         Coarse.root.length..mm._mean +
         Coarse.root.specific.length..cm.g._mean +
         Fine.root.specific.length..cm.g._mean +
         Proportion.fine.roots_mean, traits.means.comp, display = "sp") # ONLY Height, Leaf area, and RMF significant



# is it different when you only do certain traits (above- vs. below-ground)? or can you just throw everything in there and it will tell you what is significant and you can subset from there?

envfit(nmds.comp.trait ~ Height..cm._mean + 
         Leaf.Area..cm2._mean +
         RMF_mean, traits.means.comp, display = "sp") # Looks consistent with mega-model

# try different variations to see what happens (only leaf stuff)
envfit(nmds.comp.trait ~ Height..cm._mean + 
         Fresh.leaf.mass..g._mean + 
         Dry.leaf.mass..g._mean +
         LDMC_mean +
         Leaf.Area..cm2._mean +
         SLA..cm2.g._mean +
         Shoot.dry.biomass..g._mean, traits.means.comp, display = "sp")

# (only root stuff)
envfit(nmds.comp.trait ~ Root.dry.biomass..g._mean +
         Total.biomass..g._mean +
         RMF_mean +
         Root.volume..cm3._mean +
         Root.density..g.cm3._mean +
         Coarse.root.diameter..mm._mean +
         Length..mm._mean +
         Fine.root.length..mm._mean +
         Coarse.root.length..mm._mean +
         Coarse.root.specific.length..cm.g._mean +
         Fine.root.specific.length..cm.g._mean +
         Proportion.fine.roots_mean, traits.means.comp, display = "sp")

vectors.nmds.comp.trait <- envfit(nmds.comp.trait ~ Height..cm._mean + 
                                    Fresh.leaf.mass..g._mean + 
                                    Dry.leaf.mass..g._mean +
                                    LDMC_mean +
                                    Leaf.Area..cm2._mean +
                                    SLA..cm2.g._mean +
                                    Shoot.dry.biomass..g._mean +
                                    Root.dry.biomass..g._mean +
                                    Total.biomass..g._mean +
                                    RMF_mean +
                                    Root.volume..cm3._mean +
                                    Root.density..g.cm3._mean +
                                    Coarse.root.diameter..mm._mean +
                                    Length..mm._mean +
                                    Fine.root.length..mm._mean +
                                    Coarse.root.length..mm._mean +
                                    Coarse.root.specific.length..cm.g._mean +
                                    Fine.root.specific.length..cm.g._mean +
                                    Proportion.fine.roots_mean, traits.means.comp, display = "sp")

ggplot(nmds_dat_trait, aes(NMDS1, NMDS2,  fill = nut_trt, color = nut_trt)) +
  geom_point(shape = 21, size = 4, stroke = 0.5) + 
  stat_ellipse(type = "norm", size = 1, alpha = 0.5) + # Ellipses for each group
  labs(title = "Grouped by Nutrient Treatment") +
  custom_colors +
  scale_color_manual(
    values = c(
      "control" = rgb(220, 220, 220, maxColorValue = 255),      # Lighter Gray
      "+N fertilizer" = rgb(140, 40, 150, maxColorValue = 255), # Vibrant Purple
      "+compost" = rgb(220, 180, 60, maxColorValue = 255),       # Bright Golden Yellow
      "wet" = rgb(0, 0, 255, maxColorValue = 255),
      "dry" = rgb(183, 65, 14, maxColorValue = 255)
    )
  ) +
  geom_segment(data = as.data.frame(vectors.nmds.comp.trait$vectors$arrows),
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", size = 1)

# ---------------------------

# subset nmds to match species names
nmds.comp.trait <- metaMDS(matrix.bray %>%
                             dplyr::select(intersect.names), 
                           k=2, trymax = 25)
nmds.comp.trait # pretty not great stress
stressplot(nmds.comp.trait)
plot(nmds.comp.trait, type="t")

nmds1 <- as.data.frame(scores(nmds.comp.trait, choices=c(1), display=c("sites")))
nmds2 <- as.data.frame(scores(nmds.comp.trait, choices=c(2), display=c("sites")))

nmds_dat_trait <- cbind(nmds1, nmds2) %>%
  as.data.frame() %>%  # Ensure it's a data frame
  rownames_to_column(var = "site_id") %>%  # Move row names into a new column
  separate(site_id, into = c("nut_trt", "ppt_trt", "yr", "grazing_hist"), sep = "_")


# subset matrix and trait data to match names
traits.means.comp <- traits.means %>%
  filter(ID %in% intersect.names) %>%
  rename(
    Height = Height..cm._mean,
    FreshLeafMass = Fresh.leaf.mass..g._mean,
    DryLeafMass = Dry.leaf.mass..g._mean,
    LDMC = LDMC_mean, # Leaf dry matter content (LDMC, the ratio of leaf dry mass to fresh mass)
    LeafArea = Leaf.Area..cm2._mean,
    SLA = SLA..cm2.g._mean,
    ShootDryBiomass = Shoot.dry.biomass..g._mean,
    RootDryBiomass = Root.dry.biomass..g._mean,
    TotalBiomass = Total.biomass..g._mean,
    RMF = RMF_mean, # Root mass fraction (RMF) is a plant trait that measures the proportion of a plant's dry mass that is in its roots
    RootVol = Root.volume..cm3._mean,
    RootDensity = Root.density..g.cm3._mean,
    CoarseRootDiameter = Coarse.root.diameter..mm._mean,
    RootLength = Length..mm._mean,
    FineRootLength = Fine.root.length..mm._mean,
    CoarseRootLength = Coarse.root.length..mm._mean,
    CoarseRootSpecLength = Coarse.root.specific.length..cm.g._mean,
    FineRootSpecLength = Fine.root.specific.length..cm.g._mean,
    PropFineRoots = Proportion.fine.roots_mean
    
  )
# traits.means.comp2 <- traits.means.comp %>%
# dplyr::select(-ID)
rownames(traits.means.comp2) <- traits.means.comp$ID

# trait vectors onto the NMDS
# envfit(nmds.comp.trait, traits.means.comp, permutations = 999)
# example from stack overflow: envfit(ord ~ var1 + var2 + var3, dune.spec, display="sp"), 'ord' is the nmds model, 'var1' is the trait1 i think?

colnames(traits.means.comp)
vectors.nmds.comp.trait <- envfit(nmds.comp.trait ~ Height + 
                                    LeafArea +
                                    RMF, traits.means.comp, display = "sp") # ONLY Height, Leaf area, and RMF significant

# is it different when you only do certain traits (above- vs. below-ground)? or can you just throw everything in there and it will tell you what is significant and you can subset from there?

# Looks like different subsets of traits put in are the same, so I shall proceed with the mega models and verify with Lauren

trait_vectors <- as.data.frame(vectors.nmds.comp.trait$vectors$arrows)

ggplot(nmds_dat_trait, aes(NMDS1, NMDS2,  fill = nut_trt, color = nut_trt)) +
  geom_point(shape = 21, size = 4, stroke = 0.5) + 
  # stat_ellipse(type = "norm", size = 1, alpha = 0.5) + # Ellipses for each group
  labs(title = "Grouped by Nutrient Treatment") +
  custom_colors +
  scale_color_manual(
    values = c(
      "control" = rgb(220, 220, 220, maxColorValue = 255),      # Lighter Gray
      "+N fertilizer" = rgb(140, 40, 150, maxColorValue = 255), # Vibrant Purple
      "+compost" = rgb(220, 180, 60, maxColorValue = 255),       # Bright Golden Yellow
      "wet" = rgb(0, 0, 255, maxColorValue = 255),
      "dry" = rgb(183, 65, 14, maxColorValue = 255)
    )
  ) +
  geom_segment(
    data = as.data.frame(vectors.nmds.comp.trait$vectors$arrows), 
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
    inherit.aes = FALSE,  # 🚀 Prevents it from looking for 'nut_trt'
    arrow = arrow(length = unit(0.2, "cm")), 
    color = "red", 
    size = 1
  ) +
  geom_text(
    data = as.data.frame(vectors.nmds.comp.trait$vectors$arrows), 
    aes(x = NMDS1, y = NMDS2, label = rownames(vectors.nmds.comp.trait$vectors$arrows)), 
    inherit.aes = FALSE,  # Avoids unwanted inheritance
    hjust = -0.2, 
    vjust = -0.2, 
    color = "red",
    size = 5
  )


# change trait names to be more read-able for figures -----
# colnames(traits.means)
# traits.means %>%
#   rename(
#     Height = Height..cm._mean,
#     FreshLeafMass = Fresh.leaf.mass..g._mean,
#     DryLeafMass = Dry.leaf.mass..g._mean,
#     LDMC = LDMC_mean, # Leaf dry matter content (LDMC, the ratio of leaf dry mass to fresh mass)
#     LeafArea = Leaf.Area..cm2._mean,
#     SLA = SLA..cm2.g._mean,
#     ShootDryBiomass = Shoot.dry.biomass..g._mean,
#     RootDryBiomass = Root.dry.biomass..g._mean,
#     TotalBiomass = Total.biomass..g._mean,
#     RMF = RMF_mean, # Root mass fraction (RMF) is a plant trait that measures the proportion of a plant's dry mass that is in its roots
#     RootVol = Root.volume..cm3._mean,
#     RootDensity = Root.density..g.cm3._mean,
#     CoarseRootDiameter = Coarse.root.diameter..mm._mean,
#     RootLength = Length..mm._mean,
#     FineRootLength = Fine.root.length..mm._mean,
#     CoarseRootLength = Coarse.root.length..mm._mean,
#     CoarseRootSpecLength = Coarse.root.specific.length..cm.g._mean,
#     FineRootSpecLength = Fine.root.specific.length..cm.g._mean,
#     PropFineRoots = Proportion.fine.roots_mean
#     
#   )
