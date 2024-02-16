# assess coverage of traits for nat recruitment experiment

# main data source:
# > hallett lab GH traits screening (Ashley, Lina, Carmen)
# > supplement with Julie's GH trait screening for ClimVar
# > then Loralee's (in ClimVar)
# > then Butterfield if needed



# -- SETUP -----
# load needed libraries not in source veg script
library(cowplot)

source("Native-analysis/native_prep_cover.R")

# path to google drive (Julie's folder) [can make this for more users w access via googledrive package later if wish.. but this less annoying for now]
# googlepath <- "/Users/scarlet/Library/CloudStorage/GoogleDrive-cawh3971@colorado.edu/Mi\ unidad/CA_Traits_Collaboration" 
# # list all csvs in traits collab
# gfiles <- list.files(googlepath, full.names = T, pattern = "csv")
# # read all csvs to list to inspect
# climvar_traits_list <- lapply(gfiles, read_csv)
# # set file names as list element names
# names(climvar_traits_list) <- gsub(paste0(googlepath, "/"), "", gfiles)

# ashley, lina traits [cleaned by as]
asla_traits <- read.csv("/Users/scarlet/Library/CloudStorage/Dropbox/USDA-compost/Data/PlantTraits/greenhouse_traits/GH_CleanedData/GreenhouseTraits.csv")  



# -- SCREEN SPP NEEDED v. TRAIT AVAILABLE ------
# calculate from observed cover for each species:
# 1. relative prop sp represents in entire nat recruit dataset (all 4 subplots per soil x ppt treat)
# 2. relative prop in seeded plots only
# 2. average cover in all plots (4 subplots per soil x ppt trt in 2021)

# calculate using count since will be modeling negbin distr
allplots_summary <- tidysp0 %>%
  # first determine totcov across all plots
  mutate(allcov = sum(count_cover)) %>%
  # by sp
  group_by(code4) %>%
  reframe(sp_relcov = round(sum(count_cover)/allcov, 6),
          sp_meancov = round(mean(count_cover),6)) %>%
  distinct() %>%
  arrange(desc(sp_relcov))
  # calculate rel 

seededplots_summary <- tidysp0_seeded %>%
  # first determine totcov across all plots
  mutate(allcov = sum(count_cover)) %>%
  # by sp
  group_by(code4) %>%
  reframe(sp_relcov = round(sum(count_cover)/allcov, 6),
          sp_meancov = round(mean(count_cover),6)) %>%
  distinct() %>%
  arrange(desc(sp_relcov))


# -- species-screen ashley, line & co trait data ----
sort(unique(asla_traits$ID)) # drop purchased v field distinction to pair, need to 4letter some codes
# try to crosswalk SciName in spp key
asla_spp <- distinct(asla_traits, Taxon,ID) %>%
  mutate(ID = trimws(gsub("[p|f]$|pi$", "", ID)),
         Taxon = trimws(Taxon),
         # make crosswalk for code4
         xwalk_code4 = ifelse(nchar(ID)==4, ID, ID),
  # hardcode others
 xwalk_code4 = recode(xwalk_code4, `Agoseris` = spplist$code4[grepl("Agos", spplist$genus)], 
              `Leontodon` = spplist$code4[grepl("Leon.* sp", spplist$species)], 
              `AVEBAR` = spplist$code4[grepl("Av.* ba", spplist$species)], 
                   # only one crassula in spp list so code that
                   `Crassula` = spplist$code4[grepl("Crass", spplist$genus)],
              # assign Sanicula one of the SANI spp then remove species info
              `Sanicula`= spplist$code4[grepl("Sani", spplist$genus)][1],
              # fectuca myuros = vulpia myuros
              `FEMY` = "VUMY")) %>%
  # distinct to remove multiple NA rows for things that didn't pair
  distinct() %>%
  left_join(spplist, by = c("xwalk_code4" = "code4"))
# ^ some did not pair bc difference code4 or spp not in Compost comp


# pull things that still need pairing and pair by code
asla_spp_unmatched <- subset(asla_spp, is.na(unknown), select = c(Taxon:xwalk_code4)) %>%
  # try pairing on sp name
  left_join(spplist, by = c("Taxon" = "species")) %>%
  # create species col so can rbind
  mutate(species = Taxon)
# ^only bromus diandrus has match, all else not present in Compost comp (e.g., some of these spp match ClimVar spp)  

# put things that have a pair back in, then drop Taxon rows that don't
asla_sppxwalk <- subset(asla_spp, !Taxon %in% with(asla_spp_unmatched, Taxon[!is.na(unknown)])) %>%
  rbind(subset(asla_spp_unmatched, !is.na(unknown), select = names(.))) %>%
  arrange(species) %>%
  # for simplicity, drop all sppinfo cols and rejoin on 'species' so cols are in order, code4 is corred
  subset(select = c(Taxon, ID, species)) %>%
  left_join(spplist) %>%
  # create alternate code for sanicula
  mutate(code4_alternate = ifelse(grepl("Sanic", species), spplist$code4[grepl("Sani", spplist$genus)][2], 
                                  # and for leontodon (taraxacoides is the only sp in spplist)
                                  ifelse(grepl("Leon", species), spplist$code4[grepl("Leon.* tara", spplist$species)], NA)),
         # make cols to note if spp is in main composition and if in native recruit composition
         in_maincomp = code4 %in% coverlong$code4 | code4_alternate %in% coverlong$code4,
         # in any of the 4 sub-sub-subplots in 2021
         in_natambcomp = code4 %in% natambsp_cov$code4 | code4_alternate %in% natambsp_cov$code4,
         # in seeded sub-sub-subplots
         in_seededcomp = code4 %in% tidysp0_seeded$code4 | code4_alternate %in% tidysp0_seeded$code4)

# assess coverage
table(unique(natambsp_cov$code4) %in% unique(c(asla_sppxwalk$code4,asla_sppxwalk$code4_alternate))) # 53 present, 19 not
# assess coverage
table(unique(tidysp0_seeded$code4) %in% unique(c(asla_sppxwalk$code4,asla_sppxwalk$code4_alternate))) # 45 present, 11 not (probably unk and native forbs?)
# check ones that aren't
no_aslatraits_code4 <- sort(unique(natambsp_cov$code4[!natambsp_cov$code4 %in% unique(c(asla_sppxwalk$code4,asla_sppxwalk$code4_alternate))]))
no_asltatraits_code4_seeded <- sort(unique(tidysp0_seeded$code4[!tidysp0_seeded$code4 %in% unique(c(asla_sppxwalk$code4,asla_sppxwalk$code4_alternate))]))
View(subset(spplist, code4 %in% no_aslatraits_code4))
distinct(spplist[spplist$code4 %in% no_asltatraits_code4_seeded,], code4, species)
# ^I don't think Julie's 2 trait datasets or Brad's will have any of these
# > can use Galium aparine as approximation of Galium parisiesce (if collected in compost pasture, most likely parisiense -- we called in aparine for a few years in ClimVar)
# > can average hypo glabra and radi for aster spp?

# assess how much cover these spp comprise
notraits_propcov <- tidysp0_seeded %>%
  group_by(nut_trt, ppt_trt, herbicide, seedtrt) %>%
  mutate(totcov = sum(count_cover)) %>%
  subset(code4 %in% no_aslatraits_code4) %>%
  group_by(code4, totcov, .add = T) %>%
  reframe(spcov = sum(count_cover),
          propcov = round(spcov/totcov,4)) %>%
  distinct() %>%
  arrange(desc(propcov))
# ^ the only species that comprises non-neglible cover is Galium paris., 1-4%. Can use Galium aparine as approx. All else is < .01 proportionally for the treatment combo.

# assign GAAP to GAPA row
asla_sppxwalk$code4_alternate[grepl("GAAP", asla_sppxwalk$code4)] <- spplist$code4[grepl("Galium par", spplist$species)]


# -- AGGREGATE TRAIT DATA -----
# take means, sd, and n per trait per sp. keep purchased from field separated to review whether to pool or only use field (i.e., how different are they?)

# check that traits are all numeric
str(asla_traits) # date is character, but may string flatten so ok
# Root.volume..cm3., leaf mass cols are character. see why
with(asla_traits, unique(Root.volume..cm3.[grepl("[a-z]| ", Root.volume..cm3.)])) # "not scanned"
with(asla_traits, unique(Root.dry.biomass..g.[grepl("[a-z]| ", Root.dry.biomass..g.)])) # "lost", "roots lost"
with(asla_traits, unique(Dry.leaf.mass..g.[grepl("[a-z]| ", Dry.leaf.mass..g.)])) # "No sample"
# ^if convert these to numeric will make NA
sort(unique(asla_traits$Root.volume..cm3.))
asla_traits_aggregated <- group_by(asla_traits, Taxon, ID) %>%
  # convert and order dates so listed chronologically
  mutate(Date.harvest = as.character(Date.harvest)) %>%
  arrange(Taxon, Date.harvest) %>%
  mutate(Date.harvest = str_flatten_comma(unique(Date.harvest))) %>%
  ungroup() %>%
  # clean up .. in trait names before gather
  rename_all(function(x) gsub( "[.]{2}", ".", x)) %>%
  # remove period as last character
  rename_all(function(x) gsub( "[.]$", "", x)) %>%
  # gather traits to take average
  gather(trait, value, names(.)[!grepl("Tax|ID$|Rep|Date", names(.))],) %>%
  # convert to numeric (will throw warning)
  mutate(value = as.numeric(value)) %>% # warning is for notes converted to NA
  # take means, sd, n
  grouped_df(names(.)[!grepl("val|Rep", names(.))]) %>%
  reframe(avg = mean(value, na.rm = T),
          med = median(value, na.rm = T),
          sd = sd(value, na.rm = T),
          reps = sum(!is.na(value))
          ) %>%
  distinct() %>%
  # note which species are both field collected and purchased
  group_by(Taxon) %>%
  mutate(spnID = length(unique(ID))) %>%
  ungroup()

# review differences
subset(asla_traits_aggregated, spnID > 1) %>%
  ggplot(aes(Taxon, avg, col = grepl("p$|pi$", ID), group = ID)) +
  geom_jitter(width = 0.2, height = 0) +
  facet_wrap(~trait, scales = "free") +
  scale_color_manual(name = "Seed source", values = c("seagreen3", "mediumpurple2"), labels = c(`FALSE` = "Field", `TRUE`= "Purchased")) +
  coord_flip() +
  theme(legend.position = c(0.97, 0.03),
        legend.justification = c("right", "bottom"))
# ^ the species that have greatest diffs (just eyeballing) are not in natambsp; use field values for those that are (e.g., TRHI, BRHO)


# -- MAKE NATAMB TRAIT DATA DF -----
# choose just the spp we have coverage for.. can either source or write out

natamb_traits <- left_join(asla_traits_aggregated, asla_sppxwalk) %>%
  # drop purchased if there is a field version
  filter(!(grepl("[p|pi]$", ID) &spnID > 1))
# check codes present
data.frame(distinct(natamb_traits, Taxon, ID)) # looks good

# winnow to just those in the native recruitment - 2021 ambient sub-sub-subplots
natamb_traits <- subset(natamb_traits, (code4 %in% tidysp0$code4) | (code4_alternate %in% tidysp0$code4))

# need species in rows and traits in columns
natamb_avgtraits_wide <- subset(natamb_traits, select = c(code4, code4_alternate, in_maincomp:in_seededcomp, trait, avg)) %>%
  spread(trait, avg) %>% 
  # pair descriptive traits
  left_join(subset(spplist, select = c(code4, Category,Family, Duration,Growth_Habit, fxnl_grp, nativity))) %>%
  # condense duration
  mutate(duration.abbr = ifelse(grepl(", Pere|Unk", Duration), "Multi-year", # the unknown is an aster, which can be annual - peren
                                ifelse(grepl("Bien", Duration), "Biennial", Duration)))



# clean up env (remove anything not needed in source) 


