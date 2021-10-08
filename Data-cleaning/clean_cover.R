# clean compost cover datasets
# author(s): ctw (caitlin.t.white@colorado.edu)
# created: 2019-04-28 (script will be modified over project lifetime)

# script purpose:
# read in master spp list and treatment lookup table
# iterate through all compost cover datasets to:
# change species names to codes, change "T" abundance values to 0.01
# build 2 types of cover data tables:
## 1) long form (tidy)
## 2) wide form (rows = site x date, cols = species)
## > for wide: transpose species to columns, sites to rows, fill in blank abundance values with 0s
## > for long: append all species descriptive columns in species key
## > for long and wide: add date sampled, recorder, notes, and plot treatment data
# write out clean cover datasets

# notes:
# re-run whenever unknowns are identified and changed in the entered cover dataset
# 2020 adjustment: process and append spring 2020 cover separately (Nikolai collected for us bc of COVID)
# > more unknowns to deal with and potentially bin 


# *** IMPORTANT TO DO BEFORE RUNNING THIS SCRIPT*** ##
# -- >> Run build_spplist.R before running the clean_cover script! << ---
# If new species have been entered, or to be sure, re-create the Compost species list via build_spplist.R
# Species list data appends to the compiled cover data, so important the species list dataset is complete before running this script




# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed
library(tidyverse) # for dplyr, tidyr, and stringr
library(lubridate)
options(stringsAsFactors = F)
na_vals <- c(" ", "", NA, "NA")
theme_set(theme_bw())

# set path to compost data (main level)
# > varies by user
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# list files in entered data folder
vegfiles <- list.files(paste0(datpath, "Cover/Cover_EnteredData"), full.names = T, pattern = "_Cover_", ignore.case = T)

# read in master spp list (lookup table)
spplist <- read.csv(paste0(datpath, "Compost_SppList.csv"), na.strings = na_vals, strip.white = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)




# -- TIDY AND TRANSPOSE ENTERED COVER DATA -----
# loop through cover data to read in, transpose, and append to master cover dataset
# initiate master data frames
cover_master_long <- data.frame()
cover_master_wide <- data.frame()
cover_list <- list()

# loop iterates through each cover dataset and adds new spp to master spp set
# > process 2020 separately
for(i in vegfiles){
  # read in dataset
  vegdat <- read.csv(i, na.strings = na_vals, header = F, blank.lines.skip = T, strip.white = T)
  print(paste("Transposing and tidying", str_extract(i, "(?<=_)[:alnum:]+(?=.csv)"), "composition dataset"))
  
  # store rawdat in list for post-compile checks
  cover_list[[which(i == vegfiles)]] <- vegdat
  names(cover_list)[which(i == vegfiles)] <- str_extract(i, "(?<=_)[:alnum:]+(?=.csv)")
  
  # remove blank rows
  vegdat <- vegdat[!apply(vegdat, 1, function(x) all(is.na(x))),]
  # remove common name column if there
  vegdat <- vegdat[,!sapply(vegdat, function(x) any(grep("^Common name", x, ignore.case = T)))]
  # remove Grasses and Forbs section headers if present
  vegdat <- vegdat[!apply(vegdat,1, function(x) all(x %in% c("Grasses", "Forbs", NA))), ]
  # id where plot and cover data start
  plotpos <- grep("plot", vegdat[,1], ignore.case = T)
  
  
  # -- CREATE PLOT-DATE DATA TABLE -----
  # pull out recorder, sample date, and notes into separate data frame
  notes <- vegdat[1:plotpos,]
  notes <- data.frame(t(notes)) #transpose, change matrix class to data frame
  colnames(notes) <- casefold(gsub(":", "", notes[1,])) # set row 1 as column names, make lower case and remove colon
  notes <- notes[-1,] # remove row 1
  rownames(notes) <- NULL # rename rownames in numeric sequence
  
  # change date from character to date format and add year column
  notes$date <- as.Date(notes$date, format = "%m/%d/%y")
  notes$yr <- as.numeric(substr(as.character(notes$date), 1,4))
  
  # prep plot col depending on which ID used
  # if plot col is number id, convert to numeric
  if(all(!grepl("[[:alpha:]]+", notes$plot))){
    notes$plot <- as.numeric(notes$plot)
  }
  # if plot ID is block-nut-ppt ID, rename col to plotid to join correctly with trtkey
  # > 2019 we use block-nut-ppt, 2020 use plot #
  if(any(grepl("[[:alpha:]]", notes$plot))){
    notes <- rename(notes, plotid = plot)
  } 
  
  # join plot treatment info
  notes <- left_join(notes, trtkey) #left_join preserves order of sampling, merge alphabetizes plots
  #reorder cols
  notes <- notes[c("plot", "plotid", "fulltrt", "block", "nut_trt", "ppt_trt", "yr", "date", "recorder", "notes")]
  
  
  # -- PREP ABUNDANCE DATA FOR TIDYING AND TRANSPOSING -----
  # remove notes df rows from vegdat
  vegdat <- vegdat[plotpos:nrow(vegdat),] 
  # id row where species abundance data starts
  litpos <- grep("depth", vegdat[,1])
  # remove any species rows that don't have any entries for abundance value
  allNAs <- apply(vegdat[, 2:ncol(vegdat)], 1, function(x) all(is.na(x)))
  # set rows up to litter depth as FALSE to preserve them (in case 0s not entered and rock or bare never encountered)
  allNAs[1:litpos] <- FALSE 
  # allNAs <- allNAs[!grepl("NA",names(allNAs))] # remove NAs created by ignoring rows 1 through litter depth
  vegdat <- vegdat[!allNAs,] 
  # rename rownames in sequential order to preserve correct order as in datasheets
  rownames(vegdat) <- 1:nrow(vegdat)
  
  
  # replace entered species name with synonym if present
  synonyms <- vegdat$V1[vegdat$V1 %in% spplist$species[!is.na(spplist$compost_synonym)]]
  if(length(synonyms)>0){
    for(s in synonyms){
      vegdat$V1[vegdat$V1 == s] <- spplist$compost_synonym[spplist$species == s]
    }
  }
  # join species codes and descriptive info
  vegdat <- right_join(spplist[c("species", "code4")], vegdat, by = c("species" = "V1"))
  vegdat <- distinct(vegdat)
  # add codes for non-species cells (e.g. percent grass, percent bare/other)
  vegdat$code4[grepl("%", vegdat$species)] <- with(vegdat[grepl("%", vegdat$species),], 
                                                   paste0("pct_", casefold(str_extract(species, "[A-Z][a-z]+"))))
  
  # ID pct cover column headers as plot or plotid depending on whether letter present
  ## track new litpos and plotpos because something apparently has changed in 2021 (CTW hates new package versions)
  plotpos <- grep("plot", vegdat[,1], ignore.case = T)
  litpos <- grep("depth", vegdat[,1])
  vegdat$code4[plotpos] <- ifelse(grepl("[[:alpha:]]", vegdat$V3[plotpos]), "plotid", "plot")
  vegdat$code4[litpos] <- "litter_depth_cm"
  # check to be sure all rows have a code4 value before proceeding
  stopifnot(all(!is.na(vegdat$code4)))
  
  vegdat <- vegdat[,2:ncol(vegdat)] # remove species col (using codes for headers)
  vegdat <- data.frame(t(vegdat)) # transpose, change matrix to data frame
  colnames(vegdat) <- vegdat[1,]
  vegdat <- vegdat[-1,]
  rownames(vegdat) <- NULL
  
  #clean up trace values and empty cells for non-species cells (e.g. percent green, percent bare, litter depth)
  vegdat[vegdat == "T"] <- 0.01
  # make sure plot and pct cols are first (got reorganized in 2021 with some type of R/package update)
  plotnames <- names(vegdat)[!grepl("[A-Z0-9]{4}", names(vegdat))]
  sppnames <- names(vegdat)[!names(vegdat) %in% plotnames]
  # reorganize
  vegdat <- cbind(vegdat[plotnames],vegdat[sppnames])
  vegdat[plotnames] <- sapply(vegdat[plotnames], function(x) ifelse(is.na(x)|x == "n/a",0,x))
  # make all cols except plot numeric [except if plot is the numeric ID and not fulltrt]
  vegdat[,!grepl("plot", names(vegdat), ignore.case = T)] <- sapply(vegdat[,!grepl("plot", names(vegdat), ignore.case = T)], as.numeric)
  # if plot numeric ID, also make numeric so it joins with trtkey
  if(!all(grepl("[[:alpha:]]", vegdat[,grepl("plot", names(vegdat), ignore.case = T)]))){
    vegdat[,grepl("plot", names(vegdat), ignore.case = T)] <- as.numeric(vegdat[,grepl("plot", names(vegdat), ignore.case = T)])
  }
  # reorder species cols alphabetically
  vegdat <- vegdat[c(plotnames, sort(sppnames))]
  
  # join notes and cover data
  clean_vegdat <- full_join(notes,vegdat, by = names(vegdat)[1])
  
  # NOTE > can introduce more logic checks here as needed...
  
  
  # -- CREATE TIDY LONG FORM DATASET ----
  vegdat_long <- clean_vegdat %>%
    # gather species cover only (keep pct green through litter depth in their own cols in case want to drop or break out in analysis)
    gather(code4, pct_cover, all_of(sppnames), na.rm = T) %>%
    # clean up period and number additions from cols with same code
    mutate(code4 = gsub("[.][0-9]", "", code4)) %>%
    # group by code and sum pct cover [2020 has nikolai's unqiue species names, but correspond to duplicated code4]
    grouped_df(names(.)[names(.) != "pct_cover"]) %>%
    summarise(pct_cover = sum(pct_cover)) %>%
    ungroup() %>%
    # append species descriptive info
    left_join(distinct(subset(spplist, !grepl("red california|different .Navar.|tall branchy", species) & is.na(compost_synonym), select = c(species:nativity))), by = "code4") %>%
    # order by plot, date, code
    arrange(plot, date, code4) %>%
    # rearrange cols
    select(plot, plotid, fulltrt, block:ncol(.))
  
  
  # -- APPEND LONG FORM TO MASTER LONG FORM -----
  cover_master_long <- rbind(cover_master_long, vegdat_long)
  
  
  # -- CREATE MASTER WIDE FORM FROM MASTER LONG FORM -----
  cover_master_wide <- cover_master_long[,1:grep("pct_cover",colnames(cover_master_long), fixed = T)] %>%
    spread(code4, pct_cover, fill = 0) %>%
    arrange(plot, date)
  
  # if end of loop, print done
  if(i == vegfiles[length(vegfiles)]){
    print("Master cover dataset compiled! Inspect and if all looks good write out and proceed to analysis! (w00t w00t!)")
  }
}

# -- CHECK UPDATED MAY IDs ----
# if unknowns from April each year were ID'd in May, update those (on a plot basis)
# > CTW looked back through 2019 data and only things unk in April were true unknowns (too baby to ID)
# > if anything CYEC and maybe DICA need to be updated for 2021 April data, but will compare systematically

# add rowid to coverlong for checking
cover_master_long$rowid <- as.numeric(rownames(cover_master_long))

# ctw dediced to make function to check all yrs anyway
check_overlap <- function(yrnum){
  checkdat <- mutate(cover_master_long, mon = month(date)) %>%
    subset(yr == yrnum, select = c(rowid, plot, plotid, yr, mon, date, code4, species, Family, fxnl_grp, unknown, pct_cover)) %>%
    arrange(plot, code4, date) %>%
    # check how many times code appears per plot
    group_by(plot, code4) %>%
    mutate(nobs_code = length(date)) %>%
    ungroup() %>%
    # keep only spp that appear once
    filter(nobs_code == 1) %>%
    # only keep april unk and anything from may -- only unknowns are forbs so can exclude anything from may that's a grass
    filter((mon == 4 & unknown == 1) | mon == 5) %>%
    # for the sake of screening, make N-fixer forb to potentially match unknown forb
    mutate(fxnl_grp2 = gsub("N-fixer", "Forb", fxnl_grp)) %>%
    # only keep the groups that match the unk groups in Apr
    group_by(plot) %>%
    mutate(aprgrp = str_flatten(unique(fxnl_grp[mon == 4]))) %>%
    filter(fxnl_grp2 == aprgrp) %>%
    filter(length(unique(mon)) >1) %>%
    # check for family overlap
    # replace NA with Unk (for Dead forb)
    replace_na(list(Family = "Unknown")) %>%
    # > # May fam is in Apr Fams
    mutate(flag_fam = grepl(str_flatten(unique(Family[mon == 5])), str_flatten(unique(Family[mon == 4]))) |
             # > or Apr fam is in May fams
             grepl(str_flatten(unique(Family[mon == 4])), str_flatten(unique(Family[mon == 5])))) %>%
    ungroup() %>%
    # only keep what has potential match
    filter(flag_fam) %>%
    arrange(plot, mon)
  return(checkdat)
}

check19 <- check_overlap(2019) 
# > only potential correction here is hairy aster germs in apr and hair aster sp. 1 in may (make same?)
# > check plots 1, 11, 16 for asters
check20 <- check_overlap(2020) # nothing bc only sampled once
check21 <- check_overlap(2021)
# > plot 17 can combo may DICA with Apr unk. bulb, maybe also may CESO with Apr unk. aster (check if it was smooth leaves, hairy stem entry)
# > plot 18 can combo HYRA with Apr unk. hairy asters
# > plot 24 can combo TRLA with Apr. unk bulb

# manual review..
names(cover_list)
# 2019
apr19 <- cover_list[["April2019"]]
may19 <- cover_list[["May2019"]]
# add names to make it easier to work with
names(apr19) <- paste0("p", apr19[grepl("plot", apr19[,1], ignore.case = T),])
names(may19) <- paste0("p", may19[grepl("plot", may19[,1], ignore.case = T),])
View(apr19[grep("plot|1FX|2CD|2FD", names(apr19), ignore.case = T)])
View(may19[grep("plot|1FX|2CD|2FD", names(may19), ignore.case = T)])
# what was the fate of generic asters in april?
review19 <- cover_master_long %>%
  subset(yr == 2019 & Family == "Asteraceae") %>%
  mutate(mon = month(date)) %>%
  group_by(plot) %>%
  filter(length(unique(mon))>1) %>%
  ungroup() %>%
  dplyr::select(rowid, plot, yr, date, mon, code4, species, pct_cover, Family) %>%
  arrange(plot, mon) %>%
  # kick out any codes that repeat in may and april (no discrep)
  group_by(plot, code4) %>%
  mutate(nobs_code = length(unique(mon))) %>%
  filter(nobs_code == 1) %>%
  # drop filago and chamomila (not the likely spp for aster unks)
  filter(!grepl("Filago|Chamomil", species)) %>%
  ungroup() %>%
  group_by(plot) %>%
  # drop any plots that only have 1 month now
  filter(length(unique(mon))>1)

# after review drop smooth aster because only going to address hairy asters
review19 <- filter(review19, !grepl("smooth", species, ignore.case = T)) %>%
  # keep only plots that have obs in each month
  group_by(plot) %>%
  # drop any plots that only have 1 month now
  filter(length(unique(mon))>1) %>%
  ungroup() 

replace19 <- review19 %>%
  group_by(plot) %>%
  # only keep the plots that have 1 hairy aster in may for updating hairy aster in apr
  filter(length(code4[mon == 5]) == 1) %>% # these are the plots to replace 
  # drop nobs code
  dplyr::select(-nobs_code)


# 2021
apr21 <- cover_list[["April2021"]]
may21 <- cover_list[["May2021"]]
# add names to make it easier to work with
names(apr21) <- paste0("p", apr21[grepl("plot", apr21[,1], ignore.case = T),])
names(may21) <- paste0("p", may21[grepl("plot", may21[,1], ignore.case = T),])
View(apr21[c("pPlot:", "p17", "p18", "p24")]) # p17 entry is "Aster germinants (smooth)", raw data sheet has "? 1", in native cover apr smooth turned out to be CESO to feel okay about changing, esp since same cover amt each month
View(may21[c("pPlot:", "p17", "p18", "p24")])
# also check feathery germ in Apr for NAPU in May
feather <- apply(apr21[grepl("feathery", apr21$'pPlot:'),], 1, is.na)
View(apr21[,!feather]) # in plot 21 in apr, no NAPU there in may
# what was in plot 21 in may?
View(may21[c("pPlot:", "p21")]) # I don't see anything that could be a feathery forb; checked drawing on apr datasheet (IMG_3675.HEIC). Leave be.

replace21 <- subset(check21, Family %in% c("Liliaceae", "Asteraceae") & plot %in% c(17, 18, 24)) %>%
  dplyr::select(names(replace19)) %>%
  # drop anything that only has one fam per plot
  group_by(plot, Family) %>%
  filter(length(pct_cover)>1) %>%
  ungroup()

# CTW looked at all asters over 2019-2021 in plots.. mostly likely will want to code all hair asters as HYRA and smooth as HYGL to be consistent over year
# will defer to AS to do what is best for analysis. CTW/AS called hairy aster as sp 1 most often in 2019, NS called it hypochaeris + madia, AS/CE called it hypochaeris. Some agoseris sprikled in 2019 and 2021.
replacedf <- rbind(replace19, replace21)

for(y in unique(replacedf$yr)){
  # isolate to yr of interest
  yrdat <- subset(replacedf, yr == y)
  for(i in unique(yrdat$plot)){
    tempdat <- subset(yrdat, plot == i)
    for(f in unique(tempdat$Family)){
      # pull apr info
      aprtemp <- subset(tempdat, plot == i & mon == 4 & Family == f)
      # pull may info
      maytemp <- subset(tempdat, plot == i & mon == 5 & Family == f)
      # replace code and spp info
      cover_master_long$code4[aprtemp$rowid] <- cover_master_long$code4[maytemp$rowid]
      cover_master_long[aprtemp$rowid, names(spplist)[!grepl("code4|synonym", names(spplist), ignore.case = T)]] <- cover_master_long[maytemp$rowid, names(spplist)[!grepl("code4|synonym", names(spplist), ignore.case = T)]] 
      # add QA note
      aprnote <- cover_master_long$notes[aprtemp$rowid]
      stopifnot(length(aprnote) ==1)
      tempnote <- paste0("QA note: CTW updated plot ", aprtemp$plot, " Apr unknown based on May data. Original plot ", aprtemp$plot, " Apr species was: ", aprtemp$code4, ", ", aprtemp$species)
      if(is.na(aprnote)){
        cover_master_long$notes[aprtemp$rowid] <- tempnote
      }else{
        cover_master_long$notes[aprtemp$rowid] <- paste(aprnote, tempnote, sep = "; ")
      }
    }
  }
}

#review
View(subset(cover_master_long, rowid %in% replacedf$rowid))
# check no duplicate codes per plot per year per month
mutate(cover_master_long, mon = month(date)) %>%
  group_by(plot, yr, mon, code4) %>%
  summarise(nobs = length(pct_cover)) %>%
  summary() # no duplicates


# -- CLEAN UP AND REMAKE MASTER COVER DATASETS ----
# add sample event to dataset (apr = 1, may = 2, except 2020 has only 1 visit [all dates late april])
cover_master_long <- mutate(cover_master_long, sample_event = ifelse(month(date)==4, 1, 2)) %>%
  dplyr::select(plot:yr, sample_event, date:nativity)
# check
str(cover_master_long)
# check dates by sample event
sapply(split(cover_master_long$date, paste(cover_master_long$yr, cover_master_long$sample_event)), unique) # looks good

# note: for wide, need to consolidate notes because QA note added at species level creates conflict with general plot note
cover_master_wide <- cover_master_long[,1:grep("pct_cover",colnames(cover_master_long), fixed = T)] %>%
  group_by(plot, yr, sample_event) %>%
  mutate(notes = str_flatten(unique(notes[!is.na(notes)]), collapse = "; ")) %>%
  ungroup() %>%
  spread(code4, pct_cover, fill = 0) %>%
  arrange(plot, date)
# check notes
notecheck <- subset(cover_master_wide, grepl("QA note", cover_master_wide$notes), select = c(plot, yr, sample_event, notes))
notecheck$notes # fix records 3 and 12
# manual edit via gsub
cover_master_wide$notes <- gsub("\"Other\" cover includes lichen; \"Other\" cover includes lichen;","\"Other\" cover includes lichen;", cover_master_wide$notes)
cover_master_wide$notes[cover_master_wide$plot == notecheck$plot[12] & cover_master_wide$yr == notecheck$yr[12] & cover_master_wide$sample_event == notecheck$sample_event[12]] <- "QA note: CTW updated plot 17 Apr unknowns (2) based on May data. Original plot 17 Apr species were: AST3, Smooth aster sp.; UNBU, Unknown bulb - canoe tips"

# check str
str(cover_master_wide)
# are there 36 plots per yr x sample event?
group_by(cover_master_wide, yr, sample_event) %>%
  summarise(length(plot)) # yes

# check rowid removed
names(cover_master_long)
names(cover_master_wide)


# -- PRELIM PLOT -----
# quick check that cover dat looks okay
# spp richness
data.frame(cover_master_long) %>%
  mutate(mo = month(date, label = T, abbr = T),
         date_group = paste(mo, yr, sep = " "),
         date_group = factor(date_group, levels = c("Apr 2019", "May 2019", "Apr 2020", "Apr 2021", "May 2021")),
         nut_trt = gsub("N", "XC", nut_trt)) %>%
  grouped_df(c("plot", "nut_trt", "ppt_trt", "date_group")) %>%
  summarise(S = length(unique(code4))) %>%
  ungroup() %>%
  ggplot(aes(plot, S, group = plot)) +
  geom_line(aes(group = plot)) +
  geom_point(aes(fill = date_group), pch = 21, alpha = 0.6) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_viridis_d(name = "Sample date") +
  ggtitle(paste0("USDA Compost data QA (", Sys.Date() ,"):\nprelim spp comp, richness by plot, treatments, and sample date")) +
  facet_grid(nut_trt~ppt_trt) # AS and CTW sampling in 2019 vs NS sampling in 2020.. influence in spp richness
# save to QA figs
ggsave(filename = paste0(datpath, "Cover/Cover_Figures/Prelim_QAFigures/SppComp_S_PlotxYr.pdf"),
       width = 5.25, height = 4, units = "in", scale = 1.1) 

# what about cover by family?
cover_master_long %>%
  mutate(mo = month(date, label = T, abbr = T),
         date_group = paste(mo, yr, sep = " "),
         date_group = factor(date_group, levels = c("Apr 2019", "May 2019", "Apr 2020", "Apr 2021", "May 2021"))) %>%
  replace_na(list(Family = "Unknown")) %>%
  grouped_df(c("plot","date_group", "Family")) %>%
  summarise(totcov = sum(pct_cover)) %>%
  ungroup() %>%
  ggplot(aes(plot, totcov)) +
  geom_line(aes(group = plot)) +
  geom_point(aes(fill = date_group), pch = 21, alpha = 0.6) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_viridis_d(name = "Sample date") +
  ggtitle(paste0("USDA Compost data QA (", Sys.Date() ,"): prelim spp comp, Family total cover by plot and sample date")) +
  facet_wrap(~Family, scales = "free") # interesting switch with poa vs forbs in 2020 for lower blocks..
# save to QA figs
ggsave(filename = paste0(datpath, "Cover/Cover_Figures/Prelim_QAFigures/SppComp_Familycover_PlotxYr.pdf"),
       width = 9, height = 6, units = "in", scale = 1.1) 



# -- FINISHING -----
write.csv(cover_master_long, paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_LongClean.csv"), row.names = F)
write.csv(cover_master_wide, paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_WideClean.csv"), row.names = F)
