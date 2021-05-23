# build usda-compost species list
# author(s): ctw (caitlin.t.white@colorado.edu)
# date created: 2019-04-28 (will be modified over project lifetime)

# script purpose:
# add new species not present in the master compost cover species list (or start one [e.g. April 2019])
# add in-house plant codes, for ID'd species and unknowns
# append USDA Plants database info to species list through R package 'request' by Scott Chamberlain (https://cran.r-project.org/web/packages/request/request.pdf)
# > note: USDA Plants Databse does not have its own API, so Scott created package for that but hasn't maintained it since 2016 (still works for simple scraping)
# add simplified functional group and nativity columns based on USDA Plants data
# write out to main Data folder on compost dropbox so can be linked to cover and plant traits analyses (or anything that analyzes species level data)
# > note: can add in datasets to read apart from cover data if for some reason cover datasets do not capture all compost plant species of interest

# requirements to run script:
#devtools::install_github("sckott/request")

# notes 2019 apr 28: 
# 1. can be reworked a little when have more than 1 cover dataset to append if don't want to run through all cover datasets (but doesn't matter either way)
# 2. script relies on cover dataset having "_Cover_" in the name (okay if all lower case)
# can change script so reads in all files living in the Cover_EnteredData folder, but that means only entered cover data should live in that folder
# the way it's written currently is more flexible should anything other than entered cover data live in there

# notes 2020 jun 30:
# CTW still waiting to hear back from Nikolai on some of the species noted in composition.. e.g. "skinny" and "tall branchy" -- not sure if those are forbs or grams
# we may want to recode NS's "tarweed" if don't think it's actually Madia (native aster)

# notes from 1 july 2020:
# The question on Navarretia…the name in question is “Navarretia Junior” the other version of nav. Much smaller and daintier looking. Ashley said later that they were both the same but in different life stages. The 70 % cover was for the “Junior” and “regular Nav.” was 10.
# The question on invasive vs. native red brome…I figure it was invasive seeing the large seed heads and brown stems.
# I sent a photo on the “Tall Branchy”. (Silene gallica -- sent to AS on April 30)
# As for the “Skinny” I am still processing that request. (sent April 30 too? it's a forb.. don't know what, too undeveloped.. but could be unknown forb from last year)

# notes from 15 may 2021:
# Scott Chamberlain took down his RESTful API for USDA Plants in Sep 2020 per USDA saying they have an API to use..
# .. they don't (CTW is super mad)
# new workflow is to read in existing spplist compiled through summer 2020 and then use .rdata USDA PLANTS DB I found (https://uofi.app.box.com/v/usdaplants)
# import data from .rdata source for any spp that aren't already on Compost spp list
# appears info in .rdata come from 2011-2016 USDA PLANTS info, which was as current as S. Chamberlain API info (his was from 2011 I think)
# CTW is super annoyed with USDA (so mad I am documenting it. grrrrrr!!)



# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed to use USDA plants api
library(request) # to access USDA plants api
library(tidyverse) # for tibble, dplyr, stringr
options(stringsAsFactors = F)
na_vals <- c(" ", "", NA, "NA", ".")

# set path to compost data (main level)
# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# list files in entered data folder
## look at both main cover dat and native recruitment cover 
vegfiles <- list.files(paste0(datpath, "Cover/Cover_EnteredData"), full.names = T, pattern = "_Cover_", ignore.case = T)
# append native recruitment
vegfiles <- c(vegfiles, list.files(paste0(datpath, "Native/Native_EnteredData"), full.names = T))

# take out crosswalk files (to use in native prep only)
vegfiles <- vegfiles[!grepl("crosswalk", vegfiles)]

# read in USDA PLANTS .rdata from Univ. Ill.
load(paste0(datpath, "USDA_PLANTS/usdaplants.RData"))
# read in current Compost spplist
spplist_current <- read.csv(paste0(datpath, "Compost_SppList.csv"), na.strings = na_vals) 



# -- TIDY USDA PLANTS DATA -----
names(usdaplants)
# specify cols to keep (to match API process below)
rdat_vars <- c("Symbol","AcceptedSymbol","ScientificName","CommonName","State",
               "Category","Family","FamilyCommonName","Duration","GrowthHabit","NativeStatus")

tidyusda  <- usdaplants[rdat_vars]
# NA empty cells
tidyusda <- data.frame(sapply(tidyusda, function(x) ifelse(x %in% c("", " "), NA, x)))
# clean up state info
head(tidyusda$State)
head(spplist_current$State_and_Province)
tidyusda$State <- gsub("<strong>|</strong>", "", tidyusda$State)
# check results
sort(unique(tidyusda$State))
summary(unique(tidyusda$State) %in% unique(spplist_current$State_and_Province))
sort(unique(spplist_current$State_and_Province))
spplist_current$State_and_Province[spplist_current$species == "Avena barbata"]
tidyusda$State[tidyusda$ScientificName == "Avena barbata"] # will look at this more after run through species



# -- COMPILE SPP LIST -----
# loop through cover data to read in, transpose, and append to master cover dataset
# initiate master data frame
spplist_master <- data.frame()
# loop iterates through each cover dataset and adds new spp to master spp set
for(i in vegfiles){
  # read in dataset
  vegdat <- read.csv(i, na.strings = na_vals, header = F, blank.lines.skip = T, strip.white = T)
  print(paste("Adding new species from", str_extract(i, "(?<=_)[:alnum:]+(?=.csv)"), "composition survey to master species list"))
  
  # id where litter depth and cover data start
  litpos <- grep("depth", vegdat[,1], ignore.case = T)
  # id col where plots/abundance data start (2020 data start in different col than 2019)
  # > this assumes litter depth has a numeric value (which it should..)
  # > may 2021 dats don't always if just an update on april data. force to be col that corresponds to first plot if empty
  if(sum(grep("[[:digit:]]", vegdat[litpos,])) > 0){
    covpos <- min(grep("[[:digit:]]", vegdat[litpos,]))
  }else{
    # id plot pos
    plotpos <- grep("plot", vegdat[,1], ignore.case = T)
    covpos <- min(grep("[[:digit:]]", vegdat[plotpos,]))
  }
  # remove any species rows that don't have any entries for abundance value
  allNAs <- apply(vegdat[litpos+1:nrow(vegdat), covpos:ncol(vegdat)],1,function(x) all(is.na(x)))
  # pull species list
  spplist <- vegdat[litpos+1:nrow(vegdat),1]
  spplist <- spplist[!allNAs]
  # compile spp list data frame for i cover dataset, remove NAs from list
  spplist <- data.frame(species = sort(spplist[!is.na(spplist)]),
                        unknown = NA,
                        genus = NA,
                        epithet = NA,
                        code4 = NA,
                        code6 = NA)
  # identify unknowns in species list
  spplist$unknown[grepl("p[.]$|[(]|[0-9]|group", spplist$species)] <- 1
  spplist$unknown[is.na(spplist$unknown)] <- 0
  
  # fill in cols for known species
  spplist$genus[spplist$unknown == 0] <- trimws(gsub(" .*", "", spplist$species[spplist$unknown == 0]))
  spplist$epithet[spplist$unknown == 0] <- trimws(gsub("^[A-Z][a-z]+ ", "", spplist$species[spplist$unknown == 0]))
  spplist$code4[spplist$unknown == 0] <- with(spplist[spplist$unknown == 0,], casefold(paste0(substr(genus,1,2), substr(epithet, 1,2)), upper = T))
  spplist$code6[spplist$unknown == 0] <- with(spplist[spplist$unknown == 0,], casefold(paste0(substr(genus,1,3), substr(epithet, 1,3)), upper = T))
  
  # if any species have same 4-letter code or same 6-letter code, enumerate
  code4_dups <- spplist$code4[duplicated(spplist$code4) == T & !is.na(spplist$code4)]
  for(c4 in code4_dups){
    num <- 1
    for(o in which(spplist$code4 == c4)){
      spplist$code4[o] <-paste0(spplist$code4[o],num)
      num <- num+1
    }
  }
  code6_dups <- spplist$code6[duplicated(spplist$code6) == T & !is.na(spplist$code6)]
  for(c6 in code6_dups){
    num <- 1
    for(o in which(spplist$code6 == c6)){
      spplist$code6[o] <-paste0(spplist$code6[o],num)
      num <- num+1
    }
  }
  # print new species added
  if(ncol(spplist_master) == 0){
    print(paste(length(unique(spplist$species)), "new species added:"))
    print(spplist$species)
  }else{
    newspp <- spplist$species[!spplist$species %in% spplist_master$species]
    print(paste(length(newspp), "new species added:"))
    print(newspp)
  }
  # append species to master species list data frame
  spplist_master <- rbind(spplist_master, spplist)
  # remove any duplicated species
  spplist_master <- spplist_master[!duplicated(spplist_master$species),]
}

# > only need to run this part once (done)
# add generic dichelostemma and leontodon sp. (in may notes to replace apr, but not recorded/entered in data rows)
# addgeneric <- data.frame(matrix(nrow = 2, ncol = ncol(spplist_master), dimnames = list(NULL, names(spplist_master))))
# addgeneric$species <- c("Dichelostemma sp.", "Leontodon sp.") 
# spplist_master <- rbind(spplist_master, addgeneric)



# -- APPEND USDA PLANTS DATA -----
# specify vars desired from usda plants database (there are 134)
usda_plantvars <- c("Symbol","Accepted_Symbol_x","Scientific_Name_x","Common_Name","State_and_Province",
                    "Category","Family","Family_Common_Name","Duration","Growth_Habit","Native_Status")

spplist_master <- cbind(spplist_master, data.frame(matrix(nrow = nrow(spplist_master), ncol=length(usda_plantvars))))
colnames(spplist_master)[which(colnames(spplist_master) == "X1"):ncol(spplist_master)] <- usda_plantvars
# flag generic genera to pair on symbol instead of genus and epithet
genera <- spplist_master$species[spplist_master$unknown == 1 & grepl("sp[.]|spp[.]", spplist_master$species)]
# remove asters and unknowns geric grams and forbs
genera <- genera[!grepl("aster|forb|grass", genera, ignore.case = T)]

# ============= SKIP THESE SECTIONS IN 2021 =========================

# run usda plants api query to scrape species info
# NOTE!!: this will throw an error ["Client error: (400) Bad Request"] if a species is spelled incorrectly in the cover data (a good QA check)
# loop will issue warnings about cbind command providing more variables to replace than there are in 
for(p in c(spplist_master$species[spplist_master$unknown == 0], genera)){
  print(paste("Pulling USDA Plants data for",p))
  temp_genus <- spplist_master$genus[spplist_master$species == p]
  temp_epithet <- spplist_master$epithet[spplist_master$species == p] 
  
  # grab usda plants data
  if(grepl("ssp.", p)){
    temp_susbp <- gsub("^[A-Z].+ ssp. ", "", p)
    temp_epithet <- gsub(" .*", "", temp_epithet)
    templist <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Genus = eval(temp_genus), Species = eval(temp_epithet), Subspecies = eval(temp_susbp))
  }
  # special case for medusahead (any hyphenated species epithet is a special case, search function bonks with hyphen)
  if(grepl("Taen", p)){
    templist <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query(Genus = Taeniatherum, Species = `caput-medusae`)
  }
  # special case for generic genus
  if(p %in% genera){
    temp_symbol <- casefold(substr(p,1,5),upper = T)
    templist <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Symbol = temp_symbol)
  }
  # for all others
  if(!grepl(" ssp[.]|-", p) & !p %in% genera){
    templist <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Genus = eval(temp_genus), Species = eval(temp_epithet))
  }
  
  
  # isolate desired cols
  temp_df <- data.frame(templist$data)
  temp_df <- temp_df[,names(temp_df) %in% usda_plantvars]
  # rematch to updated name if accepted symbol doesn't match symbol
  if(temp_df$Symbol != temp_df$Accepted_Symbol_x){
    templist2 <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Symbol = eval(temp_df$Accepted_Symbol_x))
    update_df <- data.frame(templist2$data)
    update_df <- update_df[names(update_df) %in% usda_plantvars]
    temp_df[,which(colnames(temp_df)=="Common_Name"):ncol(temp_df)] <- update_df[,which(colnames(update_df)=="Common_Name"):ncol(update_df)]
  }
  # cleanup empty cells
  temp_df[temp_df==""] <- NA
  # append to master data frame
  #usdaplants_df <- rbind(usdaplants_df, temp_df)
  # add to spplist_master
  spplist_master[spplist_master$species == p,usda_plantvars] <- as.data.frame(temp_df)
}


str(spplist_master)
# save work
copydf <- spplist_master
# clean up environment
rm(temp_df, templist, templist2, update_df, temp_epithet, temp_genus, temp_susbp, temp_symbol, i, p)





# -- MANUAL INFILLING -----
# manual edits: fill in info for unknowns

# unknown Annual grass spp. (AVBA, CYEC, or GAPH)
spplist_master[grepl("Annual grass",spplist_master$species), c("code4", "code6")] <- c("UNGR", "UNGSPP")
# copy descriptive info from TACA since most general and applies
spplist_master[grepl("Annual grass",spplist_master$species), 
               which(colnames(spplist_master)=="Category"):ncol(spplist_master)] <- spplist_master[spplist_master$code4 == "TACA" & !is.na(spplist_master$code4), 
                                                                                                   which(colnames(spplist_master)=="Category"):ncol(spplist_master)]

# unknown Trifolium 
for(i in which(grepl("Trif.*[0-9]|Trifolium spp.", spplist_master$species))){
  num <- gsub("[^0-9]","", spplist_master$species[i])
  spplist_master[i, c("genus", "code4", "code6")] <- c("Trifolium", paste0("TRI",num), paste0("TRISP",num))
  # match descriptive info with TRHI (most generic)
  spplist_master[i, which(colnames(spplist_master)=="Category"):ncol(spplist_master)] <- spplist_master[spplist_master$code4 == "TRHI" & !is.na(spplist_master$code4), 
                                                                                                        which(colnames(spplist_master)=="Category"):ncol(spplist_master)]
}
# fix code 4 and 6 if only 1 unknown trifolium entry for multiple spp
spplist_master[grepl("Tri.* spp.",spplist_master$species), c("code4", "code6")] <- c("TRSP", "TRSPP")

# assign code 4 and 5 + genus for generic genera
for(i in genera[grepl("[[:digit:]]", genera)]){
  if(grepl("Bromus", i)){
    next
  }
  tempos <- which(spplist_master$species == i)
  spplist_master$genus[tempos] <- str_extract(i, "[A-Z][a-z]+")
  spplist_master$code4[tempos] <- paste0(casefold(substr(i,1,3), upper = T), str_extract(i,"[[:digit:]]+"))
  spplist_master$code6[tempos] <- paste0(casefold(substr(i,1,5), upper = T), str_extract(i,"[[:digit:]]+"))
  
}
# for generic genera without number, just use alpha chars
for(i in genera[!grepl("[[:digit:]]", genera)]){
  tempos <- which(spplist_master$species == i)
  spplist_master$genus[tempos] <- str_extract(i, "[A-Z][a-z]+")
  spplist_master$code4[tempos] <- paste0(casefold(substr(i,1,2), upper = T), "SP")
  spplist_master$code6[tempos] <- paste0(casefold(substr(i,1,3), upper = T), "SPP")
  
}
# assign codes for unk asters (3 different species, and then a general aster group for immature asters)
asters <- spplist_master$species[is.na(spplist_master$code4) & grepl("aster", spplist_master$species)]
# set where to start with numbering asters (will work for future additions)
plantnum <- max(as.numeric(str_extract(asters,"[[:digit:]]")), na.rm = T) + 1
for(i in asters){
  tempos <- which(spplist_master$species == i)
  # infill family info (can't assign nativity or duration since agoseris can be peren and is native, madia also native, but non-natives also mixed in here)
  spplist_master[tempos, c("Category","Family","Family_Common_Name", "Growth_Habit")] <- c("Dicot", "Asteraceae", "Aster family", "Forb/herb")
  
  if(grepl("[[:digit:]]", i)){
    spplist_master$code4[tempos] <- paste0("AST", str_extract(i,"[[:digit:]]"))
    spplist_master$code6[tempos] <- paste0("ASTER", str_extract(i,"[[:digit:]]"))
    next
  }
  if(grepl("spp[.]", i)){
    spplist_master$code4[tempos] <- "ASSP"
    spplist_master$code6[tempos] <- "ASTSPP"
    next
  }
  if(grepl("smooth", i, ignore.case = T)){
    spplist_master$code4[tempos] <- paste0("AST", plantnum)
    spplist_master$code6[tempos] <- paste0("ASTER", plantnum)
    plantnum <- plantnum + 1
  }
}


# infill agoseris (based on usda plants info agoseris spp on sfrec list)
spplist_master[grep("AGOSE", spplist_master$Symbol), c("Growth_Habit", "Native_Status")] <- c("Forb/herb", "L48 (N)")

# assign codes for grouped genera, and clean up aster fam
spplist_master$code4[grepl("group", spplist_master$species)] <- paste0(substr(spplist_master$species[grepl("group", spplist_master$species)],1,1),
                                                                       str_extract(spplist_master$species[grepl("group", spplist_master$species)], "(?<=-)[A-Z]"),
                                                                       "GR")
spplist_master$code6[grepl("group", spplist_master$species)] <- casefold(paste0(substr(spplist_master$species[grepl("group", spplist_master$species)],1,2),
                                                                                str_extract(spplist_master$species[grepl("group", spplist_master$species)], "(?<=-)[A-Z][a-z]"),
                                                                                "GR"), upper = T)
spplist_master$Category[grepl("group", spplist_master$species)] <- "Dicot" # all dicots
spplist_master$Duration[grepl("group|Hypochaeris sp|Madia sp|Geranium sp", spplist_master$species)] <- "Annual" # all of these are annual
spplist_master$Growth_Habit[grepl("group|Hypochaeris sp|Madia sp|Geranium sp", spplist_master$species)] <- "Forb/herb"
spplist_master$Native_Status[grepl("group|Hypochaeris sp|Geranium sp", spplist_master$species)] <- "L48 (I)" # they are at least this
spplist_master$Native_Status[grepl("Madia sp", spplist_master$species)] <- "L48 (N)"
# taraxacum and hypochaeris are both asters
spplist_master[grepl("Taraxacum-", spplist_master$species), c("Family", "Family_Common_Name")] <- c("Asteraceae", "Aster family")
# if duration not assigned to asters by now, is unknown
spplist_master$Duration[grepl("Aster", spplist_master$Family_Common_Name) & is.na(spplist_master$Duration)] <- "Unknown"

# name remaining forbs as unk forb #
plantnum <- 1
for(i in which(is.na(spplist_master$code4) & grepl("forb", spplist_master$species, ignore.case = T))){
  # if multiple forb (forb germinates), create generic group
  if(grepl("spp", spplist_master$species[i])){
    spplist_master[i, c("code4", "code6")] <- c("UNFR", "UNFSPP")
  }else{
    spplist_master[i, c("code4", "code6")] <- c(paste0("UNF",plantnum), paste0("UNFRB",plantnum))
    plantnum <- plantnum + 1
  }
  # assign other general annual non-native forb info
  spplist_master[i, c("Category", "Growth_Habit", "Duration", "Native_Status")] <- c("Dicot", "Forb/herb", "Annual", "L48 (I)") # is at least invasive in lower 48 (probably other parts of US too) 
  
}

# finish by adding fxnl_grp and simplified nativity col
spplist_master$fxnl_grp[grepl("Gram", spplist_master$Growth_Habit)] <- "Grass"
spplist_master$fxnl_grp[grepl("Forb", spplist_master$Growth_Habit, ignore.case = T)] <- "Forb"
spplist_master$fxnl_grp[grepl("Fabac", spplist_master$Family)] <- "N-fixer"
spplist_master$fxnl_grp[spplist_master$Growth_Habit == "Unknown"] <- "Unknown"

spplist_master$nativity[grepl("L48 .I.",spplist_master$Native_Status)] <- "Exotic"
spplist_master$nativity[grepl("L48 .N.",spplist_master$Native_Status)] <- "Native"
spplist_master$nativity[is.na(spplist_master$Native_Status)] <- "Unknown"


# ============= RESUME SCRIPT FOR 2021 =========================


# -- ADD IN NEW 2021 SPP ----
# assign full master to a different df so code in finishing runs as intended
spplist_master_base <- spplist_master

# separate new list for 2021
# > if nothing exists (because list full compiled, skip this part)
spplist_2021 <- subset(spplist_master, !species %in% spplist_current$species)
# designate unknowns that didn't get marked in code above
spplist_2021$species[!spplist_2021$species %in% tidyusda$ScientificName] #FEMI is VUMI
spplist_2021$unknown[!(spplist_2021$species %in% tidyusda$ScientificName) & !grepl("microstachys", spplist_2021$species)] <- 1
# NA cols for unknowns
spplist_2021[spplist_2021$unknown == 1, c("genus", "epithet", "code4", "code6")] <- NA

# pull USDA data for known spp
spplist_2021_known <- left_join(subset(spplist_2021, unknown == 0, select = c(species:code6)), tidyusda, by = c("species" = "ScientificName"), keep = T)
# add in data for VUMI
spplist_2021_known[grepl("microstach", spplist_2021_known$species), grep("^Symbol", names(spplist_2021_known)):ncol(spplist_2021_known)] <- subset(tidyusda, Symbol == "VUMI")


# clean up unknowns
spplist_2021 <- subset(spplist_2021, unknown == 1)
spplist_2021$compost_synonym <- NA

# > check what might overlap already with current species
# manual corrections
## ASTERS
spplist_2021[grepl("hairy", spplist_2021$species, ignore.case = T), grep("code4", names(spplist_2021)):ncol(spplist_2021)] <- subset(spplist_current, code4 == "ASSP", select = c(code4:Native_Status, species))
spplist_2021[grepl("[(]smooth[)]$", spplist_2021$species), grep("code4", names(spplist_2021)):ncol(spplist_2021)] <- subset(spplist_current, code4 == "AST3", select = c(code4:Native_Status, species))
spplist_2021[grepl("catpaw", spplist_2021$species), grep("code4", names(spplist_2021)):ncol(spplist_2021)] <- subset(spplist_current, code4 == "HYSP", select = c(code4:Native_Status, species))
spplist_2021[grepl("Leontodon sp.", spplist_2021$species), grep("^Symbol$", names(spplist_2021)):(ncol(spplist_2021)-1)] <- subset(tidyusda, ScientificName == "Leontodon")
spplist_2021[grepl("Leontodon sp.", spplist_2021$species), c("genus", "code4", "code6")] <- c("Leontodon", "LEON", "LEONSP")
spplist_2021[grepl("Leontodon sp.", spplist_2021$species), c("Duration", "Growth_Habit", "Native_Status")] <- subset(tidyusda, Symbol == "LETA", select = c(Duration:NativeStatus))

## TRIFOL
spplist_2021[grepl("smooth heart", spplist_2021$species), grep("code4", names(spplist_2021)):ncol(spplist_2021)] <- subset(spplist_current, code4 == "TRSP", select = c(code4:Native_Status, species))

## assign unk bulb (unsure if genus will be tritel. or dichelo.)
spplist_2021[grepl("bulb", spplist_2021$species, ignore.case = T), c("code4", "code6")] <- rep(c("UNBU", "UNBULB"), each = 2)
spplist_2021[grepl("bulb", spplist_2021$species), grep("^State", names(spplist_2021)):(ncol(spplist_2021)-1)] <- subset(tidyusda, Symbol == "DICA14", select = c(State:NativeStatus))
spplist_2021[grepl("Dichelostemma sp.", spplist_2021$species), grep("^Symbol$", names(spplist_2021)):(ncol(spplist_2021)-1)] <- subset(tidyusda, ScientificName == "Dichelostemma")
spplist_2021[grepl("Dichelostemma sp.", spplist_2021$species), c("genus", "code4", "code6")] <- c("Dichelostemma", "DICH", "DICHSP")
spplist_2021[grepl("Dichelostemma sp.", spplist_2021$species), c("Duration", "Growth_Habit", "Native_Status")] <- subset(tidyusda, Symbol == "DICA14", select = c(Duration:NativeStatus))

# assign remaining unk forbs (currently only 1 in spplist)
num <- sum(grepl("UNFRB[0-9]+", spplist_current$code6)) # more dynamic way of counting
for(s in spplist_2021$species[is.na(spplist_2021$code4)]){
  spplist_2021$code4[spplist_2021$species == s] <- paste0("UNF", (num+1))
  spplist_2021$code6[spplist_2021$species == s] <- paste0("UNFRB", (num+1))
  num <- num+1
}
# assign other info
spplist_2021[grepl("UNFRB", spplist_2021$code6), grep("Category", names(spplist_2021)):(ncol(spplist_2021)-1)] <- subset(spplist_current, code4 == "UNF1", select = c(Category:Native_Status))

# stack knowns to unknowns and assign fxnl group then nativity to append to master
names(spplist_2021)
names(spplist_2021_known)
spplist_2021_known$compost_synonym <- NA
names(spplist_2021_known) <- names(spplist_2021) 
spplist_2021 <- rbind(spplist_2021, spplist_2021_known)

# finish by adding fxnl_grp and simplified nativity col
spplist_2021$fxnl_grp[grepl("Gram", spplist_2021$Growth_Habit)] <- "Grass"
spplist_2021$fxnl_grp[grepl("Forb", spplist_2021$Growth_Habit, ignore.case = T)] <- "Forb"
spplist_2021$fxnl_grp[grepl("Fabac", spplist_2021$Family)] <- "N-fixer"
spplist_2021$fxnl_grp[spplist_2021$Growth_Habit == "Unknown"] <- "Unknown"

spplist_2021$nativity[grepl("L48 .I.",spplist_2021$Native_Status)] <- "Exotic"
spplist_2021$nativity[grepl("L48 .N.",spplist_2021$Native_Status)] <- "Native"
spplist_2021$nativity[is.na(spplist_2021$Native_Status)] <- "Unknown"

# create all years master
spplist_master <- rbind(spplist_current, spplist_2021)


# -- FINISHING -----
# order alphabetically by species
spplist_master <- spplist_master[order(spplist_master$species),]
# remove any spp no longer in ongoing veg datasets
spplist_master <- subset(spplist_master, species %in% spplist_master_base$species)

# qa checks
# > check for typos
sort(unique(spplist_master$Family))
sort(unique(spplist_master$Growth_Habit))
sort(unique(spplist_master$Duration))
# check who is unknown for nativity
sort(spplist_master$species[spplist_master$nativity == "Unknown"])
# check gram and forb spp look good
sort(spplist_master$species[grepl("forb", spplist_master$Growth_Habit, ignore.case = T)])
sort(spplist_master$species[grepl("gram", spplist_master$Growth_Habit, ignore.case = T)])
# check NAs
sapply(spplist_master, function(x) summary(is.na(x))) # looks okay -- NAs present for unknowns



# -- WRITE OUT -----
write_csv(spplist_master, paste0(datpath, "Compost_SppList.csv"))
