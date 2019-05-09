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


# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed to use USDA plants api
library(request) # to access USDA plants api
library(tibble) # request reads in data to tibble
library(dplyr) # request functions rely on pipes  
options(stringsAsFactors = F)
na_vals <- c(" ", "", NA, "NA")

# set path to compost data (main level)
# dependent on user, comment out path(s) that aren't pertinent to you
datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/" #ctw's path
#datpath <- "~/Dropbox/USDA-compost/Data/" # should work for LMH and AS

# list files in entered data folder
vegfiles <- list.files(paste0(datpath, "Cover/Cover_EnteredData"), full.names = T, pattern = "_Cover_", ignore.case = T)


# -- COMPILE SPP LIST -----
# loop through cover data to read in, transpose, and append to master cover dataset
# initiate master data frame
spplist_master <- data.frame()
# loop iterates through each cover dataset and adds new spp to master spp set
for(i in vegfiles){
  # read in dataset
  vegdat <- read.csv(i, na.strings = na_vals, header = F, blank.lines.skip = T, strip.white = T)
  print(paste("Adding new species from", vegdat$V2[2], "composition survey to master species list"))
  
  # id where litter depth and cover data start
  litpos <- grep("depth", vegdat[,1], ignore.case = T)
  # remove any species rows that don't have any entries for abundance value
  allNAs <- apply(vegdat[litpos+1:nrow(vegdat), 2:ncol(vegdat)],1,function(x) all(is.na(x)))
  # pull species list
  spplist <- vegdat[litpos+1:nrow(vegdat),1]
  spplist <- spplist[allNAs==FALSE]
  # compile spp list data frame for i cover dataset, remove NAs from list
  spplist <- data.frame(species = sort(spplist[!is.na(spplist)]),
                        unknown = NA,
                        genus = NA,
                        epithet = NA,
                        code4 = NA,
                        code6 = NA)
  # identify unknowns in species list
  spplist$unknown[grepl("p[.]$|[(]|[0-9]", spplist$species)] <- 1
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


# -- APPEND USDA PLANTS DATA -----
# specify vars desired from usda plants database (there are 134)
usda_plantvars <- c("Symbol","Accepted_Symbol_x","Scientific_Name_x","Common_Name","State_and_Province",
                    "Category","Family","Family_Common_Name","Duration","Growth_Habit","Native_Status")

spplist_master <- cbind(spplist_master, data.frame(matrix(nrow = nrow(spplist_master), ncol=length(usda_plantvars))))
colnames(spplist_master)[which(colnames(spplist_master) == "X1"):ncol(spplist_master)] <- usda_plantvars
# run usda plants api query to scrape species info
# NOTE!!: this will throw an error ["Client error: (400) Bad Request"] if a species is spelled incorrectly in the cover data (a good QA check)
# loop will issue warnings about cbind command providing more variables to replace than there are in 
for(p in spplist_master$species[spplist_master$unknown == 0]){
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
  if(!grepl(" ssp[.]|-", p)){
    templist <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Genus = eval(temp_genus), Species = eval(temp_epithet))
  }
  # isolate desired cols
  temp_df <- templist$data[1,colnames(templist$data) %in% usda_plantvars]
  # rematch to updated name if accepted symbol doesn't match symbol
  if(temp_df$Symbol != temp_df$Accepted_Symbol_x){
    templist2 <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Symbol = eval(temp_df$Accepted_Symbol_x))
    update_df <- templist2$data[1,colnames(templist2$data) %in% usda_plantvars]
    temp_df[,which(colnames(temp_df)=="Common_Name"):ncol(temp_df)] <- update_df[,which(colnames(update_df)=="Common_Name"):ncol(update_df)]
  }
  # cleanup empty cells
  temp_df[temp_df==""] <- NA
  # append to master data frame
  #usdaplants_df <- rbind(usdaplants_df, temp_df)
  # add to spplist_master
  spplist_master[spplist_master$species == p,usda_plantvars] <- as.data.frame(temp_df)
}


# -- FINISHING -----
# manual edits: fill in info for unknowns
# unknown Avena sp.
spplist_master[grepl("Avena sp",spplist_master$species), c("genus", "code4", "code6")] <- c("Avena", "AVSPP", "AVESPP")
spplist_master[grepl("Avena sp",spplist_master$species), 
               which(colnames(spplist_master)=="Category"):ncol(spplist_master)] <- spplist_master[spplist_master$code4 == "AVBA" & !is.na(spplist_master$code4), 
                                                                                                   which(colnames(spplist_master)=="Category"):ncol(spplist_master)]
# unknown Trifolium 
for(i in which(grepl("Trif.*[0-9]", spplist_master$species))){
  num <- gsub("[^0-9]","", spplist_master$species[i])
  spplist_master[i, c("genus", "code4", "code6")] <- c("Trifolium", paste0("TRI",num), paste0("TRISP",num))
  # match descriptive info with TRHI (most generic)
  spplist_master[i, which(colnames(spplist_master)=="Category"):ncol(spplist_master)] <- spplist_master[spplist_master$code4 == "TRHI" & !is.na(spplist_master$code4), 
                                                                                                        which(colnames(spplist_master)=="Category"):ncol(spplist_master)]
}
# name remaining forbs as unk forb #
num2 <- 5 #last number assigned for unknowns is 4
for(i in which(is.na(spplist_master$genus))){
  num <- gsub("[^0-9]","", spplist_master$species[i])
  # if no number already assigned, use num2 count
  if(nchar(num)==0){
    num <- num2
    # increase num2 by 1 for next unknown species
    num2 <- num2+1
  }
  spplist_master[i, c("code4", "code6", "Growth_Habit")] <- c(paste0("UNKF",num), paste0("UNKFRB",num), "Forb/herb")
  # assign asteraceae info if an unknown aster
  if(grepl(" aster ", spplist_master$species[i])){
    spplist_master[i, c("Category", "Family", "Family_Common_Name")] <- c("Dicot", "Asteraceae", "Aster family")
  }
}

# finish by adding fxnl_grp and simplified nativity col
spplist_master$fxnl_grp[grepl("Gram", spplist_master$Growth_Habit)] <- "Grass"
spplist_master$fxnl_grp[grepl("Forb", spplist_master$Growth_Habit, ignore.case = T)] <- "Forb"
spplist_master$fxnl_grp[spplist_master$Family == "Fabaceae"] <- "N-fixer"

spplist_master$nativity[grepl("L48 .I.",spplist_master$Native_Status)] <- "Exotic"
spplist_master$nativity[grepl("L48 .N.",spplist_master$Native_Status)] <- "Native"


# -- WRITE OUT -----
write.csv(spplist_master, paste0(datpath, "Compost_SppList.csv"), row.names = F)
