# read in, compile, and clean all anpp dats for analysis
# author(s): ctw
# date create: may 2021

# notes: as of 5/14/2021, will still need to add in clip #2 for 2021. 
## CTW discovered 2 plots weren't clipped in june 2020. if those samples turn up, can add in
## also may 2019 scanned datasheets are not posted to raw files. if/when those add, can update "file" in dataset as well.


# -- SETUP --
library(tidyverse)
library(lubridate)
# modify default settings
options(stringsAsFactors = F)
na_vals <- c("" , " ", "NA", NA)

# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# list files in entered data folder
anpp_files <- list.files(paste0(datpath, "ANPP/ANPP_EnteredData"), full.names = T, pattern = "20([0-9]{2})", ignore.case = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)


# -- COMPILE ANPP ----
#initiate df for all data
anpp_master <- data.frame()

for(i in anpp_files){
  temp_anpp <- read.csv(i, na.strings = na_vals)
  if(grepl(2019, i)){
    # clean up
    temp_anpp <- unite(temp_anpp, plotid, plot, subplot, sep = "") %>%
      mutate_at(c("date_clip", "date_sort", "date_weigh"), as.Date, format = "%m/%d/%y") %>%
      # add yr, month, and yday
      mutate(yr = year(date_clip), mon = month(date_clip), doy = yday(date_clip)) %>%
      # remove empty rows
      subset(!(is.na(dry_wgt_g) & is.na(notes))) %>%
      # write out fxnl grp
      mutate(fxnl_grp = ifelse(grepl("G", fxnl_grp, ignore.case = T), "Grass", "Forb")) %>%
      # join trt key
      left_join(trtkey) %>%
      # reorganize cols
      select(file:date_clip, yr:doy, names(trtkey), fxnl_grp:ncol(.)) %>%
      arrange(page, line)
  }else{
    # all other dats are entered in the same way
    temp_anpp <- temp_anpp %>%
      # remove empty rows
      subset(!(is.na(dry_wgt_g) & is.na(notes))) %>%
      mutate_at(c("date_clip", "date_sort", "date_weigh"), as.Date, format = "%m/%d/%y") %>%
      # add yr, month, and yday
      mutate(yr = year(date_clip), mon = month(date_clip), doy = yday(date_clip)) %>%
      # write out fxnl grp
      mutate(fxnl_grp = ifelse(grepl("G", fxnl_grp, ignore.case = T), "Grass", "Forb")) %>%
      # join trt key
      left_join(trtkey) %>%
      # reorganize cols
      select(file:date_clip, yr:doy, names(trtkey), fxnl_grp:ncol(.)) %>%
      arrange(page, line)
  }
  #append to master
  anpp_master <- rbind(anpp_master, temp_anpp)
}

# tidy dataset by date and add yrly clip order
anpp_master <- group_by(anpp_master, yr) %>%
  mutate(clip_order = as.numeric(as.factor(date_clip)),
         # 1 = target forb, 2 = target grams
         clip_event= as.numeric(as.factor(mon))) %>%
  ungroup() %>%
  select(file:date_clip, clip_order, clip_event, yr:ncol(.)) %>%
  arrange(yr, date_clip, page, line)

# check each plot has value for each fxnl grp per yr per clip event (sampling trip)
length(unique(anpp_master$plot)) # should be 36 nobs per group
group_by(anpp_master, yr, clip_event, fxnl_grp) %>%
  summarise(nobs = length(plot)) # two plots missed in 2nd clip of 2020



# -- QUALITY SCREEN ----
summary(anpp_master)
summary(is.na(anpp_master)) #infill missing file until it gets posted, clarify "not recorded" for anpp processing info
anpp_master <- replace_na(anpp_master, list(file = "raw file missing", sort_init = "not recorded", weigh_init = "not recorded")) %>%
  # add in not recorded for sort and weigh dates -- need to convert to character class to work
  mutate_at(c("date_sort", "date_weigh"), function(x) ifelse(is.na(x), "not recorded", as.character(x)))
# check unique values  
sapply(anpp_master, function(x) sort(unique(x))) #looks ok

# quick plot to screen for wonky vals
ggplot(anpp_master, aes(factor(clip_event), dry_wgt_g, fill = factor(yr), group = paste(clip_event, yr))) +
  geom_boxplot(alpha = 0.5, col = "grey50") +
  geom_point(aes(shape = fxnl_grp, group = paste(clip_event, yr)), size = 2, position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.3, )) +
  facet_grid(nut_trt~ppt_trt) +
  scale_shape_manual(values= c("F", "G")) +
  theme_test() # CD missing a data point in 2020 clip 2

with(subset(anpp_master, yr == 2020 & clip_event == 2), sapply(split(plot, fulltrt), unique)) # no plot 4 in NW, also missing 11 plot for CD trt
# check other years to be sure
group_by(anpp_master, yr, clip_event, fulltrt) %>%
  summarise(nobs = length(unique(plot))) %>%
  filter(nobs < 4) # just 2020. CTW double-checked datasheets and plots 4 and 11 are not there for june clip

# save plot to Dropbox with annotations in title
ggplot(anpp_master, aes(factor(clip_event), dry_wgt_g, fill = factor(yr), group = paste(clip_event, yr))) +
  geom_boxplot(alpha = 0.5, col = "grey50") +
  geom_point(aes(shape = fxnl_grp, group = paste(clip_event, yr)), size = 1.5, position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.3, )) +
  facet_grid(nut_trt~ppt_trt) +
  scale_shape_manual(values= c("F", "G")) +
  labs(x = "Annual clip event", y = "ANPP (g)",
       title = "USDA Compost 2019-2021, quickplot for QA",
       subtitle = "Plots 4 (NW) and 11 (CD) missing for 2020 event 2 (only 3 data points)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))
ggsave(paste0(datpath, "ANPP/ANPP_Figures/ANPP_prelimQA.pdf"), width = 6, height = 4, scale = 1.2)

# add in records for missing plots to indicate they are missing
missing_recs <- subset(anpp_master, yr == 2020 & clip_event == 2 & fulltrt %in% c("NW", "CD")) %>%
  group_by(fulltrt) %>%
  filter(block == min(block)) %>%
  ungroup() %>%
  # edit recs to create observations for plots 4 and 11
  mutate(plot = ifelse(fulltrt == "CD", 11, 4)) %>%
  # drop block and plotid then add back in to be sure
  dplyr::select(-c(block, plotid)) %>%
  left_join(trtkey) %>%
  # na most cols
  mutate_at(vars(file:line, dry_wgt_g:weigh_init), function(x) x <-NA) %>%
  dplyr::select(names(anpp_master)) %>%
  # add QA note to notes col
  mutate(notes = "QA note: not collected, missing")

# rbind to master and re-sort
anpp_master <- rbind(anpp_master, missing_recs) %>%
  arrange(yr, date_clip, page, line)
# one last check
summary(anpp_master) # looks okay



# -- FINISHING ----
# write out dataset to cleaned folder
write_csv(anpp_master, paste0(datpath, "ANPP/ANPP_CleanedData/Compost_ANPP_clean.csv"))

