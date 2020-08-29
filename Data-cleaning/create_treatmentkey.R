# create compost treatment lookup table
# author(s): ctw
# created: april 2019
# last update: aug 2020



# -- SETUP -----
# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/"
}

# read in updated treatment key in site-map dropbox folder (created in sep 2019)
trtkey <- read.csv(paste0(datpath, "site-map/treatment_plots.csv"), stringsAsFactors = F)



# -- CREATE KEY -----
# setup table of block, nutrient ammendment treatment, and precipitation treatment:
## 4 blocks: 1, 2, 3, 4
block <- c(rep(1,9), rep(2,9), rep(3,9), rep(4,9))
## 3 nutrient treatments: +Compost ("C"), +Fertilizer ("F"), None ("N")
nut_trt <- rep(c(rep("C",3), rep("F", 3), rep("N", 3)), 4)
## 3 precip treatments: Drought [shelters] ("D"), Wet [sprinklers] ("W"), Control ("XC")
ppt_trt <- rep(rep(c("D", "W", "XC"),3),4)
# put together
trt_df <- data.frame(block, nut_trt, ppt_trt)
# add full treatment: concatenate nutrient treatment and precipitation treatment
trt_df$plotid <- paste0(trt_df$block, trt_df$nut_trt, trt_df$ppt_trt)
# add plot id: concatenate block, nutrient treatment, and precipitation treatment
trt_df$fulltrt <- paste0(trt_df$nut_trt, trt_df$ppt_trt)
# add in plot numbers
trt_df <- merge(trt_df, trtkey)
# order by plotnum
trt_df <- trt_df[order(trt_df$plot),]
# move plot to first column
trt_df <- trt_df[c("plot", "plotid", "fulltrt", "block", "nut_trt", "ppt_trt")]



# -- FINISHING ----
# write out to compost dropbox main data folder level
write.csv(trt_df, paste0(datpath, "Data/Compost_TreatmentKey.csv"), row.names = F)

             