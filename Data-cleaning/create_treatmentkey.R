# create compost treatment lookup table
# author(s): ctw
# created: april 2019

# note: update table when add in seeding treatment in fall 2019


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
# add plot name: concatenate block, nutrient treatment, and precipitation treatment
trt_df$plot <- paste0(trt_df$block, trt_df$nut_trt, trt_df$ppt_trt)
# move plot to first column
trt_df <- trt_df[c("plot", "block", "nut_trt", "ppt_trt")]


# -- FINISHING ----
# set path to compost data (main level)
# dependent on user, comment out path(s) that aren't pertinent to you
datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/" #ctw's path
#datpath <- "~/Dropbox/USDA-compost/Data/" # should work for LMH and AS

# write out to compost dropbox main data folder level
write.csv(trt_df, paste0(datpath, "Compost_TreatmentKey.csv"), row.names = F)

             