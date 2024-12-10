# supplemental figures
# ctw

# notes:
# just making rain summary for SLM presentation 2022 Apr 5 for now
# more later...


# -- SETUP ----
# load needed libraries
library(tidyverse)
library(lubridate)
library(SPEI)
library(cimir)
library(cowplot)
# for stats
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(broom)

# modify default settings
options(stringsAsFactors = F)
theme_set(theme_test())
na_vals <- c("" , " ","NA", NA)
# plot theme for ggplots
theme_ctw <- function(){
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
}

# specify dropbox pathway (varies by user)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# climvar datpath for native experiment ther
climvar_path <- "~/Dropbox/ClimVar/DATA/"

# raw cimis data for daily ppt and temp summaries over time
# > there is some temporal overlap here so can stitch together
compost_cimis <- read.csv(paste0(datpath, "CIMIS/Download_RawData/CIMIS_084BrownsValley_hourly_20211021.csv"), blank.lines.skip = T, na.strings = na_vals)
climvar_cimis <- read.csv(paste0(climvar_path, "Met\ Station\ CIMIS/cimis_brownsvalley_daily_20190320.csv"), blank.lines.skip = T, na.strings = na_vals)

# read in cleaned soil moisture and rainfall prepped for soil moisture
compost_sm <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMoisture_all_clean.csv"), na.strings = na_vals)
compost_sm_ppt <- read.csv(paste0(datpath, "/SoilMoisture/SoilMoisture_DataQA/CIMIS_084BrownsValley_pptforSoilMoisture.csv"), na.strings = na_vals)
str(compost_sm)

# try getting cimis data via API
# > doesn't work using browns valley station id, so just use SFREC zip code
cimis_2010 <- cimis_data(targets = 95918, start.date = as.Date("2009-01-01"), end.date = as.Date("2009-01-01")+1749, measure.unit = "M", 
                    items = c("day-precip", "day-air-tmp-max", "day-air-tmp-min", "day-air-tmp-avg"))
cimis_2016 <- cimis_data(targets = 95918, start.date = max(cimis_2010$Date), end.date =  max(cimis_2010$Date)+1749, measure.unit = "M", 
                         items = c("day-precip", "day-air-tmp-max", "day-air-tmp-min", "day-air-tmp-avg"))
cimis_2021 <- cimis_data(targets = 95918, start.date = max(cimis_2016$Date), end.date =  max(cimis_2016$Date)+1749, measure.unit = "M", 
                         items = c("day-precip", "day-air-tmp-max", "day-air-tmp-min", "day-air-tmp-avg"))
cimis_2023 <- cimis_data(targets = 95918, start.date = max(cimis_2021$Date), end.date =  Sys.Date()-1, measure.unit = "M", 
                         items = c("day-precip", "day-air-tmp-max", "day-air-tmp-min", "day-air-tmp-avg"))
all_cimis <- rbind(cimis_2010, cimis_2016, cimis_2021,cimis_2023) %>%
  distinct()
cimis_missing_dates <- unique(all_cimis$Date[is.na(all_cimis$Value)])
# try missing data at station 195 (auburn, ca)
station_195.2009 <- cimis_data(targets = 195, start.date = as.Date("2009-12-01"), end.date = as.Date("2009-12-31"), measure.unit = "M", 
                            items = c("day-precip", "day-air-tmp-max", "day-air-tmp-min", "day-air-tmp-avg"))
station_195.1 <- cimis_data(targets = 195, start.date = cimis_missing_dates[2], end.date = cimis_missing_dates[7], measure.unit = "M", 
                          items = c("day-precip", "day-air-tmp-max", "day-air-tmp-min", "day-air-tmp-avg"))
# second period missing
station_195.2 <- cimis_data(targets = 195, start.date = cimis_missing_dates[8], end.date = cimis_missing_dates[length(cimis_missing_dates)], measure.unit = "M", 
                          items = c("day-precip", "day-air-tmp-max", "day-air-tmp-min", "day-air-tmp-avg"))

# stack station 195
station_195 <- rbind(station_195.2009, station_195.1, station_195.2)

# -- DATA PREP ----
# want daily ppt and temp to show trend over time
# maybe also make spei
# add water year to raw cimis dats

# infill missing data in cimis dat from api
# change station 195 value name to differentiate from bv cimis
station_195 <- rename(station_195, "Value_195" = "Value")
all_cimis <- left_join(all_cimis, station_195[c("Date", "Julian", "Item","Value_195")]) %>%
  mutate(infill_values = ifelse(!is.na(Value), Value, Value_195))
# plot to see how it looks
ggplot(all_cimis, aes(Date, infill_values)) +
  geom_line() +
  geom_point(data = subset(all_cimis, is.na(Value)), col = "purple") +
  facet_wrap(~Item, scales = "free_y") # seems ok. tmin is super weird for bv pre-2016
summary(all_cimis) # no more missing values using auburn-infilled data

make_wy <- function(dat, datecol, datform){
  # remove any blank rows
  dat <- dat[!apply(dat,1, function(x) all(is.na(x))),]
  dat[[datecol]] <- as.Date(dat[[datecol]], format = datform)
  # add yr, mon, day
  dat$yr <- year(dat[[datecol]])
  dat$mon <- month(dat[[datecol]])
  dat$doy <- yday(dat[[datecol]])
  
  # make waterYear and dowy
  dat$waterYear <- with(dat, ifelse(mon >= 10, yr+1, yr))
  dat$wymon <- with(dat, ifelse(mon >= 10, mon-9, mon+3))
  # add dowy by whether leap yr (no remainder dividin calyr by 4)
  dat <- group_by(dat, yr) %>%
    mutate(dowy = ifelse(mon %in% 10:12 & yr%%4 != 0, doy-273, #nonleap (365 days/yr)
                         ifelse(mon %in% 10:12 & yr%%4 == 0, doy-274, # leap yr (366 days/yr)
                                doy+92))) # all else
}



# be sure timestamps are correctly classed
climvar_cimis <- make_wy(climvar_cimis, "Date", "%m/%d/%Y")
compost_cimis <- make_wy(compost_cimis, "Date", "%m/%d/%Y")
all_cimis <- make_wy(all_cimis, datecol = "Date", datform = "%Y-%m-%d")

# aggregate compost to daily
compost_cimis_daily <- group_by(compost_cimis, Date) %>%
  summarise(ppt.mm = sum(Precip..mm., na.rm = T),
         tmax.c = max(Air.Temp..C., na.rm = T),
         tmin.c = min(Air.Temp..C., na.rm = T),
         tmean.c = mean(Air.Temp..C., na.rm = T)) %>%
  ungroup() %>%
  # if anything infinite or NaN, make NA
  mutate_if(names(.) != "Date", function(x) ifelse(is.infinite(x), NA, 
                                ifelse(is.nan(x), NA, x))) %>%
  make_wy(., "Date", "%m/%d/%Y") %>%
  mutate(source = "compost")


# subset just dat cols and temp mets of interest
sfrec_climate <- climvar_cimis %>%
  subset(select = names(.)[grep("Date|yr|mon|^wat|doy|Air.Temp|Preci", names(.))]) %>%
  rename_if(grepl("^Min", names(.)), function(x) x <- "tmin.c") %>%
  rename_if(grepl("^Avg", names(.)), function(x) x <- "tmean.c") %>%
  rename_if(grepl("^Max", names(.)), function(x) x <- "tmax.c") %>%
  rename_if(grepl("^Pre", names(.)), function(x) x <- "ppt.mm") %>%
  mutate(source = "climvar") %>%
  rbind(subset(compost_cimis_daily, select = names(.))) %>%
  group_by(Date) %>%
  mutate(nobs = length(source)) %>%
  # overlaps 2018-07-01 to 2019-03-19, ppt same but temps <1 degree diff. unsure why and not important for trends fig, choose compost
  filter(nobs == 1 | (nobs == 2 & source == "compost")) %>%
  ungroup()

# need to infill temp NAs for SPEI to run
sfrec_climate <- mutate(sfrec_climate, test = ifelse(is.na(tmin.c), zoo::rollmean(tmin.c, 7, fill = NA, align = "left"), tmin.c))

# make 1mo SPEI (coords of SFREC: 39.25170524042971, -121.31351343863173)
sfrec_PET <- thornthwaite(sfrec_climate$tmean.c[sfrec_climate$yr < 2016], lat = 39.25170524042971, na.rm = T)
sfrec_BAL <- sfrec_climate$ppt.mm[sfrec_climate$yr < 2016] - sfrec_PET
sfrec_SPEI.1 <- spei(sfrec_BAL, 1, na.rm = T)
plot(sfrec_SPEI.1) # missing data are making full spei calculations funky

# make 1mo SPEI (coords of SFREC: 39.25170524042971, -121.31351343863173)
# make monthly times series first
cimis_monthly <- group_by(all_cimis, mon, yr, waterYear, wymon) %>%
  mutate(mon_ppt = sum(infill_values[Item == "DayPrecip"]),
         mon_tmin = min(infill_values[Item == "DayAirTmpMin"]),
         mon_tmax = max(infill_values[Item == "DayAirTmpMax"]),
         mon_tmean = mean(infill_values[Item == "DayAirTmpAvg"]),
         nobs = length(infill_values[Item == "DayAirTmpAvg"])) %>%
  ungroup() %>%
  subset(select = c(mon, yr, waterYear, wymon, mon_ppt:ncol(.))) %>%
  distinct() %>%
  # drop anything that's not a full month
  subset(!nobs < 28)
# check first and final month-year available
cimis_monthly[1,] # 2009-01
cimis_monthly[nrow(cimis_monthly),] # 2024-11
# create time series for SPEI
cimis_ts_tmean <- ts(all_cimis$infill_values[all_cimis$Item == "DayAirTmpAvg"], start = c(2009,01), freq = 365)
cimis_PET <- thornthwaite(ts(cimis_monthly$mon_tmean, 
                             start = c(2009,01), 
                             end = c(2024, 11), 
                             freq = 12), lat = 39.25170524042971)
cimis_BAL <- ts(cimis_monthly$mon_ppt, start = c(2009, 01), end = c(2024, 11), freq=12) - cimis_PET
cimis_SPEI.1 <- spei(cimis_BAL, 1)
cimis_SPEI.3 <- spei(cimis_BAL, 3)
plot(cimis_SPEI.1)

ggplot(all_cimis, aes(Date, infill_values)) +
  geom_line() +
  facet_wrap(~Item, scales = "free_y")


# detrend temp
tmean_trend <- decompose(cimis_ts_tmean)
tmax_trend = decompose(ts(all_cimis$infill_values[all_cimis$Item == "DayAirTmpMax"], start = c(2009,01), freq = 365))
tmin_trend = decompose(ts(all_cimis$infill_values[all_cimis$Item == "DayAirTmpMin"], start = c(2009,01), freq = 365))
ppt_trend = decompose(ts(all_cimis$infill_values[all_cimis$Item == "DayPrecip"], start = c(2009,01), freq = 365))
ppt_trend_monthly <- decompose(ts(cimis_monthly$mon_ppt, start = c(2009,01), freq = 12))
plot(ppt_trend_monthly)
plot(tmax_trend)

# make daily temp trend
temptrends <- data.frame(date = unique(all_cimis$Date),
                         tmax = as.numeric(tmax_trend$trend),
                         tmax_deseasoned = as.numeric(tmax_trend$trend + tmax_trend$random),
                         tmin = as.numeric(tmin_trend$trend),
                         tmin_deseasoned = as.numeric(tmin_trend$trend + tmin_trend$random),
                         tmean = as.numeric(tmean_trend$trend),
                         tmean_deseasoned = as.numeric(tmean_trend$trend + tmean_trend$random)) %>%
  gather(temp, trend, tmax:ncol(.))
ggplot(temptrends, aes(date, trend, col = temp)) +
  geom_line()

# look at trend -- calculate for 2010-2024 WY
summary(lm(trend ~ date, data = subset(temptrends, temp == "tmin" & date >= as.Date("2009-10-01") & date < as.Date("2024-10-01")))) # tmin and tmean are both about +1deg per decade like tmax
summary(lm(trend ~ date, data = subset(temptrends, temp == "tmean" & date >= as.Date("2009-10-01") & date < as.Date("2024-10-01"))))
summary(lm(trend ~ date, data = subset(temptrends, temp == "tmax"& date >= as.Date("2009-10-01") & date < as.Date("2024-10-01"))))
# tmin: date  1.707e-04  4.275e-06   39.92   <2e-16 ***
# tmean: date  1.588e-04  5.358e-06   29.64   <2e-16 ***
# tmax: date  1.541e-04  6.594e-06   23.38   <2e-16 ***

1.707e-04 * 365 *10  # tmin increase at 0.623 C per decade for WY 2010-2024
1.588e-04 * 365 *10  # tmean increase at 0.579 C per decade
1.541e-04 * 365 *10  # tmax increase at 0.562 C per decade

with(subset(temptrends, !is.na(trend)), range(date))
temptrends_datespres <- subset(temptrends, !is.na(trend))
max(temptrends_datespres$date) - min(temptrends_datespres$date)


# -- build plots ----
# tmin, tmean, and tmin trends
tmin_plot <- ggplot() +
  #geom_line(data = subset(all_cimis, grepl("Min", Item)), aes(Date, Value), alpha = 0.8, col = "lightblue") +
  geom_line(data = subset(temptrends, temp == "tmin_deseasoned"), aes(date, trend), col= "lightblue", alpha = 0.7) +
  geom_line(data = subset(temptrends, temp == "tmin"), aes(date, trend), col= "dodgerblue2", lwd = 1, alpha = 0.7) +
  labs(y = "Daily Tmin (°C)", x = "Date") +
  theme(plot.margin = margin(0,0,0,0, "pt"))
tmean_plot <- ggplot() +
  #geom_line(data = subset(all_cimis, grepl("Avg", Item)), aes(Date, Value), alpha = 0.8, col = "orchid3") +
  geom_line(data = subset(temptrends, temp == "tmean_deseasoned"), aes(date, trend), col= "orchid3", alpha = 0.7) +
  geom_line(data = subset(temptrends, temp == "tmean"), aes(date, trend), col= "purple4", lwd = 1, alpha = 0.8) +
  labs(y = "Daily Tmean (°C)", x = NULL) +
  theme(axis.text.x = element_blank(),
        plot.margin = margin(0,0,0,0, "pt"))
tmax_plot <- ggplot() +
  #geom_line(data = subset(all_cimis, grepl("Max", Item)), aes(Date, Value), alpha = 0.8, col = "indianred1") +
  geom_line(data = subset(temptrends, temp == "tmax_deseasoned"), aes(date, trend), col= "indianred1", alpha = 0.7) +
  geom_line(data = subset(temptrends, temp == "tmax"), aes(date, trend), col= "coral4", lwd = 1, alpha = 0.8) +
  labs(y = "Daily Tmax (°C)", x = NULL) +
  theme(axis.text.x = element_blank(),
        plot.margin = margin(0,0,0,0, "pt"))
plot_grid(tmax_plot, tmean_plot, tmin_plot, nrow = 3, align = "hv")

# remmake it ggplot panel
ggtemp <- separate(temptrends, temp,into = c("temp", "type"), sep = "_", remove = T) %>%
  mutate(type = ifelse(!is.na(trend) & is.na(type), "trend", type),
         temp = gsub("t", "T", temp)) 

ggtemp_fig <- ggplot(subset(ggtemp, !is.na(trend) & date >= as.Date("2009-10-01") & date <as.Date("2024-10-01")), aes(date, trend)) +
  geom_vline(aes(xintercept = as.Date("2018-10-01")), lty = 2) +
  geom_ribbon(aes(xmin = as.Date("2020-10-01"), xmax = as.Date("2021-06-01")), fill = "grey80") +
  geom_line(aes(linewidth = type, group = type, col = type, alpha = type)) +
  # study period
  geom_vline(aes(xintercept = as.Date("2020-10-01")), lty = 2) +
  geom_vline(aes(xintercept = as.Date("2021-06-01")), lty = 2) +
  scale_color_manual(name = "Temperature", 
                     values = c("black", "red")) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  scale_x_date(expand = c(0.01,0.01), date_breaks = "2 years", date_labels = "%Y") +
  scale_linewidth_discrete(name = "Temperature", range = c(0.5,1.5)) +
  scale_alpha_discrete(range= c(0.8, 1)) +
  labs(y = "Daily Temperature (°C)", x = NULL) +
  theme_ctw() +
  facet_wrap(~temp, nrow = 3, scales = "free_y") +
  theme(legend.position = "none")

monthdat <- cbind(cimis_monthly, 
                  spei3 = as.numeric(cimis_SPEI.3$fitted),
                  ppt_trend = as.numeric(ppt_trend_monthly$trend)
                  ) %>%
  mutate(date = as.Date(paste(yr, mon, 01, sep = "-"))) %>%
  left_join(distinct(all_cimis[c("Date", "doy", "dowy")]), by = c("date"= "Date")) %>%
  distinct()

speiplot <- monthdat  %>% #[c("date", "spei3")] %>%
  mutate(nondrought = ifelse(spei3>=0, spei3, 0),
         drought = ifelse(spei3<0, spei3, 0)) %>%
  subset(waterYear < 2025 & waterYear > 2009) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), lty = 1) +
  # add ribbon of study period
  geom_ribbon(aes(y = spei3, xmin = as.Date("2020-10-01"), xmax = as.Date("2021-06-01")), fill = "grey80") + 
  # add vlines at experiment start and end
  geom_vline(aes(xintercept = as.Date("2018-10-01")), lty = 2, col = "black") +
  geom_area(aes(date, nondrought), fill = "royalblue1", col = "grey30", alpha = 0.9) +
  geom_area(aes(date, drought), fill = "salmon3", col = "grey30", alpha = 0.9) +
  # study period
  geom_vline(aes(xintercept = as.Date("2020-10-01")), lty = 2, col = "black") +
  geom_vline(aes(xintercept = as.Date("2021-06-01")), lty = 2, col = "black") +
  scale_x_date(expand = c(0.01,0.01), date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  labs(y = "3-month SPEI", x = NULL)+
  #theme(plot.margin = margin(0,0,0,0, "pt")) +
  theme_ctw()

pptplot <- ggplot(subset(monthdat, waterYear < 2025 & waterYear > 2009)) +
  # add ribbon of study period
  geom_ribbon(aes(date, mon_ppt, xmin = as.Date("2020-10-01"), xmax = as.Date("2021-06-01")), fill = "grey80") + 
  # vert lines to denote full study period
  geom_vline(aes(xintercept = as.Date("2018-10-01")), lty = 2, col = "black") +
  geom_col(aes(date, mon_ppt), fill = "deepskyblue1") +
  #geom_point(aes(date, mon_ppt), pch = 21, col = "slateblue2", fill = "deepskyblue2", alpha = 0.7, size = 3) +
  geom_line(aes(date, ppt_trend), col = "steelblue", lwd = 1.5) +
  # study period
  geom_vline(aes(xintercept = as.Date("2020-10-01")), lty = 2, col = "black") +
  geom_vline(aes(xintercept = as.Date("2021-06-01")), lty = 2, col = "black") +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  scale_x_date(expand = c(0.01,0.01), date_breaks = "2 years", date_labels = "%Y") +
  labs(y = "Monthly precip (mm)", x = NULL) +
  theme_ctw() +
  theme(axis.text.x = element_blank()
        #plot.margin = margin(0,0,0,0, "pt")
        )


plot_grid(plot_grid(ggtemp_fig, labels= "A", label_size = 18), 
          plot_grid(pptplot, speiplot, nrow = 2, 
                    axis = "l",
                    align = "vh", 
                    labels = c("B", "C"),
                    label_size = 18), 
          ncol = 2)

# plot for esa 2023
monthdat[c("date", "spei3")] %>%
  mutate(nondrought = ifelse(spei3>=0, spei3, 0),
         drought = ifelse(spei3<0, spei3, 0)) %>%
  #subset(!is.na(spei3)) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), lty = 1) +
  # add vlines at experiment start and end
  geom_ribbon(aes(y = (spei3), xmin = as.Date("2020-10-01"), xmax = as.Date("2021-06-01")), fill = "grey80", alpha = 0.5) +
  # geom_vline(aes(xintercept = as.Date("2018-10-01")), lty = 1, col = "grey50" ) +
  geom_vline(aes(xintercept = as.Date("2020-10-01")), lty = 2, col = "grey50") +
  geom_vline(aes(xintercept = as.Date("2021-06-01")), lty = 2, col = "grey50") +
  geom_area(aes(date, nondrought), fill = "royalblue1", col = "grey30", alpha = 0.8) +
  # geom_vline(aes(xintercept = as.Date("2028-10-01")), lty = 2, col = "grey30") +
  # geom_vline(aes(xintercept = as.Date("2021-06-01")), lty = 2, col = "grey30") +
  scale_x_date(expand = c(0.005,0.005), date_labels = "%Y", 
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2023-01-01"), by = (365.25*4)) ) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  geom_area(aes(date, drought), fill = "salmon3", col = "grey30", alpha = 0.8) +
  labs(caption = "data source: CIMIS", 
       y = "3-month SPEI", x = NULL, 
       subtitle = "Browns Valley, California") +
  theme(axis.text = element_text(size = 11),
        plot.caption = element_text(size = 10),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12))

# write out for esa 2023
# figpath <- "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/esa2023_presfigs/"
# ggsave(paste0(figpath,"sfrec_3monspei_20092023.pdf"), units = "in", width = 7, height = 3)
# 
# # no experiment timeline markings
# ggsave(paste0(figpath,"sfrec_3monspei_20092023_nomarks.pdf"), units = "in", width = 7, height = 3)
# # native exp only
# ggsave(paste0(figpath,"sfrec_3monspei_20092023_natexo_marks.pdf"), units = "in", width = 7, height = 3)



# -- SOIL MOISTURE v. PREPPED PPT ------
sm_summary <- subset(compost_sm, waterYear == 2021 & !mon %in% 7:9) %>%
  # make any sub 0 vwc = 0
  mutate(vwc = ifelse(vwc <0, 0, vwc)) %>%
  group_by(cleanorder, clean_datetime, fulltrt, nut_trt, ppt_trt) %>%
  summarise(mean_vwc = mean(vwc, na.rm = T),
            nobs = length(vwc[!is.na(vwc)]),
            se = sd(vwc, na.rm = T)/sqrt(nobs),
            .groups = "drop_last") %>%
  ungroup() %>%
  # convert NaN to NA
  mutate(mean_vwc = ifelse(is.nan(mean_vwc), NA, mean_vwc),
         upper = mean_vwc + se,
         lower = mean_vwc - se,
         ppt_trt = factor(ppt_trt, levels = c("W", "XC", "D")),
         clean_datetime = as.POSIXct(clean_datetime, tz = "UTC"))

meanvwc_plot <- ggplot(sm_summary, aes(clean_datetime, mean_vwc, col = nut_trt)) +
  geom_line(lwd = 0.75) +
  #geom_line(aes(cleanorder, upper, lty = nut_trt)) +
  #geom_line(aes(cleanorder, lower, lty = nut_trt)) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  facet_wrap(~ppt_trt, nrow = 3) +
  labs(x = NULL, y = "Mean VWC") +
  scale_color_brewer(name = NULL, palette = "Accent") +
  theme(legend.position = c(0.95, 0.97),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.key.height = unit(1, "pt")) +
  theme_ctw()

vwc_ppt_plot <-  ggplot(subset(compost_sm_ppt, waterYear == 2021 & !mon %in% 7:9), 
       aes(as.POSIXct(clean_datetime, tz = "UTC"), ppt_mm)) +
  geom_col(col = "dodgerblue3", lwd = 0.75) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  labs(x = "Native restoration study (Oct 2020 - Jun 2021)", y = "Precip (mm)") +
  theme_ctw()

plot_grid(meanvwc_plot, vwc_ppt_plot, nrow = 2, 
          axis = "l",
          align = "hv", rel_heights = c(1,0.5))

# trying adding difference for WC back when drops to 1 sensor. these are average values across all sensors available, and the sensor available was dry
# 2024-12-03: this is a holdover until I have time to correct properly
sm_summary_wc <- subset(sm_summary, nut_trt == "C" & ppt_trt == "W") %>%
  mutate(mon = month(clean_datetime), 
         diff = mean_vwc - lag(mean_vwc),
         diff_nobs = nobs - lag(nobs))
subset(sm_summary_wc, diff_nobs != 0) # these are the dates where it gets off

# manually code corrections
fix1 <- c(10248:11046) # 0.06
fix2 <- c(11087:11112) # 0.07
fix3 <- c(11128:11144) # 0.07
fix4 <- c(11405:max(sm_summary_wc$cleanorder)) # 0.05

sm_summary_wc$mean_vwc_adj <- sm_summary_wc$mean_vwc
sm_summary_wc$mean_vwc_adj[sm_summary_wc$cleanorder %in% fix1] <- sm_summary_wc$mean_vwc[sm_summary_wc$cleanorder %in% fix1] + 0.06
sm_summary_wc$mean_vwc_adj[sm_summary_wc$cleanorder %in% c(fix2, fix3)] <- sm_summary_wc$mean_vwc[sm_summary_wc$cleanorder %in% c(fix2, fix3)] + 0.07
sm_summary_wc$mean_vwc_adj[sm_summary_wc$cleanorder %in% fix4] <- sm_summary_wc$mean_vwc[sm_summary_wc$cleanorder %in% fix4] + 0.05

ggplot(sm_summary_wc, aes(clean_datetime, mean_vwc, col = nobs)) +
  geom_line() +
  geom_line(aes(clean_datetime, mean_vwc_adj), col = "pink") # this seems better, there is a still a sharp drop in spring that seems suspicious

# calculate diffs again
sm_summary_wc$diff2 <- with(sm_summary_wc, mean_vwc_adj - lag(mean_vwc_adj))

ggplot(subset(sm_summary_wc, mon %in% 4:5), aes(clean_datetime, mean_vwc, col = nobs)) +
  geom_line() +
  geom_line(aes(clean_datetime, mean_vwc_adj), col = "pink") # seems like its around april 15 to late may

subset(sm_summary_wc, diff2 == min(diff2, na.rm = T)) # april 16, dropped 0.03, dropped -0.0348
subset(sm_summary_wc, mon == 5 & diff2 == max(diff2[mon == 5], na.rm = T)) # to may 28, gained 0.0127
mean(c(0.0348, 0.0127)) # 0.023

ggplot(subset(sm_summary_wc, cleanorder %in% c(10500:11100)), aes(clean_datetime, mean_vwc, col = nobs)) +
  geom_line() +
  geom_line(aes(clean_datetime, mean_vwc_adj), col = "pink") +
  geom_line(data = subset(sm_summary_wc, cleanorder %in% c(10543:11046)), aes(clean_datetime, mean_vwc_adj+ 0.023), col = "orange") # add mean diff

sm_summary_wc$mean_vwc_adj[sm_summary_wc$cleanorder %in% c(10543:11046)] <- sm_summary_wc$mean_vwc[sm_summary_wc$cleanorder %in% c(10543:11046)] + 0.06 + 0.015 #0.023
# > switched to more modest adjustment. 
# add back to sm_summary
sm_summary_adj <- subset(sm_summary, fulltrt != "CW") %>%
  mutate(mean_vwc_adj = mean_vwc) %>%
  rbind(sm_summary_wc[names(.)])

# remake plot with adjusted values
meanvwc_plot_adj <- ggplot(sm_summary_adj, aes(clean_datetime, mean_vwc_adj, col = nut_trt)) +
  geom_line(lwd = 0.75) +
  #geom_line(aes(cleanorder, upper, lty = nut_trt)) +
  #geom_line(aes(cleanorder, lower, lty = nut_trt)) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  facet_wrap(~ppt_trt, nrow = 3) +
  labs(x = NULL, y = "Mean VWC") +
  scale_color_brewer(name = NULL, palette = "Accent") +
  theme(legend.position = c(0.99, 0.99),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.key.height = unit(1, "pt")) +
  theme_ctw()
meanvwc_plot_adj


# -- CHAPTER FIGS ----
# show it by nut_trt instead of ppt
meanvwc_plot_adjnut <- mutate(sm_summary_adj, nut_trt = factor(nut_trt, levels = c("C", "F", "N"), labels = c("Compost", "Fertilizer", "Ambient (Control)"))) %>%
  ggplot(aes(clean_datetime, mean_vwc_adj, col = ppt_trt)) +
  geom_line(lwd = 0.75) +
  #geom_line(aes(cleanorder, upper, lty = nut_trt)) +
  #geom_line(aes(cleanorder, lower, lty = nut_trt)) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  facet_wrap(~nut_trt, nrow = 3, labeller = ) +
  labs(x = NULL, y = "Mean VWC") +
  #scale_color_brewer(name = NULL, palette = "Blues", direction = -1) +
  scale_color_manual(name = NULL,values = c("W" = "blue2", "XC" = "steelblue2", "D" = "darkgoldenrod4"), 
                     labels = c("W" = "Wet", "XC" = "Control", "D"= "Drought")) +
  theme(legend.position = c(0.99, 0.99),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.key.height = unit(1, "pt")) +
  theme_ctw()
meanvwc_plot_adjnut

# make temp plot for the year
study_temp <- subset(all_cimis, waterYear == 2021 & grepl("Air", Item)) %>%
  mutate(Item = gsub("DayAirTmp", "T", Item),
         Item = factor(Item, levels = c("TMax", "TAvg", "TMin"), labels = c("Tmax", "Tmean", "Tmin"))) %>%
  # calculate rollmean for each
  group_by(Item) %>%
  mutate(rolling = zoo::rollmean(infill_values, k = 5, na.pad = T))

study_temp_fig <- subset(study_temp, !mon %in% 7:9) %>%
  ggplot(aes(Date, infill_values,  col = Item)) +
  geom_point(pch = 1) +
  geom_line(aes(Date, rolling), lwd = 0.75) +
  geom_smooth(aes(fill = Item), show.legend = F) +
  labs(y = "Daily Temp (°C)", x = NULL) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  scale_color_manual(name = NULL, values = c("indianred3", "mediumpurple", "steelblue3")) +
  scale_fill_manual(name = NULL, values = c("indianred1", "mediumpurple1", "steelblue1")) +
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme_ctw() +
  theme(legend.position = "inside",
        legend.direction = "horizontal",
        legend.justification.inside = c(0.01,0.99))

# make study SPEI fig
glimpse(monthdat[c("date", "spei3")])
study_spei <- monthdat[c("date", "spei3")] %>%
  filter(date >= as.Date("2020-10-01") & date < as.Date("2021-07-01")) %>%
  ggplot(aes(date, spei3)) +
  geom_hline(aes(yintercept = 0)) +
  geom_area(xmin = as.Date("2020-10-01"), xmax = as.Date("2021-07-01"), fill = "sienna") +
  scale_x_date(expand = c(0,0))

plot_grid(study_temp_fig, meanvwc_plot_adjnut, vwc_ppt_plot, nrow = 3, 
          axis = "l",
          align = "hv", rel_heights = c(0.4, 1.2,0.4))

plot_grid(plot_grid(study_temp_fig, vwc_ppt_plot, nrow = 2, 
          axis = "l",
          align = "hv"), meanvwc_plot_adjnut, ncol = 2)


# -- compost VWC check. plot all years to see if/when effect attenuated ----
sm_summary_all <- subset(compost_sm) %>%
  # make any sub 0 vwc = 0
  mutate(vwc = ifelse(vwc <0, 0, vwc)) %>%
  group_by(cleanorder, clean_datetime, waterYear, fulltrt, nut_trt, ppt_trt) %>%
  reframe(mean_vwc = mean(vwc, na.rm = T),
            nobs = length(vwc[!is.na(vwc)]),
            se = sd(vwc, na.rm = T)/sqrt(nobs)) %>%
  #          .groups = "drop_last") %>%
  #ungroup() %>%
  # convert NaN to NA
  mutate(mean_vwc = ifelse(is.nan(mean_vwc), NA, mean_vwc),
         upper = mean_vwc + se,
         lower = mean_vwc - se,
         ppt_trt = factor(ppt_trt, levels = c("W", "XC", "D")),
         nut_trt = factor(nut_trt, levels = c("C", "F", "N"), labels = c("Comp", "Fert", "Control")),
         clean_datetime = as.POSIXct(clean_datetime, tz = "UTC"))

ggplot(sm_summary_all, aes(clean_datetime, mean_vwc, group = fulltrt, col = nut_trt)) +
  geom_line(aes(alpha = nut_trt), lwd = .75) +
  scale_alpha_manual(values = c(1,0.8, 0.8)) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(ppt_trt~paste("WY", waterYear), scales = "free_x") +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  theme_ctw()


ggplot(subset(sm_summary_all, nut_trt == "Comp"), aes(clean_datetime, mean_vwc, group = fulltrt)) +
  geom_line(aes(col = nobs)) +
  geom_point(aes(size = nobs,col = nobs),pch = 1) +
  scale_size_continuous(range = c(8,1)) +
  scale_alpha_manual(values = c(1,0.8, 0.8)) +
  scale_color_viridis_c() +
  facet_grid(ppt_trt~paste("WY", waterYear), scales = "free_x") +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  theme_ctw()

ggplot(subset(sm_summary_all, nut_trt == "Control"), aes(clean_datetime, mean_vwc, group = fulltrt)) +
  geom_line(aes(col = nobs)) +
  geom_point(aes(size = nobs,col = nobs),pch = 1) +
  scale_size_continuous(range = c(8,1)) +
  scale_alpha_manual(values = c(1,0.8, 0.8)) +
  scale_color_viridis_c() +
  facet_grid(ppt_trt~paste("WY", waterYear), scales = "free_x") +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  theme_ctw()

ggplot(subset(sm_summary_all, nut_trt == "Fert"), aes(clean_datetime, mean_vwc, group = fulltrt)) +
  geom_line(aes(col = nobs)) +
  geom_point(aes(size = nobs, col = nobs), pch = 1) +
  scale_size_continuous(range = c(6,1)) +
  scale_alpha_manual(values = c(1,0.8, 0.8)) +
  scale_color_viridis_c() +
  facet_grid(ppt_trt~paste("WY", waterYear), scales = "free_x") +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b", expand = c(0.01, 0.01)) +
  theme_ctw() +
  theme(panel.grid.major.y = element_line(color = "grey70")
        #panel.grid.minor.y = element_line(color = "grey70")
        )


# -- SUPPLEMENTAL STATS ----

# total rainfall oct - may wy 2021
with(subset(compost_sm_ppt, waterYear == 2021 & !mon %in% 6:9), sum(ppt_mm))  # 259.6 mm
with(subset(compost_sm_ppt, !mon %in% 6:9), sapply(split(ppt_mm, waterYear), sum))
#   2019   2020   2021 
# 1004.7  412.5  259.6

totppt_gs <- sort(with(subset(all_cimis, waterYear %in% 2010:2024 & grepl("Precip", Item) & !mon %in% 6:9), sapply(split(infill_values, waterYear), sum)))
maxtemp_gs <- sort(with(subset(all_cimis, waterYear %in% 2010:2024 & grepl("TmpMax", Item) & !mon %in% 6:9), sapply(split(infill_values, waterYear), max)))
meantemp_gs <- sort(with(subset(all_cimis, waterYear %in% 2010:2024 & grepl("TmpMax", Item) & !mon %in% 6:9), sapply(split(infill_values, waterYear), mean)))
sort(totppt_gs)  
median(totppt_gs)  
sort(maxtemp_gs) 
sort(meantemp_gs) 

# for comparing VWC, cumulative VCW I don't think will be possible since sensors went offline at different times
# can I compare daily VWC Jan - April 15 with glmm on dates where every treatment has at least 1 sensor?
compost_sm_wy2021 <- subset(compost_sm, waterYear == 2021) %>%
  mutate(ppt_trt = factor(ppt_trt, levels = c("XC", "W", "D"), labels = c("Control", "Wet", "Drought")),
         nut_trt = factor(nut_trt, levels = c("N", "C", "F"), labels = c("Control", "Compost", "Fert")),
         clean_datetime = as.POSIXct(clean_datetime, tz = "UTC")) %>%
  # compare just in the window where treatments applied
  mutate(trtwindow = mon %in% c(1:4))
# may need to correct compost data again in CW for comparisons
str(compost_sm_wy2021)
compare_vwcx_all <- lmerTest::lmer(vwc ~ ppt_trt * nut_trt + (1|portid/block), data = subset(compost_sm_wy2021))
compare_vwcx_all_ja <- lmerTest::lmer(vwc ~ ppt_trt * nut_trt + (1|portid/block), data = subset(compost_sm_wy2021, trtwindow))
# failed
compare_vwcx_all_ja <- lmerTest::lmer(vwc ~ ppt_trt * nut_trt + (1|portid) + (1|block), data = subset(compost_sm_wy2021, trtwindow))
summary(compare_vwcx_all_ja)
compare_vwc_all_ja <- lmerTest::lmer(vwc ~ ppt_trt + nut_trt + (1|portid/block), data = subset(compost_sm_wy2021, trtwindow))
# CAN nest if additive model only
summary(compare_vwc_all_ja)

sm_summary_adj_mod <- mutate(sm_summary_adj, 
                         ppt_trt = factor(ppt_trt, levels = c("XC", "W", "D"), labels = c("Control", "Wet", "Drought")),
                         nut_trt = factor(nut_trt, levels = c("N", "C", "F"), labels = c("Control", "Compost", "Fert")),
                         mon = month(clean_datetime)) %>% 
  # join water year info
  left_join(distinct(testsm[c("clean_datetime", "waterYear", "dowy")]))


# want to compare full growing season Nov - May to understand total reduction
# if comparing just daily VWC, can use cleaned up average daily. compare with all ports mod
compare_vwcx_means <- lmerTest::lmer(mean_vwc ~ ppt_trt * nut_trt + (1|nobs) + (1|dowy), 
                                     #contrasts = list(ppt_trt = "contr.sum"),
                                     data = subset(sm_summary_adj_mod, mon %in% c(11:12, 1:5)))
summary(compare_vwcx_means)
summary(multcomp::glht(compare_vwcx_means,
               linfct = c("ppt_trtWet + nut_trtCompost + ppt_trtWet:nut_trtCompost = 0",
                          "ppt_trtWet + nut_trtFert + ppt_trtWet:nut_trtFert = 0",
                          "ppt_trtWet = 0",
                          "nut_trtCompost = 0",
                          "nut_trtFert = 0",
                          "ppt_trtDrought + nut_trtCompost + ppt_trtDrought:nut_trtCompost = 0",
                          "ppt_trtDrought + nut_trtFert + ppt_trtDrought:nut_trtFert = 0",
                          "ppt_trtDrought = 0"
                          )))

plot(emmeans(compare_vwcx_means, trt.vs.ctrl ~ppt_trt | nut_trt)) + theme_ctw()
compare_vwc_means <- lmerTest::lmer(mean_vwc ~ ppt_trt + nut_trt + (1|nobs) + (1|dowy), 
                                     data = subset(sm_summary_adj_mod, mon %in% c(11:12, 1:5)))
summary(compare_vwc_means)
anova(compare_vwcx_means, compare_vwc_means)
vwc_results <- broom.mixed::tidy(compare_vwcx_means)
compare_vwcx_means_glmm <- glmmTMB::glmmTMB(mean_vwc ~ ppt_trt * nut_trt + (1|nobs) + (1|dowy), 
                                     data = subset(sm_summary_adj_mod, mon %in% 1:4))
summary(compare_vwcx_means_glmm)
plotResiduals(compare_vwcx_means_glmm)
compare_vwcx_means_glmm_ord <- glmmTMB(mean_vwc ~ ppt_trt * nut_trt + (1|nobs) + (1|dowy), family = "ordbeta",
                                            data = subset(sm_summary_adj_mod, mon %in% 1:4))
summary(compare_vwcx_means_glmm_ord)
plotResiduals(compare_vwcx_means_glmm_ord)
performance::check_distribution(compare_vwcx_means_glmm) # residuals are 56% normal; response is 53% beta
performance::compare_performance(compare_vwcx_means_glmm, compare_vwcx_means_glmm_ord)

compare_vwc_means_glmm <- glmmTMB::glmmTMB(mean_vwc ~ ppt_trt + nut_trt + (1|nobs) + (1|mon), 
                                            data = subset(sm_summary_adj_mod, mon %in% 1:4))
summary(compare_vwc_means_glmm)
DHARMa::plotQQunif(compare_vwcx_means_glmm_ord)
performance::compare_performance(compare_vwcx_means_glmm, compare_vwcx_means_glmm_ord)
anova(compare_vwcx_all_ja)

compare_vwc_means_lmer <- lmerTest::lmer(mean_vwc ~ ppt_trt * nut_trt + (1|nobs) + (1|mon), 
                                           data = subset(sm_summary_adj_mod, mon %in% 1:4))
summary(compare_vwc_means_lmer)
plot(emmeans(compare_vwc_means_lmer, ~ ppt_trt | nut_trt))
anova(compare_vwc_means_lmer)
performance::compare_performance(compare_vwc_means_lmer, compare_vwcx_means)
plotQQunif(compare_vwc_means_lmer)

# calculate dailies, mean and max
compost_sm_ppt$clean_datetime <- as.POSIXct(compost_sm_ppt$clean_datetime, tz = "UTC")
daily_sm <- mutate(sm_summary_adj_mod, dt = date(clean_datetime)) %>%
  left_join(compost_sm_ppt) %>%
  group_by(dt, mon, ppt_trt, nut_trt) %>%
  reframe(mean_daily_adj = mean(mean_vwc_adj),
          max_daily_adj = max(mean_vwc_adj),
          mean_nobs = mean(nobs),
          max_nobs = max(nobs),
          min_nobs = min(nobs),
          daily_ppt_mm = sum(ppt_mm)) %>%
  group_by(ppt_trt, nut_trt) %>%
  # diff max to track wetup
  mutate(diff_max = round(max_daily_adj - lag(max_daily_adj),5),
         diff_mean =round(mean_daily_adj - lag(mean_daily_adj),5),
         rain = daily_ppt_mm >= 1, # at least 1mm
         # check wetup
         wetup_max = rain & diff_max > 0,
         wetup_mean = rain & diff_mean > 0,
         dry_wetup_max = rain<=0.1 & diff_max > 0.01,
         dry_wetup_mean = rain <=0.1 & diff_mean > 0.01) %>%
  # add dowy
  group_by(ppt_trt, nut_trt) %>%
  mutate(dowy = 1:length(dt),
         wymon = ifelse(mon %in% 10:12, mon-9, mon+3),
         wymonfac = factor(wymon)) %>%
  ungroup()

plot(daily_sm$mean_daily_adj, daily_sm$max_daily_adj) # may use max to represent SM to capture soil wetness that day. can also detect wetup easier

# -- daily vwc model -----
# try daily model -- holdover model until can figure out how to properly model time series (e.g, time autocor, correct family)
hist(scale(daily_sm$mean_daily_adj))
daily_vwc_mod <- glmmTMB(mean_daily_adj ~ nut_trt * ppt_trt + (1|mean_nobs) + (1|dowy), data = daily_sm, family = "tweedie")
plotResiduals(daily_vwc_mod)
summary(daily_vwc_mod)

# for simple interp, gaussian
daily_vwc_mod_norm <- glmmTMB(mean_daily_adj ~ nut_trt * ppt_trt + (1|wymonfac) + (1|mean_nobs) + (1|dowy), data = subset(daily_sm, mon%in% c(11,12,1:5)))
plotResiduals(daily_vwc_mod_norm)
plotConventionalResiduals(daily_vwc_mod_norm) # still not great, but as good as I can get for now
car::Anova(daily_vwc_mod_norm) # all terms signif

# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
#  Response: mean_daily_adj
#  Chisq Df Pr(>Chisq)    
#  nut_trt           60.576  2  7.017e-14 ***
#  ppt_trt         1131.960  2  < 2.2e-16 ***
#  nut_trt:ppt_trt   72.196  4  7.803e-15 ***

summary(daily_vwc_mod_norm)
plot(emmeans(daily_vwc_mod_norm, ~ ppt_trt | nut_trt)) + theme_ctw() + ggtitle("WY2021 Nov-May daily VWC marginal means") +
  labs(x = "Estimated marginal mean VWC (proportion water)")
emvwc <- emmeans(daily_vwc_mod_norm, ~ nut_trt * ppt_trt, type = "response")
contrast(emvwc, "trt.vs.ctrl", simple = list("nut_trt"))
contrast(emvwc, "trt.vs.ctrl", simple = list("ppt_trt", "nut_trt"), combine = T) 
pairs(emvwc)
# how much does it wetup when it rains

testsm <- left_join(sm_summary_adj_mod, compost_sm_ppt)
ggplot(testsm, aes(ppt_mm, mean_vwc_adj, col = mon)) +
  geom_point() +
  facet_grid(mon ~ ppt_trt, scales = "free")
ggplot(daily_sm, aes(daily_ppt_mm, mean_daily_adj, col = mon)) +
  geom_point() +
  facet_grid(mon ~ ppt_trt, scales = "free")


# try to find wetup signature
subset(testsm, mon %in% 2) %>%
  ggplot(aes(clean_datetime, mean_vwc_adj)) +
  #geom_point(aes(color = ppt_mm), size = 2)+
  geom_point(aes(col = log(ppt_mm+0.00001)), size = 2)+
  geom_line(alpha = 0.75)+ 
  #scale_color_viridis_c()+
  scale_color_distiller(name = "ln(ppt)", palette = "Blues", direction = 1) +
  scale_x_datetime(date_labels = "%d", expand = c(0,0)) +
  labs(x = "February 2021", y = "VWC") +
  facet_wrap(nut_trt~ppt_trt) +
  theme_ctw()
# irrigation began mid Dec and ends mid April
# drought shelters deployed late January (got Jan 4 rainfall, but not storms Jan 25)
# drought got 2 storms in feb, no storms april (only 1 happened)
# by may dry for all 
subset(daily_sm, mon %in% 3) %>%
  ggplot(aes(dt, mean_daily_adj)) +
  geom_line() +
  geom_point(aes(col = daily_ppt_mm))+
  facet_wrap(nut_trt~ppt_trt)

# can track consecutive positive difference or use a minimum threshold wetup
testsm2 <- arrange(testsm, nut_trt, ppt_trt, cleanorder) %>%
  group_by(nut_trt, ppt_trt) %>%
  mutate(diff_mean = round(mean_vwc_adj - lag(mean_vwc_adj),5),
         diffmean2 = diff_mean - lag(diff_mean))


subset(testsm2, mon %in% 11) %>% # it seems like if the change is at least 0.01 it will pick up the wetup event but not diurnal fluctuation
  ggplot(aes(clean_datetime, mean_vwc)) +
  #geom_point(aes(color = ppt_mm), size = 2)+
  geom_point(aes(col = log(ppt_mm+0.00001)), size = 2)+
  geom_line(alpha = 0.75)+ 
  #scale_color_viridis_c()+
  scale_color_distiller(name = "ln(ppt)", palette = "Blues", direction = 1) +
  scale_x_datetime(date_labels = "%d", expand = c(0,0)) +
  #scale_y_continuous(breaks = seq(-0.1,1, 0.02)) +
  labs(x = "February 2021", y = "VWC") +
  facet_wrap(nut_trt~ppt_trt) +
  theme_ctw() +
  theme(panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey"))


wetups <- subset(testsm2, diff_mean >= 0.05) %>%
  group_by(nut_trt, ppt_trt, date) %>%
  filter(cleanorder == min(cleanorder)) %>%
  ungroup() %>%
  group_by(nut_trt, ppt_trt) %>%
  # screen for consec cleanorders
  mutate(diff_order = cleanorder - lag(cleanorder))


# -- track wetups -----
# loop through each treatment to mark wetups
track_wetup <- data.frame()
for(f in unique(testsm2$fulltrt)){
  temp <- subset(testsm2, fulltrt == f)
  # track wetups
  temp_run <- list2DF(rle(temp$diff_mean >= 0.001)) %>%
    mutate(runlen = cumsum(lengths),
           start = runlen - lengths + 1) %>%
    group_by(values) %>%
    # assign event names
    mutate(ct = 1:length(lengths),
           ct = ifelse(!values, paste0("dry_",ct), paste("wet_", ct))) %>%
    ungroup()
  # create column
  temp$wetevent <- "" # nothing for first entry
  for(r in 1:nrow(temp_run)){
    temp$wetevent[temp_run$start[r]:temp_run$runlen[r]] <- temp_run$ct[r]
  }
  temp$wetevent[1] <- ifelse(temp$ppt_mm[1] <= 0.03, "dry_1", "wet_1")
  # go through and ignore any even that doesn't have a diff of at least 0.01
  temp <- group_by(temp, wetevent) %>%
    mutate(screen = any(diff_mean > 0.01))
  temp$type <- ifelse(grepl("wet", temp$wetevent) & temp$screen, "wetup", "dry")
  # id peak
  temp <- group_by(temp, wetevent, type) %>%
    mutate(peak = mean_vwc == max(mean_vwc, na.rm = T)) %>%
    ungroup() %>%
    group_by(date, type) %>%
    mutate(peak2 = mean_vwc == max(mean_vwc, na.rm = T))
  # stack to track_wetup
  track_wetup <- rbind(track_wetup, temp)
}

# subset to wetup events only
track_wetup2 <- subset(track_wetup, peak & grepl("wet", type) & mon %in% c(1:5, 11:12)) %>% #peak2
  group_by(nut_trt, ppt_trt) %>%
  mutate(difforder = cleanorder - lead(cleanorder)) %>%
  ungroup() %>%
  filter(difforder < -6 | is.na(difforder)) %>%
  group_by(nut_trt, ppt_trt) %>%
  mutate(final_event = 1:length(wetevent))

subset(testsm2, cleanorder %in% 8700:8800) %>% # it seems like if the change is at least 0.01 it will pick up the wetup event but not diurnal fluctuation
  ggplot(aes(cleanorder, mean_vwc)) +
  #geom_point(aes(color = ppt_mm), size = 2)+
  geom_point(aes(col = log(ppt_mm+0.00001)), size = 2)+
  geom_point(data = subset(wetups, cleanorder %in% 8700:8800), pch = 8, col = "red") +
  geom_point(data = subset(track_wetup, grepl("wet", type) & cleanorder %in% 8700:8800 & peak), pch = 8, col = "red") +
  geom_line(alpha = 0.75)+ 
  #scale_color_viridis_c()+
  scale_color_distiller(name = "ln(ppt)", palette = "Blues", direction = 1) +
  #scale_x_datetime(date_labels = "%m-%d", expand = c(0,0)) +
 #scale_y_continuous(breaks = seq(-0.1,1, 0.02)) +
  labs(x =NULL, y = "VWC") +
  facet_wrap(nut_trt~ppt_trt) +
  theme_ctw() +
  theme(panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey"))



subset(testsm2, mon %in% 2) %>% # it seems like if the change is at least 0.01 it will pick up the wetup event but not diurnal fluctuation
  ggplot(aes(clean_datetime, mean_vwc)) +
  #geom_point(aes(color = ppt_mm), size = 2)+
  geom_point(aes(col = log(ppt_mm+0.00001)), size = 2)+
  geom_point(data = subset(track_wetup2, mon %in% 2), pch = 8, col = "red") +
  geom_line(alpha = 0.75)+ 
  #scale_color_viridis_c()+
  scale_color_distiller(name = "ln(ppt)", palette = "Blues", direction = 1) +
  scale_x_datetime(date_labels = "%d", expand = c(0,0)) +
  #scale_y_continuous(breaks = seq(-0.1,1, 0.02)) +
  labs(x = "February 2021", y = "VWC") +
  facet_wrap(nut_trt~ppt_trt) +
  theme_ctw()
#theme(panel.grid.major = element_line(color = "grey"),
#      panel.grid.minor = element_line(color = "grey")) # i am satisfied with this. not perfect, but captures more or less single wetup events


# tally wetups by trt
count_wetups <- group_by(track_wetup2, ppt_trt, nut_trt) %>%
  reframe(totevents = max(final_event)) %>%
  group_by(ppt_trt) %>%
  mutate(mean_events = mean(totevents),
         se_events = sd(totevents)/sqrt(3)) %>%
  ungroup()

# there are only 3 obs per ppt trt (one for each nut trt, so just report summary stats in count_wetups)
summary(lm(totevents ~ -1 + ppt_trt,data = count_wetups))
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#   ppt_trtControl   22.000      2.434   9.037 0.000103 ***
#   ppt_trtWet       51.667      2.434  21.224 7.13e-07 ***
#   ppt_trtDrought   12.000      2.434   4.930 0.002632 ** 

tidy(lm(totevents ~ ppt_trt,data = count_wetups), conf.int = T)
# A tibble: 3 × 7
#  term           estimate std.error statistic  p.value conf.low conf.high
# <chr>             <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
# 1 (Intercept)        22        2.43      9.04 0.000103     16.0     28.0 
# 2 ppt_trtWet         29.7      3.44      8.62 0.000134     21.2     38.1 
# 3 ppt_trtDrought    -10        3.44     -2.90 0.0272      -18.4     -1.58



# time between events
timebetwn_wet <- arrange(track_wetup2, fulltrt, clean_datetime, final_event) %>%
  group_by(nut_trt, ppt_trt) %>%
  mutate(timebtwn = clean_datetime - lag(clean_datetime), # puts on seconds scale
         lagorder = cleanorder - lag(cleanorder),
         # put on days scale
         lagorder = (lagorder*2)/24) %>% 
  ungroup() %>%
  group_by(ppt_trt, nut_trt) %>%
  mutate(non_na_obs = length(lagorder[!is.na(lagorder)]),
         pn_mean_timebtwn = mean(timebtwn, na.rm = T),
         pn_mean_dayslag = mean(lagorder, na.rm= T),
         pn_se_dayslag = sd(lagorder, na.rm = T)/sqrt(non_na_obs)) %>%
  ungroup() %>%
  group_by(ppt_trt) %>%
  mutate(ppt_non_na = length(lagorder[!is.na(lagorder)]),
         ppt_mean_timebtwn = mean(timebtwn, na.rm = T),
         ppt_mean_dayslag = mean(lagorder, na.rm= T),
         ppt_se_dayslag = sd(lagorder, na.rm = T)/sqrt(non_na_obs)) %>%
  dplyr::select(grep("_trt|_mean_|_se|non_na", names(.))) %>%
  distinct()




# -- climate stats ----
# ppt distribution
mon2021gs <- subset(monthdat, waterYear == 2021) %>%
  # subset further to months of study
  subset(wymon %in% 1:8) %>%
  mutate(ppt_pct = round(100*(mon_ppt/sum(mon_ppt)),2))


#  -- anpp -----
anpp <- read_csv(paste0(datpath,"ANPP/ANPP_CleanedData/Compost_ANPP_clean.csv"))
str(anpp) 
anpp$nut_trt <- factor(anpp$nut_trt, levels = c("N", "C", "F"), labels = c("Control", "Compost", "Fert"))
anpp$ppt_trt <- factor(anpp$ppt_trt, levels = c("XC", "D", "W"), labels = c("Control", "Drought", "Wet"))
anpp$block <- factor(anpp$block)

# sum total ANPP (ignore fxnl grp), scale to g per m2 (clipped at 0.25m^2)
# there are 2 sampling events per year (april, may) and biomass is separated into grams and forbs
anpp_tot_sum <- subset(anpp, select = c(yr, clip_event, block, plotid, nut_trt, ppt_trt, fxnl_grp, dry_wgt_g)) %>%
  spread(fxnl_grp, dry_wgt_g) %>% # year 2020 has anpp missing for clip event 2 in 1NW and 2CD
  group_by(yr, block, plotid, nut_trt, ppt_trt) %>%
  # scales: ANPP clipped at 0.25 x 0.25cm, once in april, once in may
  # try sum april and may anpp
  reframe(anppF_0.5 = sum(Forb, na.rm = T), # sum = 0.25 x 0.5m (scale by 8 to get 1m2)
          anppG_0.5 = sum(Grass, na.rm = T),
          # as holdover, double anpp in obs that were missing for 2020
          anppF_0.5 = ifelse(plotid %in% c("1NW", "2CD") & yr == 2020, anppF_0.5*2, anppF_0.5),
          anppG_0.5 = ifelse(plotid %in% c("1NW", "2CD") & yr == 2020, anppG_0.5*2, anppF_0.5),
          # scale to 1m2
          anppF_1m = anppF_0.5*8,
          anppG_1m = anppG_0.5*8,
          # sum for totanpp
          anpp_1m = anppF_1m + anppG_1m
          )
hist(anpp_tot_sum$anpp_1m) # not normal dist
hist(scale(anpp_tot_sum$anpp_1m))
hist(log(anpp_to_sumt$anpp_1m))
hist(sqrt(anpp_tot_sum$anpp_1m))

anpp_tot <- subset(anpp, select = c(yr, clip_event, block, plotid, nut_trt, ppt_trt, fxnl_grp, dry_wgt_g)) %>%
  spread(fxnl_grp, dry_wgt_g) %>% # year 2020 has anpp missing for clip event 2 in 1NW and 2CD
  mutate(anpp_1m = (Forb + Grass)*8) # scale 0.25clip by 8 for 1m2

# look at data
ggplot(subset(anpp, yr == 2021), aes(fxnl_grp, dry_wgt_g, pch = factor(clip_event), col = block)) +
  geom_jitter(size = 3, width = 0.25, height = 0) +
  facet_wrap(nut_trt ~ ppt_trt)

ggplot(anpp_tot, aes(block, anpp_1m, col = factor(yr), group = paste(yr, clip_event))) +
  geom_point(aes(shape = factor(clip_event)),size = 3, position = position_dodge(width = 0.25)) +
  facet_grid(nut_trt ~ ppt_trt) +
  scale_color_manual(values = c("skyblue3", "chocolate2", "sienna")) +
  scale_shape_manual(values = c(1,8)) +
  theme_ctw()

# means by clip event
ggplot(anpp_tot, aes(factor(yr), anpp_1m, col = factor(clip_event), group = paste(yr, clip_event))) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), pch = 21) +
  stat_summary(position = position_dodge(width = 0.5)) +
  facet_grid(nut_trt ~ ppt_trt) +
  scale_color_manual(name = "Harvest", values= c("grey50", "grey10"), labels = c("April", "May")) +
  #scale_shape_manual(values = c(1,8)) +
  labs(x = NULL, y = "Total ANPP (g per m^2)") +
  theme_ctw()

# set factor groups
mod_anpp2021 <- glmmTMB(sqrt(dry_wgt_g) ~ ppt_trt + nut_trt + fxnl_grp + (1|block) + (1|clip_event), data = subset(anpp, yr == 2021))
summary(mod_anpp2021)
DHARMa::plotResiduals(mod_anpp2021) # residuals are normal with sqrt trans but i don't know what interp is..
# compare across years
summary(glmmTMB(dry_wgt_g ~ ppt_trt + nut_trt + fxnl_grp + (1|block) + (1|clip_event), data = subset(anpp, yr == 2019)))
summary(glmmTMB(dry_wgt_g ~ ppt_trt + nut_trt + fxnl_grp + (1|block) + (1|clip_event), data = subset(anpp, yr == 2020)))
mod_anpp_all <- glmmTMB(dry_wgt_g ~ ppt_trt + nut_trt + factor(yr) + fxnl_grp + (1|block) + (1|clip_event), data = anpp, family = "tweedie")
summary(mod_anpp_all)
DHARMa::plotResiduals(mod_anpp_all)                 
confint(mod_anpp_all)
broom.mixed::tidy(mod_anpp_all)

# use 1m2 data
# try total ANPP (one point per treat per year)
norm_totanpp_1m2 <- glmmTMB(anpp_1m ~ ppt_trt + nut_trt + factor(yr) + (1|block), data = anpp_tot_sum)
plotResiduals(norm_totanpp_1m2) # residuals not good
plotConventionalResiduals(norm_totanpp_1m2)
performance::check_distribution(norm_totanpp_1m2)
# residuals: cauchy (75%), tweedie (16%)
# response dist: lognormal (44%), tweedie (16%)


# try by clip event
anpp_tot$clip_event <- factor(anpp_tot$clip_event)
norm_clipanpp_1m2 <- glmmTMB(anpp_1m ~ ppt_trt + nut_trt + factor(yr) + (1|block) + (1|clip_event), data = anpp_tot) # # ideally would nest by plot id but cannot accomodate that error
plotResiduals(norm_clipanpp_1m2) # totally fine.. use this
hist(anpp_tot$anpp_1m) # separating into forb v grass is what creates problem
summary(norm_clipanpp_1m2)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     115.725     13.562   8.533  < 2e-16 ***
# ppt_trtDrought  -29.279      6.729  -4.351 1.35e-05 ***
# ppt_trtWet       -1.285      6.729  -0.191  0.84854    
# nut_trtCompost   20.328      6.753   3.010  0.00261 ** 
# nut_trtFert      -5.357      6.729  -0.796  0.42596    
# factor(yr)2020   12.514      6.754   1.853  0.06390 .  
# factor(yr)2021   44.876      6.705   6.693 2.18e-11 ***
anpp_tot$yrfac <- factor(anpp_tot$yr)
# compare with interactive model
norm_clipanppx_1m2 <- glmmTMB(anpp_1m ~ ppt_trt + nut_trt + yrfac + (1|block) + (1|clip_event), data = anpp_tot) # convergence issues if include interaction, convergence issues when include clip_event as well but leaving in..
summary(norm_clipanppx_1m2)
# pull marginal means from additive model
anpp_margmean <- emmeans(norm_clipanpp_1m2, specs = ~ nut_trt + ppt_trt | yr)
anpp_margmean
contrast(anpp_margmean, "trt.vs.ctrl") # effect of nutrient is the same across all years. not realistic

anpp_margmean_df <- data.frame(anpp_margmean)#broom.mixed::tidy(anpp_margmean)
ggplot(anpp_margmean_df, aes(ppt_trt, emmean, col = nut_trt)) + #col = ppt_trt, group = ppt_trt)) +
  geom_point(position = position_dodge(width = 0.25), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, position = position_dodge(width = 0.25)) +
  facet_grid(~yr) +
  labs(y = "ANPP marginal mean (g per m2) with 95% CI") +
  scale_color_brewer(palette = "Set2") +
  theme_ctw() +
  theme(legend.position = "inside",
        legend.justification.inside = c(0.99, 0.01))


# anpp models by year for comparison
norm_clipanpp2019_1m2 <- glmmTMB(anpp_1m ~ ppt_trt + nut_trt + (1|block), data = subset(anpp_tot, yr == 2019)) # convergence issue if include clip_event
summary(norm_clipanpp2019_1m2)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)     207.743     28.821   7.208 5.68e-13 ***
#   ppt_trtDrought  -27.127     21.187  -1.280   0.2004    
#   ppt_trtWet       14.277     21.187   0.674   0.5004    
#   nut_trtCompost   58.947     21.187   2.782   0.0054 ** 
#   nut_trtFert      -6.167     21.187  -0.291   0.7710    

plot(emmeans(norm_clipanpp2019_1m2, ~ nut_trt | ppt_trt)) + theme_ctw() + ggtitle("Marginal means: ANPP 2019 model (additive effects only)")
emmeans(norm_clipanpp2019_1m2, trt.vs.ctrl ~ nut_trt + ppt_trt)

# $contrasts
# contrast                          estimate   SE df t.ratio p.value
# Compost Control - Control Control    58.95 21.2 65   2.782  0.0455 * this is the only one signif diff
# Fert Control - Control Control       -6.17 21.2 65  -0.291  0.9969
# Control Drought - Control Control   -27.13 21.2 65  -1.280  0.6739
# Compost Drought - Control Control    31.82 30.0 65   1.062  0.7989
# Fert Drought - Control Control      -33.29 30.0 65  -1.111  0.7728
# Control Wet - Control Control        14.28 21.2 65   0.674  0.9483
# Compost Wet - Control Control        73.22 30.0 65   2.444  0.1026 * marginal comp wet greater than control
# Fert Wet - Control Control            8.11 30.0 65   0.271  0.9976

# 
# P value adjustment: dunnettx method for 8 tests 

# 2021 only
norm_clipanpp2021_1m2 <- glmmTMB(anpp_1m ~ ppt_trt + nut_trt + (1|block), data = subset(anpp_tot, yr == 2021)) # (1|clip_event)
norm_clipanpp2021x_1m2 <- glmmTMB(anpp_1m ~ ppt_trt * nut_trt + (1|block), data = subset(anpp_tot, yr == 2021)) # (1|clip_event)
anova(norm_clipanpp2021_1m2, norm_clipanpp2021x_1m2) # additive model is better
# Models:
#   norm_clipanpp2021_1m2: anpp_1m ~ ppt_trt + nut_trt + (1 | block), zi=~0, disp=~1
# norm_clipanpp2021x_1m2: anpp_1m ~ ppt_trt * nut_trt + (1 | block), zi=~0, disp=~1
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# norm_clipanpp2021_1m2   7 846.00 861.94 -416.00   832.00                         
# norm_clipanpp2021x_1m2 11 851.83 876.88 -414.92   829.83 2.1713      4     0.7043
summary(norm_clipanpp2021_1m2)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)     331.889     22.030  15.065  < 2e-16 ***
#   ppt_trtDrought  -77.290     22.156  -3.488 0.000486 ***
#   ppt_trtWet        4.733     22.156   0.214 0.830832    
#   nut_trtCompost   23.713     22.156   1.070 0.284496    
#   nut_trtFert     -14.410     22.156  -0.650 0.515447    

confint(norm_clipanpp2021_1m2)
#                               2.5 %    97.5 %   Estimate
# (Intercept)                288.710757 375.06698 331.888871
# ppt_trtDrought            -120.715462 -33.86435 -77.289906
# ppt_trtWet                 -38.692208  48.15890   4.733348
# nut_trtCompost             -19.712264  67.13885  23.713292
# nut_trtFert                -57.835561  29.01555 -14.410005
# Std.Dev.(Intercept)|block    4.104954  74.28927  17.462932

plotResiduals(norm_clipanpp2021_1m2) # ok
plotConventionalResiduals(norm_clipanpp2021_1m2)

# get pairwise differences
plot(emmeans(norm_clipanpp2021_1m2, ~ nut_trt | ppt_trt)) + 
  theme_ctw() + ggtitle("ANPP 2021 model: ANPP ~ nut_trt + ppt_trt + (1|block)") +
  labs(x = "Estimated marginal mean ANPP (g per m^2) with 95% CI",
       y = "Soil amendment")
emm1 <- emmeans(norm_clipanpp2021_1m2, ~ nut_trt + ppt_trt)
emm1
tidy(contrast(emm1, "trt.vs.ctrl"), conf.int = T)
# contrast                          estimate   SE df lower.CL upper.CL
# Compost Control - Control Control    23.71 22.2 65    -37.1    84.54
# Fert Control - Control Control      -14.41 22.2 65    -75.2    46.41
# Control Drought - Control Control   -77.29 22.2 65   -138.1   -16.47
# Compost Drought - Control Control   -53.58 31.3 65   -139.6    32.44
# Fert Drought - Control Control      -91.70 31.3 65   -177.7    -5.68
# Control Wet - Control Control         4.73 22.2 65    -56.1    65.56
# Compost Wet - Control Control        28.45 31.3 65    -57.6   114.47
# Fert Wet - Control Control           -9.68 31.3 65    -95.7    76.34
# 
# Confidence level used: 0.95 
# Conf-level adjustment: dunnettx method for 8 estimates 

contrast(emm1, "trt.vs.ctrl", simple = "nut_trt") %>%
   confint() # no diff between
contrast(emm1, "trt.vs.ctrl", by = "ppt_trt") %>%
  confint()
  

# -- functions to grab model stats ----
# > fuss with table summary packages later

# function to compare model stats and confidence intervals
getmod_effects <- function(mod, modlab = "anpp 2021"){
  coefdf <- tidy(mod, conf.int = T)
  coefdf$call <- str_flatten(as.character(mod$call))
  coefdf$mod <- modlab
  coefdf <- cbind(coefdf, glance(mod)) # redundant info but so it's there
  return(coefdf)
}

# function to compile emmeans, contrasts, with confidence intervals
getemmeans <- function(emobj, conobj, modlab = "anpp 2021"){
  marg <- tidy(emobj, conf.int = T)
  marg$mod <- modlab
  emdf <- tidy(conobj, conf.int = T)
  emdf$mod <- modlab
  tmplist <- list(marg, emdf)
  names(tmplist) <- c("marginal_means", "comparisons") 
  return(tmplist)  
}


# -- stats tables -----
# vwc
vwc2021_modsum <- getmod_effects(daily_vwc_mod_norm, modlab = "daily mean vwc WY2021 Nov-May")
emvwc <- emmeans(daily_vwc_mod_norm, ~ nut_trt * ppt_trt, type = "response")

vwc2021_emmeans <- getemmeans(emvwc, contrast(emvwc, "trt.vs.ctrl"))
vwc2021_pairwise <- getemmeans(emvwc, pairs(emvwc))
# get vwc diffs within soil treatments
ppt_diffs_bynut <- tidy(contrast(emvwc, "pairwise", simple = "nut_trt"), conf.int = T) # P value adjustment: tukey method for comparing a family of 3 estimates 
ppt_diffs_bynut$mod <- "daily mean vwc WY2021 Nov-May"

nut_diffs_byppt <- tidy(contrast(emvwc, "pairwise", simple = "ppt_trt"), conf.int = T)
nut_diffs_byppt$mod <- "daily mean vwc WY2021 Nov-May"

write_csv(vwc2021_modsum,"Native-analysis/outputs/vwc2021_coeffs.csv")
write_csv(vwc2021_emmeans$marginal_means,"Native-analysis/outputs/vwc2021_margmeans.csv")
write_csv(vwc2021_emmeans$comparisons,"Native-analysis/outputs/vwc2021_contrasts.csv")
write_csv(vwc2021_pairwise$comparisons,"Native-analysis/outputs/vwc2021_pairwise.csv")
write_csv(ppt_diffs_bynut, "Native-analysis/outputs/vwc2021_within_nuttrt.csv")
write_csv(nut_diffs_byppt, "Native-analysis/outputs/vwc2021_within_ppttrt.csv")
write_csv(count_wetups, "Native-analysis/outputs/wetups_by_ppt.csv")
write_csv(timebetwn_wet, "Native-analysis/outputs/wetups_returnint.csv")


# anpp
anpp2021_modsum <- getmod_effects(norm_clipanpp2021_1m2)
anpp2021_emmeans <- getemmeans(emm1, contrast(emm1, "trt.vs.ctrl"))

write_csv(anpp2021_modsum,"Native-analysis/outputs/anpp2021_coeffs.csv")
write_csv(anpp2021_emmeans$marginal_means,"Native-analysis/outputs/anpp2021_margmeans.csv")
write_csv(anpp2021_emmeans$comparisons,"Native-analysis/outputs/anpp2021_contrasts.csv")


# save mods for convenience and datasets
save(norm_clipanpp2021_1m2, daily_vwc_mod_norm, all_cimis, sm_summary_adj_mod, 
     file = "Native-analysis/outputs/supplemental_mods.rdata")

