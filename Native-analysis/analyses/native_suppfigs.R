# supplemental figures
# ctw

# notes:
# just making rain summary for SLM presentation 2022 Apr 5 for now
# more later...


# -- SETUP ----
library(tidyverse)
library(lubridate)
library(SPEI)
library(cimir)
# modify default settings
options(stringsAsFactors = F)
theme_set(theme_test())
na_vals <- c("" , " ","NA", NA)

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
  
cimis_ts_tmean <- ts(all_cimis$infill_values[all_cimis$Item == "DayAirTmpAvg"], start = c(2009,01), freq = 365)
cimis_PET <- thornthwaite(ts(cimis_monthly$mon_tmean, start = c(2009,01), end = c(2023, 7), freq = 12), lat = 39.25170524042971)
cimis_BAL <- ts(cimis_monthly$mon_ppt, start = c(2009, 01), end = c(2023, 7), freq=12) - cimis_PET
cimis_SPEI.1 <- spei(cimis_BAL, 1)
cimis_SPEI.3 <- spei(cimis_BAL, 3)
plot(cimis_SPEI.1) # missing data are making full spei calculations funky

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
# look at trend
summary(lm(trend ~ date, data = subset(temptrends, temp == "tmin"))) # tmin and tmean are both about +1deg per decade like tmax
with(subset(temptrends, !is.na(trend)), range(date))
temptrends_datespres <- subset(temptrends, !is.na(trend))
max(temptrends_datespres$date) - min(temptrends_datespres$date)


#build plot
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

monthdat <- cbind(cimis_monthly, 
                  spei3 = as.numeric(cimis_SPEI.3$fitted),
                  ppt_trend = as.numeric(ppt_trend_monthly$trend)
                  ) %>%
  mutate(date = as.Date(paste(yr, mon, 01, sep = "-"))) %>%
  left_join(distinct(all_cimis[c("Date", "doy", "dowy")]), by = c("date"= "Date")) %>%
  distinct()

speiplot <- monthdat[c("date", "spei3")] %>%
  mutate(nondrought = ifelse(spei3>=0, spei3, 0),
         drought = ifelse(spei3<0, spei3, 0)) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), lty = 1) +
  # add vlines at experiment start and end
  geom_vline(aes(xintercept = as.Date("2018-10-01")), lty = 2, col = "black") +
  geom_vline(aes(xintercept = as.Date("2021-06-01")), lty = 2, col = "black") +
  geom_area(aes(date, nondrought), fill = "royalblue1", col = "grey30", alpha = 0.9) +
  geom_area(aes(date, drought), fill = "salmon3", col = "grey30", alpha = 0.9) +
  labs(y = "3-month SPEI", x = "Date")+
  theme(plot.margin = margin(0,0,0,0, "pt"))

pptplot <- ggplot(monthdat) +
  geom_point(aes(date, mon_ppt), pch = 21, col = "slateblue2", fill = "deepskyblue2", alpha = 0.7) +
  geom_line(aes(date, ppt_trend), col = "steelblue3", lwd = 1.5, alpha = 0.8) +
  labs(y = "Monthly precip (mm)", x = NULL) +
  theme(axis.text.x = element_blank(),
        plot.margin = margin(0,0,0,0, "pt"))

plot_grid(pptplot, speiplot, nrow = 2, align = "vh")


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
figpath <- "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/esa2023_presfigs/"
ggsave(paste0(figpath,"sfrec_3monspei_20092023.pdf"), units = "in", width = 7, height = 3)

# no experiment timeline markings
ggsave(paste0(figpath,"sfrec_3monspei_20092023_nomarks.pdf"), units = "in", width = 7, height = 3)
# native exp only
ggsave(paste0(figpath,"sfrec_3monspei_20092023_natexo_marks.pdf"), units = "in", width = 7, height = 3)


# -- SOIL MOISTURE v. PREPPED PPT ---
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
  geom_line() +
  #geom_line(aes(cleanorder, upper, lty = nut_trt)) +
  #geom_line(aes(cleanorder, lower, lty = nut_trt)) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  facet_wrap(~ppt_trt, nrow = 3) +
  labs(x = NULL, y = "Mean VWC") +
  scale_color_brewer(name = NULL, palette = "Accent") +
  theme(legend.position = c(0.95, 0.97),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.key.height = unit(1, "pt"))

vwc_ppt_plot <-  ggplot(subset(compost_sm_ppt, waterYear == 2021 & !mon %in% 7:9), 
       aes(as.POSIXct(clean_datetime, tz = "UTC"), ppt_mm)) +
  geom_col(col = "dodgerblue3") +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Native experiment timeframe (Oct 2020 - Jun 2021)", y = "Precip (mm)")

plot_grid(meanvwc_plot, vwc_ppt_plot, nrow = 2, align = "hv", rel_heights = c(1,0.5))
