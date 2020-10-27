print(paste0("Running Creating Terrestrial Targets at ", Sys.time()))

renv::restore()
Sys.setenv("NEONSTORE_HOME" = "/efi_neon_challenge/neonstore")

library(neonUtilities)
library(neonstore)
library(tidyverse)
library(lubridate)
library(contentid)

run_full_workflow <- TRUE
generate_null_daily <- TRUE

# Terrestrial
#DP4.00200.001 & DP1.00094.001
sites <- c("BART", "KONZ", "SRER", "OSBS")
start_date <- NA

print("Downloading: DP4.00200.001")
new_data1 <- neonstore::neon_download(product = "DP4.00200.001", site = sites, type = "basic", start_date = start_date, .token = Sys.getenv("NEON_TOKEN"))
print("Downloading: DP1.00094.001")
new_data2 <- neonstore::neon_download(product = "DP1.00094.001", site = sites, type = "basic", start_date = start_date, .token = Sys.getenv("NEON_TOKEN"))

if(!is.null(new_data1) | !is.null(new_data2) | run_full_workflow){
  
  source("02_terrestrial_targets.R")
  
  print(paste0("Completed Target at ", Sys.time()))
  
  if(generate_null_daily){
    
    print(paste0("Running daily Null at ", Sys.time()))
    source("03_terrestrial_flux_daily_null.R")
    
  }
}

########

# neon_store(table = "SWS_30_minute", n = 500) 
# d2 <- neon_read(table = "sensor_positions") 
# sm30 <- neon_table(table = "SWS_30_minute")
# sensor_positions <- neon_table(table = "sensor_positions")
# 
# sensor_positions <- sensor_positions %>% 
#   mutate(horizontalPosition = str_sub(sensor_positions$HOR.VER, 1, 3),
#          verticalPosition = str_sub(HOR.VER, 5, 7),
#          siteID = str_sub(file, 10, 13)) %>% 
#   rename(sensorDepths = zOffset) %>% 
#   filter(siteID %in% c("KONZ", "BART", "OSBS", "SRER")) %>% 
#   select(sensorDepths, horizontalPosition, verticalPosition, siteID)
# 
# sm30 <- left_join(sm30, sensor_positions, by = c("siteID", "verticalPosition", "horizontalPosition"))
#   
# 
# sm30 <- sm30 %>% 
#   select(startDateTime, endDateTime, VSWCMean, siteID, horizontalPosition, verticalPosition, VSWCFinalQF) %>% 
#   mutate(VSWCMean = as.numeric(VSWCMean)) %>% 
#   filter(VSWCFinalQF == 0,
#          VSWCMean > 0)
# 
# 
# ## sd has Sensor Depths
# sensor_depth <-read.csv("inputs/SWC_depths.csv")
# 
# # multiply negative depth values by -1 for easier math
# sensor_depth$sensorDepths <- sensor_depth$sensorDepths * -1
# # only use the first sensor per site (typically the closest to each tower)
# sensor_depth <- sensor_depth[sensor_depth$plot == 1, ]
# # site-plot-meas. Used for matching up sensor depths
# sensor_depth$spd <-
#   paste(
#     sensor_depth$site,
#     sensor_depth$plot,
#     sensor_depth$measurementLevel
#   )
# 
# ### The 4 chunks beow calculate sensor widths to calculating their proportion of the depth profile (i.e. weight)
# ## it is repetitive, and a function could be used to replace it. But it is functional
# # BART
# bart_sd <- sensor_depth[sensor_depth$site == "BART", ]
# bart_sd$width <- 0 # made a new column, a numerical placeholderto be filled with sensor widths
# bart_sd[1, 9] <- ((bart_sd[2, 5] - bart_sd[1, 5]) / 2) + bart_sd[1, 5]
# bart_sd[2, 9] <- ((bart_sd[2, 5] - bart_sd[1, 5]) / 2) + ((bart_sd[3, 5] - bart_sd[2, 5]) / 2)
# bart_sd[3, 9] <- ((bart_sd[3, 5] - bart_sd[2, 5]) / 2) + ((bart_sd[4, 5] - bart_sd[3, 5]) / 2)
# bart_sd[4, 9] <- ((bart_sd[4, 5] - bart_sd[3, 5]) / 2) + ((bart_sd[5, 5] - bart_sd[4, 5]) / 2)
# bart_sd[5, 9] <- ((bart_sd[5, 5] - bart_sd[4, 5]) / 2) + ((bart_sd[6, 5] - bart_sd[5, 5]) / 2)
# bart_sd[6, 9] <- ((bart_sd[6, 5] - bart_sd[5, 5]) / 2) + ((bart_sd[7, 5] - bart_sd[6, 5]) / 2)
# bart_sd[7, 9] <- ((bart_sd[7, 5] - bart_sd[6, 5]) / 2) + ((bart_sd[8, 5] - bart_sd[7, 5]) / 2)
# bart_sd[8, 9] <- ((bart_sd[8, 5] - bart_sd[7, 5]) / 2)
# bart_sd$weight <- bart_sd$width / bart_sd[8, 5] # divide by the lowest sensor depth
# 
# # OSBS
# osbs_sd <- sensor_depth[sensor_depth$site == "OSBS", ]
# osbs_sd$width <- 0 # made a new column, a numerical placeholderto be filled with sensor widths
# osbs_sd[1, 9] <- ((osbs_sd[2, 5] - osbs_sd[1, 5]) / 2) + osbs_sd[1, 5]
# osbs_sd[2, 9] <- ((osbs_sd[2, 5] - osbs_sd[1, 5]) / 2) + ((osbs_sd[3, 5] - osbs_sd[2, 5]) / 2)
# osbs_sd[3, 9] <- ((osbs_sd[3, 5] - osbs_sd[2, 5]) / 2) + ((osbs_sd[4, 5] - osbs_sd[3, 5]) / 2)
# osbs_sd[4, 9] <- ((osbs_sd[4, 5] - osbs_sd[3, 5]) / 2) + ((osbs_sd[5, 5] - osbs_sd[4, 5]) / 2)
# osbs_sd[5, 9] <- ((osbs_sd[5, 5] - osbs_sd[4, 5]) / 2) + ((osbs_sd[6, 5] - osbs_sd[5, 5]) / 2)
# osbs_sd[6, 9] <- ((osbs_sd[6, 5] - osbs_sd[5, 5]) / 2) + ((osbs_sd[7, 5] - osbs_sd[6, 5]) / 2)
# osbs_sd[7, 9] <- ((osbs_sd[7, 5] - osbs_sd[6, 5]) / 2) + ((osbs_sd[8, 5] - osbs_sd[7, 5]) / 2)
# osbs_sd[8, 9] <- ((osbs_sd[8, 5] - osbs_sd[7, 5]) / 2)
# osbs_sd$weight <- osbs_sd$width / osbs_sd[8, 5] # divide by the lowest sensor depth
# 
# # SRER
# srer_sd <- sensor_depth[sensor_depth$site == "SRER", ]
# srer_sd$width <- 0 # made a new column, a numerical placeholderto be filled with sensor widths
# srer_sd[1, 9] <- ((srer_sd[2, 5] - srer_sd[1, 5]) / 2) + srer_sd[1, 5]
# srer_sd[2, 9] <- ((srer_sd[2, 5] - srer_sd[1, 5]) / 2) + ((srer_sd[3, 5] - srer_sd[2, 5]) / 2)
# srer_sd[3, 9] <- ((srer_sd[3, 5] - srer_sd[2, 5]) / 2) + ((srer_sd[4, 5] - srer_sd[3, 5]) / 2)
# srer_sd[4, 9] <- ((srer_sd[4, 5] - srer_sd[3, 5]) / 2) + ((srer_sd[5, 5] - srer_sd[4, 5]) / 2)
# srer_sd[5, 9] <- ((srer_sd[5, 5] - srer_sd[4, 5]) / 2) + ((srer_sd[6, 5] - srer_sd[5, 5]) / 2)
# srer_sd[6, 9] <- ((srer_sd[6, 5] - srer_sd[5, 5]) / 2) + ((srer_sd[7, 5] - srer_sd[6, 5]) / 2)
# srer_sd[7, 9] <- ((srer_sd[7, 5] - srer_sd[6, 5]) / 2) + ((srer_sd[8, 5] - srer_sd[7, 5]) / 2)
# srer_sd[8, 9] <- ((srer_sd[8, 5] - srer_sd[7, 5]) / 2)
# srer_sd$weight <- srer_sd$width / srer_sd[8, 5] # divide by the lowest sensor depth
# 
# # KONZ
# konz_sd <- sensor_depth[sensor_depth$site == "KONZ", ]
# konz_sd$width <- 0 # made a new column, a numerical placeholderto be filled with sensor widths
# konz_sd[1, 9] <- ((konz_sd[2, 5] - konz_sd[1, 5]) / 2) + konz_sd[1, 5]
# konz_sd[2, 9] <- ((konz_sd[2, 5] - konz_sd[1, 5]) / 2) + ((konz_sd[3, 5] - konz_sd[2, 5]) / 2)
# konz_sd[3, 9] <- ((konz_sd[3, 5] - konz_sd[2, 5]) / 2) + ((konz_sd[4, 5] - konz_sd[3, 5]) / 2)
# konz_sd[4, 9] <- ((konz_sd[4, 5] - konz_sd[3, 5]) / 2) + ((konz_sd[5, 5] - konz_sd[4, 5]) / 2)
# konz_sd[5, 9] <- ((konz_sd[5, 5] - konz_sd[4, 5]) / 2) + ((konz_sd[6, 5] - konz_sd[5, 5]) / 2)
# konz_sd[6, 9] <- ((konz_sd[6, 5] - konz_sd[5, 5]) / 2) + ((konz_sd[7, 5] - konz_sd[6, 5]) / 2)
# konz_sd[7, 9] <- ((konz_sd[7, 5] - konz_sd[6, 5]) / 2) + ((konz_sd[8, 5] - konz_sd[7, 5]) / 2)
# konz_sd[8, 9] <- ((konz_sd[8, 5] - konz_sd[7, 5]) / 2)
# konz_sd$weight <- konz_sd$width / konz_sd[8, 5] # divide by the lowest sensor depth
# 
# weighted_sensors <- rbind(bart_sd, osbs_sd, srer_sd, konz_sd)
# weighted_sensors
# 
# # Use the half hourly data
# # turn into numeric for later matching with sensor depth measurements
# sm30$horizontalPosition <- as.numeric(sm30$horizontalPosition)
# sm30$verticalPosition <- as.numeric(sm30$verticalPosition)
# # subtract 500 to get rid of the leading '50_' in the vertical position codes.
# sm30$verticalPosition <- sm30$verticalPosition - 500
# # add day of year
# sm30$day <- lubridate::yday(sm30$startDateTime)
# 
# sensor_depth <- sensor_depth %>% 
#   rename(horizontalPosition = plot,
#          verticalPosition = measurementLevel,
#          siteID = site) %>% 
#   select(siteID, verticalPosition, horizontalPosition, sensorDepths)
# 
# sm30 <- left_join(sm30, sensor_depth, by = c("siteID", "verticalPosition", "horizontalPosition"))
# 
# 
# View(sm30 %>% 
#        filter(startDateTime > as_datetime("2019-08-20 23:30:00"), 
#                 startDateTime < as_datetime("2019-08-30 23:30:00"), 
#                 horizontalPosition == "001") %>% 
#        select(startDateTime, siteID, sensorDepths, VSWCMean) %>% 
#        pivot_wider(names_from = startDateTime, values_from = VSWCMean))
# 
# sm30 %>% filter(siteID == "OSBS" & horizontalPosition == "001") %>% 
#   ggplot(aes(x = startDateTime, y = VSWCMean)) +
#   geom_point()+
#   facet_wrap(~factor(sensorDepths))
# 
# 
# d <- sm30 %>% group_by(siteID, sensorDepths, startDateTime) %>% summarise(moisture = mean(VSWCMean, na.rm = TRUE))
# 
# ggplot(d, aes(x = startDateTime, y = moisture, color = factor(sensorDepths))) +
#   geom_point() +
#   facet_grid(~siteID)
# 
# 
# d <- sm30 %>% filter(horizontalPosition == "001" & startDateTime == as_datetime("2020-08-05 21:00:00") & siteID == "KONZ") %>% 
#   select(startDateTime, VSWCMean, siteID, sensorDepths) %>% 
#   drop_na()
# 
# d <- sm30 %>% filter(horizontalPosition == "001") %>% 
#   select(startDateTime, VSWCMean, siteID, sensorDepths) %>% 
#   drop_na()
# 
# date_site <- d %>% group_by(startDateTime,siteID) %>% summarise()
# 
# 
# tmp_var <- NULL
# for(i in 1:nrow(date_site)){
#   d2 <- d %>% filter(startDateTime == date_site$startDateTime[i] & siteID == date_site$siteID[i])
#   tmp_var <- c(tmp_var, mean(d2$VSWCMean * diff(c(0, d2$sensorDepths))/mean(diff(c(0, d2$sensorDepths)))))
# }
# 
# vswc_30m <- cbind(date_site, tmp_var) 
# 
# names(vswc_30m) <- c("time", "siteID", "vswc")
# 
# 
# vswc_30m <- sensor1 %>% 
#   group_by(siteID, startDateTime) %>%
#   summarise(moisture = sum(weighted_moisture, na.rm = FALSE), .groups = "drop")
# 
# vswc_30m %>% 
#   ggplot(aes(x = time, y = vswc)) +
#   geom_line() +
#   facet_wrap(~siteID)
# 
# vswc_daily <- sensor1 %>% 
#   group_by(siteID, startDateTime) %>%
#   summarise(moisture = sum(weighted_moisture, na.rm = FALSE), .groups = "drop") %>% 
#   mutate(date = as_date(startDateTime)) %>% 
#   group_by(siteID, date) %>%
#   summarise(moisture = mean(moisture, na.rm = TRUE)) %>% 
#   ggplot(aes(x = date, y = moisture)) +
#   geom_point() +
#   facet_wrap(~siteID)


