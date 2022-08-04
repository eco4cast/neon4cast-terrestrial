print(paste0("Running Creating Terrestrial Targets at ", Sys.time()))

Sys.setenv("NEONSTORE_HOME" = "/home/rstudio/data/neonstore")
#Sys.setenv("NEONSTORE_DB" = "/home/rstudio/data/neonstore")
#Sys.setenv("NEONSTORE_DB")
pecan_flux_uncertainty <- "../pecan/modules/uncertainty/R/flux_uncertainty.R"
readRenviron("~/.Renviron") # compatible with littler

non_store_dir <- "/home/rstudio/data/neon_flux_data"
use_5day_data <- TRUE

source(pecan_flux_uncertainty)

library(neonUtilities)
library(neonstore)
library(tidyverse)
library(lubridate)
library(contentid)

sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-terrestrial/master/Terrestrial_NEON_Field_Site_Metadata_20210928.csv")

site_names <- sites$field_site_id

#print("Downloading: DP4.00200.001")
#neonstore::neon_download(product = "DP4.00200.001", site = site_names, type = "basic")
#neon_store(product = "DP4.00200.001") 
#print("Downloading: DP1.00094.001")
#x <- neonstore::neon_download(product = "DP1.00094.001", site = site_names, type = "basic")
#neon_store(table = "SWS_30_minute", n = 50) 

# Terrestrial

#Get the published files on the portal

#DP4.00200.001 & DP1.00094.001
#neon_store(product = "DP4.00200.001") 
flux_data <- neon_table(table = "nsae-basic", site = site_names) %>% 
  mutate(timeBgn = as_datetime(timeBgn),
         timeEnd = as_datetime(timeEnd))

#Get the current unpublished flux data (5-day latency)

if(use_5day_data){
  files <- readr::read_csv("https://storage.googleapis.com/neon-sae-files/ods/sae_files_unpublished/sae_file_url_unpublished.csv")
 #Convert old S3 links to GCS
  files$url <- base::gsub("https://s3.data.neonscience.org", "https://storage.googleapis.com", files$url)
  files <- files %>% 
    filter(site %in% site_names) %>% 
    mutate(file_name = basename(url)) |> 
    filter(date > max(lubridate::as_date(flux_data$timeBgn)))
  
  if(!dir.exists(file.path(non_store_dir,"current_month"))){
    dir.create(file.path(non_store_dir,"current_month"), recursive = TRUE)
  }
  
  for(i in 1:nrow(files)){
    destfile <- file.path(non_store_dir,"current_month",files$file_name[i])
    if(!(file.exists(destfile) | file.exists(tools::file_path_sans_ext(destfile)))){
      download.file(files$url[i], destfile = destfile)
      R.utils::gunzip(destfile)
    }
  }
  
  fn <- list.files(file.path(non_store_dir,"current_month"), full.names = TRUE)
  
  #remove files that are no longer in the unpublished s3 bucket
  
  for(i in 1:length(fn)){
    if(!(paste0(basename(fn[i]),".gz") %in% files$file_name)){
      unlink(fn[i])
    }
  }
  
  message(paste0("reading in ", length(fn), " non-NEON portal files"))
  
  flux_data_curr <- purrr::map_dfr(1:length(fn), function(i, fn){
    message(paste0("reading file ",fn[i]))
    neonstore::neon_read(files = fn[i])
  },
  fn = fn)
  
  #flux_data_curr <- neonstore::neon_read(files = fn)
  
  #remove any files unpublished data that has been published
  
  #Combined published and unpublished
  
  flux_data <- bind_rows(flux_data, flux_data_curr)
}

flux_data <- flux_data %>% 
  mutate(time = as_datetime(timeBgn))

co2_data <- flux_data %>% 
  filter(qfqm.fluxCo2.turb.qfFinl == 0 & data.fluxCo2.turb.flux > -50 & data.fluxCo2.turb.flux < 50) %>% 
  select(time,data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux) %>% 
  rename(nee = data.fluxCo2.turb.flux) %>% 
  mutate(nee = ifelse(siteID == "OSBS" & year(time) < 2019, NA, nee),
         nee = ifelse(siteID == "SRER" & year(time) < 2019, NA, nee)) %>% 
  rename(le = data.fluxH2o.turb.flux) %>% 
  mutate(le = ifelse(siteID == "OSBS" & year(time) < 2019, NA, le),
         le = ifelse(siteID == "SRER" & year(time) < 2019, NA, le))

co2_data %>% 
ggplot(aes(x = time, y = nee)) +
  geom_point() +
  facet_wrap(~siteID)

earliest <- min(as_datetime(c(co2_data$time)), na.rm = TRUE)
latest <- max(as_datetime(c(co2_data$time)), na.rm = TRUE)


full_time_vector <- seq(min(c(co2_data$time), na.rm = TRUE), 
                 max(c(co2_data$time), na.rm = TRUE), 
                 by = "30 min")

full_time <- NULL
for(i in 1:length(site_names)){
  df <- tibble(time = full_time_vector,
               siteID = rep(site_names[i], length(full_time_vector)))
  full_time <- bind_rows(full_time, df)
  
}

flux_target_30m <- left_join(full_time, co2_data, by = c("time", "siteID"))

#flux_target_30m %>% 
#  mutate(pass = ifelse(is.na(nee), 0, 1)) %>% 
#  group_by(siteID) %>% 
#  summarize(all = n(),
#            pass = sum(pass)) %>% 
#  mutate(prop = pass/all)

valid_dates_nee <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  filter(!is.na(nee)) %>% # & !is.na(le)) %>% 
  group_by(date, siteID) %>% 
  summarise(count = n()) %>% 
  filter(count >= 24)

#valid_dates_nee %>% 
#  group_by(siteID) %>% 
#  count() %>% 
#  mutate(all = length(unique(as_date(flux_target_30m$time))),
#         prop = n/all)

valid_dates_le <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  filter(!is.na(le)) %>% # & !is.na(le)) %>% 
  group_by(date, siteID) %>% 
  summarise(count = n()) %>% 
  filter(count >= 24)

flux_target_daily <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  group_by(date, siteID) %>% 
  summarize(nee = mean(nee, na.rm = TRUE),
            le = mean(le, na.rm = TRUE)) %>% 
  mutate(nee = ifelse(date %in% valid_dates_nee$date, nee, NA),
         le = ifelse(date %in% valid_dates_le$date, le, NA),
         nee = ifelse(is.nan(nee), NA, nee),
         le = ifelse(is.nan(le),NA, le)) %>% 
  rename(time = date) %>% 
  mutate(nee = (nee * 12 / 1000000) * (60 * 60 * 24))

flux_target_daily %>% 
  filter(year(time) > 2021) %>% 
  ggplot(aes(x = time, y = nee)) + 
  geom_point() +
  facet_wrap(~siteID)


nee_intercept <- rep(NA, length(site_names))
nee_sd_slopeP <- rep(NA, length(site_names))
nee_sd_slopeN <- rep(NA, length(site_names)) 
le_intercept <- rep(NA, length(site_names))
le_sd_slopeP <- rep(NA, length(site_names))
le_sd_slopeN <- rep(NA, length(site_names)) 

for(s in 1:length(site_names)){
  
  temp <- flux_target_30m %>%
    filter(siteID == site_names[s])
  
  unc <- flux.uncertainty(measurement = temp$nee, 
                          QC = rep(0, length(temp$nee)),
                          bin.num = 30)
  
  nee_intercept[s] <- unc$intercept
  nee_sd_slopeP[s] <- unc$slopeP
  nee_sd_slopeN[s] <- unc$slopeN
  
  unc <- flux.uncertainty(measurement = temp$le, 
                          QC = rep(0, length(temp$le)),
                          bin.num = 30)
  
  le_intercept[s] <- unc$intercept
  le_sd_slopeP[s] <- unc$slopeP
  le_sd_slopeN[s] <- unc$slopeN
  
  
}

nee_sd_slopeN[which(nee_sd_slopeN > 0)] <- 0
nee_sd_slopeP[which(nee_sd_slopeP < 0)] <- 0
le_sd_slopeN[which(le_sd_slopeN > 0)] <- 0
le_sd_slopeP[which(le_sd_slopeP < 0)] <- 0

site_uncertainty <- tibble(siteID = site_names,
                           nee_sd_intercept = nee_intercept,
                           nee_sd_slopeP = nee_sd_slopeP,
                           nee_sd_slopeN = nee_sd_slopeN,
                           le_sd_intercept = le_intercept,
                           le_sd_slopeP = le_sd_slopeP,
                           le_sd_slopeN = le_sd_slopeN)


flux_target_30m <- left_join(flux_target_30m, site_uncertainty, by = "siteID")

#flux_target_30m %>% duplicates(index = time, key = siteID)

# ########
# 
# #neon_store(table = "SWS_30_minute", n = 50) 
# #neon_store(table = "sensor_positions", n = 50) 
# #d2 <- neon_read(table = "sensor_positions") 
# sm30 <- neon_table(table = "SWS_30_minute")
# sensor_positions <- neon_read(table = "sensor_positions", product = "DP1.00094.001", keep_filename = TRUE)
# 
# sensor_positions <- sensor_positions %>%  
#   filter(str_detect(referenceName, "SOIL")) %>% 
#   mutate(horizontalPosition = str_sub(HOR.VER, 1, 3),
#          verticalPosition = str_sub(HOR.VER, 5, 7),
#          verticalPosition = as.numeric(verticalPosition),
#          siteID = str_sub(file, 10, 13),
#          horizontalPosition = as.numeric(horizontalPosition),
#          zOffset = as.numeric(zOffset)) %>% 
#   rename(sensorDepths = zOffset) %>% 
#   filter(siteID %in% c("KONZ", "BART", "OSBS", "SRER")) %>% 
#   select(sensorDepths, horizontalPosition, verticalPosition, siteID)
# # 
# 
# 
# sm30 <- sm30 %>% 
#   mutate(horizontalPosition = as.numeric(horizontalPosition),
#          verticalPosition = as.numeric(verticalPosition))
# sm3_combined <- left_join(sm30, sensor_positions, by = c("siteID", "verticalPosition", "horizontalPosition"))
# 
# 
# 
# 
# sm3_combined <- sm3_combined %>% 
#   select(startDateTime, endDateTime, VSWCMean, siteID, horizontalPosition, verticalPosition, VSWCFinalQF, sensorDepths, VSWCExpUncert) %>% 
#   mutate(VSWCMean = as.numeric(VSWCMean)) %>% 
#   filter(VSWCFinalQF == 0,
#          VSWCMean > 0) %>% 
#   mutate(sensorDepths = -sensorDepths)
# 
# #sm3_combined %>% duplicates(index = time, key = siteID)
# 
# sm30 %>% 
#   ggplot(aes(x = startDateTime, y = VSWCMean, color = factor(verticalPosition))) +
#   geom_point() +
#   facet_grid(siteID~horizontalPosition)
# 
# sm3_combined <- sm3_combined %>% 
#   filter(horizontalPosition == 3 & sensorDepths < 0.30) %>% 
#   mutate(depth = NA,
#          depth = ifelse(sensorDepths <= 0.07, 0.05, depth),
#          depth = ifelse(sensorDepths > 0.07 & sensorDepths < 0.20, 0.15, depth),
#          depth = ifelse(sensorDepths > 0.20 & sensorDepths < 0.30, 0.25, depth)) %>% 
#   filter(depth == 0.15) %>% 
#   rename(time = startDateTime)
# 
# 
# #sm3_combined %>% duplicates(index = time, key = siteID)
# 
# earliest <- min(as_datetime(c(flux_target_30m$time)), na.rm = TRUE)
# latest <- max(as_datetime(c(sm3_combined$time)), na.rm = TRUE)
# 
# 
# full_time <- seq(earliest, 
#                  latest, 
#                  by = "30 min")
# 
# full_time <- tibble(time = rep(full_time, 4),
#                     siteID = c(rep("BART", length(full_time)),
#                                rep("KONZ", length(full_time)),
#                                rep("OSBS", length(full_time)),
#                                rep("SRER", length(full_time))))
# 
# sm3_combined <- left_join(full_time, sm3_combined, by = c("time", "siteID"))
# 
# 
# 
# sm3_combined %>% 
#   ggplot(aes(x = time, y = VSWCMean, color = factor(depth))) + 
#   geom_point() +
#   facet_wrap(~siteID, scale = "free")
# 
# sm30_target <- sm3_combined %>% 
#   mutate(time = lubridate::as_datetime(time)) %>% 
#   rename(vswc = VSWCMean,
#          vswc_sd = VSWCExpUncert) %>% 
#   select(time, siteID, vswc, vswc_sd) 
# 
# 
# sm_daily_target <- sm3_combined %>% 
#   select(time, siteID, VSWCMean, VSWCExpUncert) %>% 
#   mutate(time = lubridate::as_date(time)) %>% 
#   group_by(time, siteID) %>% 
#   summarise(vswc = mean(VSWCMean, na.rm = TRUE),
#             count = sum(!is.na(VSWCMean)),
#             vswc_sd = mean(VSWCExpUncert)/sqrt(count)) %>% 
#   dplyr::mutate(vswc = ifelse(count > 23, vswc, NA)) %>% 
#   dplyr::select(time, siteID, vswc, vswc_sd)
# 
# terrestrial_target_30m <- full_join(flux_target_30m, sm30_target)
# 
# ggplot(terrestrial_target_30m, aes(x = time, y = vswc)) +
#   geom_point() +
#   facet_wrap(~siteID)
# 
# terrestrial_target_daily <- full_join(flux_target_daily, sm_daily_target)
# 
# ggplot(terrestrial_target_daily, aes(x = time, y = vswc)) +
#   geom_point() +
#   facet_wrap(~siteID)

write_csv(flux_target_30m, "terrestrial_30min-targets.csv.gz")
write_csv(flux_target_daily, "terrestrial_daily-targets.csv.gz")

## Publish the targets to EFI.  Assumes aws.s3 env vars are configured.
source("../challenge-ci/R/publish.R")
publish(code = c("02_terrestrial_targets.R"),
        data_out = c("terrestrial_30min-targets.csv.gz"),
        prefix = "terrestrial_30min/",
        bucket = "neon4cast-targets",
        registries = "https://hash-archive.carlboettiger.info")

source("../challenge-ci/R/publish.R")
publish(code = c("02_terrestrial_targets.R"),
        data_out = c("terrestrial_daily-targets.csv.gz"),
        prefix = "terrestrial_daily/",
        bucket = "neon4cast-targets",
        registries = "https://hash-archive.carlboettiger.info")

unlink("terrestrial_30min-targets.csv.gz")
unlink("terrestrial_daily-targets.csv.gz")

