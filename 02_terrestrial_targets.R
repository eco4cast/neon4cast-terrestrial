print(paste0("Running Creating Terrestrial Targets at ", Sys.time()))

#renv::restore()
Sys.setenv("NEONSTORE_HOME" = "/efi_neon_challenge/neonstore")
pecan_flux_uncertainty <- "../pecan/modules/uncertainty/R/flux_uncertainty.R"

source(pecan_flux_uncertainty)

#remotes::install_github("cboettig/neonstore")


library(neonUtilities)
library(neonstore)
library(tidyverse)
library(lubridate)
library(contentid)

site_names <- c("BART","KONZ","OSBS","SRER")

# Terrestrial
#DP4.00200.001 & DP1.00094.001
neon_store(product = "DP4.00200.001") 
flux_data <- neon_table(table = "nsae-basic")

flux_data <- flux_data %>% 
  mutate(time = as_datetime(timeBgn))

co2_data <- flux_data %>% 
  filter(qfqm.fluxCo2.turb.qfFinl == 0 & data.fluxCo2.turb.flux > -50 & data.fluxCo2.turb.flux < 50 & data.fluxMome.turb.veloFric >= 0.2) %>% 
  select(time,data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux) %>% 
  rename(nee = data.fluxCo2.turb.flux) %>% 
  mutate(nee = ifelse(siteID == "OSBS" & year(time) < 2019, NA, nee),
         nee = ifelse(siteID == "SRER" & year(time) < 2019, NA, nee)) %>% 
  rename(le = data.fluxH2o.turb.flux) %>% 
  mutate(le = ifelse(siteID == "OSBS" & year(time) < 2019, NA, le),
                  le = ifelse(siteID == "SRER" & year(time) < 2019, NA, le))

#ggplot(co2_data, aes(x = time, y = nee)) +
#  geom_point() +
#  facet_wrap(~siteID)

#Filter by qfFinal flow storage

#co2_data <- flux_data %>% 
#  filter(qfqm.fluxCo2.nsae.qfFinl == 0 & abs(data.fluxCo2.nsae.flux) < 50) %>% 
#  select(time,data.fluxCo2.nsae.flux, siteID) %>% 
#  rename(nee = data.fluxCo2.nsae.flux) %>% 
#  mutate(nee = ifelse(siteID == "OSBS" & year(time) < 2019, NA, nee),
#         nee = ifelse(siteID == "SRER" & year(time) < 2019, NA, nee))

#ggplot(co2_data, aes(x = time, y = nee)) +
#  geom_point() +
#  facet_wrap(~siteID)



#le_data <-  flux_data %>% 
#  filter(qfqm.fluxH2o.turb.qfFinl == 0) %>% 
#  select(time,data.fluxH2o.turb.flux, siteID)%>% 
#  rename(le = data.fluxH2o.turb.flux) %>% 
#  mutate(le = ifelse(siteID == "OSBS" & year(time) < 2019, NA, le),
#         le = ifelse(siteID == "SRER" & year(time) < 2019, NA, le))

#ggplot(le_data, aes(x = time, y = le)) +
#  geom_point() +
#  facet_wrap(~siteID)

earliest <- min(as_datetime(c(co2_data$time)), na.rm = TRUE)
latest <- max(as_datetime(c(co2_data$time)), na.rm = TRUE)


full_time <- seq(min(c(co2_data$time), na.rm = TRUE), 
                 max(c(co2_data$time), na.rm = TRUE), 
                 by = "30 min")

full_time <- tibble(time = rep(full_time, 4),
                    siteID = c(rep("BART", length(full_time)),
                               rep("KONZ", length(full_time)),
                               rep("OSBS", length(full_time)),
                               rep("SRER", length(full_time))))


flux_target_30m <- left_join(full_time, co2_data, by = c("time", "siteID"))

#ggplot(flux_target_30m, aes(x = time, y = nee)) +
#  geom_point() +
#  facet_wrap(~siteID)


valid_dates_nee <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  filter(!is.na(nee)) %>% # & !is.na(le)) %>% 
  group_by(date, siteID) %>% 
  summarise(count = n()) %>% 
  filter(count >= 44)

valid_dates_le <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  filter(!is.na(le)) %>% # & !is.na(le)) %>% 
  group_by(date, siteID) %>% 
  summarise(count = n()) %>% 
  filter(count >= 44)

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

#ggplot(flux_target_daily, aes(x = time, y = nee)) + 
#  geom_point() +
#  facet_wrap(~siteID)


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
nee_sd_slopeP[which(nee_sd_slopeN < 0)] <- 0
le_sd_slopeN[which(le_sd_slopeN > 0)] <- 0
le_sd_slopeP[which(le_sd_slopeN < 0)] <- 0

site_uncertainty <- tibble(siteID = site_names,
                           nee_sd_intercept = nee_intercept,
                           nee_sd_slopeP = nee_sd_slopeP,
                           nee_sd_slopeN = nee_sd_slopeN,
                           le_sd_intercept = le_intercept,
                           le_sd_slopeP = le_sd_slopeP,
                           le_sd_slopeN = le_sd_slopeN)


flux_target_30m <- left_join(flux_target_30m, site_uncertainty, by = "siteID")

########

neon_store(table = "SWS_30_minute", n = 50) 
neon_store(table = "sensor_positions", n = 50) 
d2 <- neon_read(table = "sensor_positions") 
sm30 <- neon_table(table = "SWS_30_minute")
sensor_positions <- neon_table(table = "sensor_positions")

sensor_positions <- sensor_positions %>%  
  filter(str_detect(referenceName, "SOIL")) %>% 
  mutate(horizontalPosition = str_sub(HOR.VER, 1, 3),
         verticalPosition = str_sub(HOR.VER, 5, 7),
         verticalPosition = as.numeric(verticalPosition),
         siteID = str_sub(file, 10, 13),
         horizontalPosition = as.numeric(horizontalPosition)) %>% 
  rename(sensorDepths = zOffset) %>% 
  filter(siteID %in% c("KONZ", "BART", "OSBS", "SRER")) %>% 
  select(sensorDepths, horizontalPosition, verticalPosition, siteID)
# 



sm30 <- sm30 %>% 
  mutate(horizontalPosition = as.numeric(horizontalPosition),
         verticalPosition = as.numeric(verticalPosition))
sm3_combined <- left_join(sm30, sensor_positions, by = c("siteID", "verticalPosition", "horizontalPosition"))
   
 
sm3_combined <- sm3_combined %>% 
   select(startDateTime, endDateTime, VSWCMean, siteID, horizontalPosition, verticalPosition, VSWCFinalQF, sensorDepths, VSWCExpUncert) %>% 
   mutate(VSWCMean = as.numeric(VSWCMean)) %>% 
   filter(VSWCFinalQF == 0,
          VSWCMean > 0) %>% 
  mutate(sensorDepths = -sensorDepths)
 
sm3_combined <- sm3_combined %>% 
   filter(horizontalPosition == 3 & sensorDepths < 0.30) %>% 
  mutate(depth = NA,
    depth = ifelse(sensorDepths <= 0.07, 0.05, depth),
          depth = ifelse(sensorDepths > 0.07 & sensorDepths < 0.20, 0.15, depth),
          depth = ifelse(sensorDepths > 0.20 & sensorDepths < 0.30, 0.25, depth)) %>% 
  filter(depth == 0.15) %>% 
  rename(time = startDateTime)


earliest <- min(as_datetime(c(flux_target_30m$time)), na.rm = TRUE)
latest <- max(as_datetime(c(sm3_combined$time)), na.rm = TRUE)


full_time <- seq(earliest, 
                 latest, 
                 by = "30 min")

full_time <- tibble(time = rep(full_time, 4),
                    siteID = c(rep("BART", length(full_time)),
                               rep("KONZ", length(full_time)),
                               rep("OSBS", length(full_time)),
                               rep("SRER", length(full_time))))

sm3_combined <- left_join(full_time, sm3_combined, by = c("time", "siteID"))



sm3_combined %>% 
  ggplot(aes(x = time, y = VSWCMean, color = factor(depth))) + 
  geom_point() +
  facet_wrap(~siteID, scale = "free")

sm30_target <- sm3_combined %>% 
  mutate(time = lubridate::as_datetime(time)) %>% 
  rename(vswc = VSWCMean,
         vswc_sd = VSWCExpUncert) %>% 
  select(time, siteID, vswc, vswc_sd) 
  

sm_daily_target <- sm3_combined %>% 
  select(time, siteID, VSWCMean, VSWCExpUncert) %>% 
  mutate(time = lubridate::as_date(time)) %>% 
  group_by(time, siteID) %>% 
    summarise(vswc = mean(VSWCMean, na.rm = TRUE),
              count = sum(!is.na(VSWCMean)),
              vswc_sd = mean(VSWCExpUncert)/sqrt(count)) %>% 
    dplyr::mutate(vswc = ifelse(count > 23, vswc, NA)) %>% 
    dplyr::select(time, siteID, vswc, vswc_sd)
    
terrestrial_target_30m <- full_join(flux_target_30m, sm30_target)

ggplot(terrestrial_target_30m, aes(x = time, y = vswc)) +
  geom_point() +
  facet_wrap(~siteID)

terrestrial_target_daily <- full_join(flux_target_daily, sm_daily_target)

#ggplot(terrestrial_target_daily, aes(x = time, y = nee)) +
#  geom_point() +
#  facet_wrap(~siteID)

write_csv(terrestrial_target_30m, "terrestrial-30min-targets.csv.gz")
write_csv(terrestrial_target_daily, "terrestrial-daily-targets.csv.gz")

## Publish the targets to EFI.  Assumes aws.s3 env vars are configured.
source("../neon4cast-shared-utilities/publish.R")
publish(code = c("02_terrestrial_targets.R"),
        data_out = c("terrestrial-30min-targets.csv.gz","terrestrial-daily-targets.csv.gz"),
        prefix = "terrestrial/",
        bucket = "targets")

