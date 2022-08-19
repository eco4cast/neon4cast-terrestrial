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

sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(terrestrial == 1)

site_names <- sites$field_site_id

#Sys.unsetenv("AWS_DEFAULT_REGION")
#Sys.unsetenv("AWS_S3_ENDPOINT")
#Sys.setenv("AWS_EC2_METADATA_DISABLED"="TRUE")
#neon <- arrow::s3_bucket("neon4cast-targets/neon",
#                         endpoint_override = "data.ecoforecast.org",
#                         anonymous = TRUE)


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
  files <- readr::read_csv("https://storage.googleapis.com/neon-sae-files/ods/sae_files_unpublished/sae_file_url_unpublished.csv", show_col_types = FALSE)
 #Convert old S3 links to GCS
  files$url <- base::gsub("https://s3.data.neonscience.org", "https://storage.googleapis.com", files$url)
  files <- files %>% 
    filter(site %in% site_names) %>% 
    mutate(file_name = basename(url)) |> 
    filter(date > max(lubridate::as_date(flux_data$timeBgn)))
  
  fs::dir_create(file.path(non_store_dir,"current_month"), recurse = TRUE)
  fs::dir_create(file.path(non_store_dir,"current_month_parquet"), recurse = TRUE)
  
  for(i in 1:nrow(files)){
    destfile <- file.path(non_store_dir,"current_month",files$file_name[i])
    parquet_file <- file.path(non_store_dir,"current_month_parquet",paste0(tools::file_path_sans_ext(files$file_name[i]),".parquet"))
    if(!(file.exists(parquet_file))){
      download.file(files$url[i], destfile = destfile)
      R.utils::gunzip(destfile)
    }
  }
  
  fn_parquet <- list.files(file.path(non_store_dir,"current_month_parquet"))
  
  #remove files that are no longer in the unpublished s3 bucket because they are now in NEON portal
  
  for(i in 1:length(fn_parquet)){
    if(!(paste0(tools::file_path_sans_ext(basename(fn_parquet[i])),".gz") %in% files$file_name)){
      unlink(fn_parquet[i])
    }
  }
  
  fn <- list.files(file.path(non_store_dir,"current_month"), full.names = TRUE)
  #fn_parquet <- file.path(non_store_dir,"current_month_parquet", paste0(basename(fn[i]),".parquet"))
  
  message(paste0("reading in ", length(fn), " non-NEON portal files"))
  
  #future::plan("future::multisession", workers = 4)
  
  #new_files <- fn[which(!(basename(fn) %in% tools::file_path_sans_ext(basename(fn_parquet))))]
  
  purrr::walk(1:length(fn), function(i, fn){
    message(paste0(i, " of ", length(fn), " reading file ",fn[i]))
    df <- neonstore::neon_read(files = fn[i])
    arrow::write_parquet(x = df, file.path(non_store_dir,"current_month_parquet", paste0(basename(fn[i]),".parquet")))
    unlink(fn[i])
  },
  fn = fn)
  
  s3 <- arrow::SubTreeFileSystem$create(file.path(non_store_dir,"current_month_parquet"))
  
  flux_data_curr <- arrow::open_dataset(s3) |> 
    collect()

  
  #Combined published and unpublished
  
  flux_data <- bind_rows(flux_data, flux_data_curr)
}

flux_data |> 
  group_by(siteID) |> 
  summarize(min = min(timeBgn),
            max = max(timeBgn)) |> 
  arrange(min) |> 
  print(n = 50) 

flux_data <- flux_data %>% 
  mutate(time = as_datetime(timeBgn))

co2_data <- flux_data %>% 
  filter(qfqm.fluxCo2.turb.qfFinl == 0 & data.fluxCo2.turb.flux > -50 & data.fluxCo2.turb.flux < 50) %>% 
  select(time,data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux) %>% 
  rename(nee = data.fluxCo2.turb.flux,
         le = data.fluxH2o.turb.flux,
         site_id = siteID) |> 
  mutate(nee = ifelse(site_id == "OSBS" & year(time) < 2019, NA, nee),
         nee = ifelse(site_id == "SRER" & year(time) < 2019, NA, nee),
         nee = ifelse(site_id == "BARR" & year(time) < 2019, NA, nee),
         le = ifelse(site_id == "OSBS" & year(time) < 2019, NA, le),
         le = ifelse(site_id == "SRER" & year(time) < 2019, NA, le),
         le = ifelse(site_id == "BARR" & year(time) < 2019, NA, le)) |> 
  pivot_longer(-c("time","site_id"), names_to = "variable", values_to = "observed")

co2_data %>% 
filter(variable == "nee") |> 
ggplot(aes(x = time, y = observed)) +
  geom_point() +
  facet_wrap(~site_id)

earliest <- min(as_datetime(c(co2_data$time)), na.rm = TRUE)
latest <- max(as_datetime(c(co2_data$time)), na.rm = TRUE)


full_time_vector <- seq(min(c(co2_data$time), na.rm = TRUE), 
                 max(c(co2_data$time), na.rm = TRUE), 
                 by = "30 min")

full_time <- NULL
for(i in 1:length(site_names)){
  df_nee <- tibble(time = full_time_vector,
               site_id = rep(site_names[i], length(full_time_vector)),
               variable = "nee")
  df_le <- tibble(time = full_time_vector,
               site_id = rep(site_names[i], length(full_time_vector)),
               variable = "le")
  full_time <- bind_rows(full_time, df_nee, df_le)
  
}

flux_target_30m <- left_join(full_time, co2_data, by = c("time", "site_id", "variable"))

valid_dates <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  filter(!is.na(observed)) %>%
  group_by(date, site_id, variable) %>% 
  summarise(count = n(), .groups = "drop")

flux_target_daily <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  group_by(date, site_id, variable) %>% 
  summarize(observed = mean(observed, na.rm = TRUE)) |> 
  left_join(valid_dates, by = c("date","site_id", "variable")) |> 
  mutate(observed = ifelse(count > 24, observed, NA),
         observed = ifelse(is.nan(observed), NA, observed)) %>% 
  rename(time = date) %>% 
  select(-count) |> 
  mutate(observed = ifelse(variable == "nee", (observed * 12 / 1000000) * (60 * 60 * 24), observed))

flux_target_daily %>% 
  filter(year(time) > 2021) %>% 
  ggplot(aes(x = time, y = observed)) + 
  geom_point() +
  facet_grid(variable~site_id, scale = "free")

# Adding observational uncertainity to the 30 minute fluxes
# 
# nee_intercept <- rep(NA, length(site_names))
# nee_sd_slopeP <- rep(NA, length(site_names))
# nee_sd_slopeN <- rep(NA, length(site_names)) 
# le_intercept <- rep(NA, length(site_names))
# le_sd_slopeP <- rep(NA, length(site_names))
# le_sd_slopeN <- rep(NA, length(site_names)) 
# 
# for(s in 1:length(site_names)){
#   
#   temp <- flux_target_30m %>%
#     filter(site_id == site_names[s],
#            variable == "nee")
#   
#   unc <- flux.uncertainty(measurement = temp$observed, 
#                           QC = rep(0, length(temp$observed)),
#                           bin.num = 30)
#   
#   nee_intercept[s] <- unc$intercept
#   nee_sd_slopeP[s] <- unc$slopeP
#   nee_sd_slopeN[s] <- unc$slopeN
#   
#   temp <- flux_target_30m %>%
#     filter(site_id == site_names[s],
#            variable == "le")
#   
#   unc <- flux.uncertainty(measurement = temp$observed, 
#                           QC = rep(0, length(temp$observed)),
#                           bin.num = 30)
#   
#   le_intercept[s] <- unc$intercept
#   le_sd_slopeP[s] <- unc$slopeP
#   le_sd_slopeN[s] <- unc$slopeN
#   
#   
# }
# 
# nee_sd_slopeN[which(nee_sd_slopeN > 0)] <- 0
# nee_sd_slopeP[which(nee_sd_slopeP < 0)] <- 0
# le_sd_slopeN[which(le_sd_slopeN > 0)] <- 0
# le_sd_slopeP[which(le_sd_slopeP < 0)] <- 0
# 
# nee_site_uncertainty <- tibble(site_id = site_names,
#                            variable = "nee",
#                            sd_intercept = nee_intercept,
#                            sd_slopeP = nee_sd_slopeP,
#                            sd_slopeN = nee_sd_slopeN)
# 
# le_site_uncertainty <- tibble(site_id = site_names,
#                                variable = "le",
#                                sd_intercept = le_intercept,
#                                sd_slopeP = le_sd_slopeP,
#                                sd_slopeN = le_sd_slopeN)
# 
# site_uncertainty <- bind_rows(nee_site_uncertainty, le_site_uncertainty)
# 
# 
# flux_target_30m <- left_join(flux_target_30m, site_uncertainty, by = c("site_id", "variable"))

flux_target_30m <- flux_target_30m |> 
  select(time, site_id, variable, observed)

flux_target_daily <- flux_target_daily |> 
  select(time, site_id, variable, observed)

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

