print(paste0("Running Creating Daily Terrestrial Forecasts at ", Sys.time()))

renv::restore()

library(tidyverse)
library(lubridate)
library(rjags)
library(tidybayes)
library(modelr)
library(aws.s3)
library(prov)
library(EFIstandards)
library(EML)
library(jsonlite)

set.seed(329)

team_list <- list(list(individualName = list(givenName = "Quinn", surName = "Thomas"), 
                       id = "https://orcid.org/0000-0003-1282-7825"),
                  list(individualName = list(givenName = "Alex",  surName ="Young")),
                  list(individualName = list(givenName = "George",  surName ="Burba")),
                  list(individualName = list(givenName = "Jamie",  surName ="Cleverly")),
                  list(individualName = list(givenName = "Ankur",  surName ="Desai")),
                  list(individualName = list(givenName = "Mike",  surName ="Dietze")),
                  list(individualName = list(givenName = "Andy",  surName ="Fox")),
                  list(individualName = list(givenName = "William",  surName ="Hammond")),
                  list(individualName = list(givenName = "Danica",  surName ="Lombardozzi"))
)

team_name <- "pers_null_daily"
forecast_project_id <- "efi_null"

download.file("https://data.ecoforecast.org/targets/terrestrial/terrestrial-daily-targets.csv.gz",
              "terrestrial-daily-targets.csv.gz")

terrestrial_targets <- read_csv("terrestrial-daily-targets.csv.gz", guess_max = 10000)

site_names <- c("BART","KONZ","OSBS","SRER")

RandomWalk = "
model{

  #### Priors
  x[1] ~ dnorm(x_ic,tau_obs)
  tau_add ~ dgamma(0.1,0.1)
  tau_obs ~ dgamma(0.1,0.1)

  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
    x_obs[t] ~ dnorm(x[t],tau_obs)
  }

  #### Data Model
  for(i in 1:nobs){
    y[i] ~ dnorm(x[y_wgaps_index[i]], tau_obs)
  }

}
"

# NEE Model

forecast_saved_nee <- NULL

for(s in 1:length(site_names)){

  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max_time
  # This is key here - I added 35 days on the end of the data for the forecast period
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "1 day"))
  
  site_data_var <- left_join(full_time, site_data_var)
  

  
  #Full time series with gaps
  y_wgaps <- site_data_var$nee
  time <- c(site_data_var$time)
  #Remove gaps
  y_nogaps <- y_wgaps[!is.na(y_wgaps)]
  #Indexes of full time series with gaps
  y_wgaps_index <- 1:length(y_wgaps)
  #keep indexes to reference the gappy time series
  y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
  
  init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
  
  data <- list(y = y_nogaps,
               y_wgaps_index = y_wgaps_index,
               nobs = length(y_wgaps_index),
               n = length(y_wgaps),
               x_ic = 0.0)
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(tau_add = 1/var(diff(y_nogaps)),
                      tau_obs = 10 * 1/var(diff(y_nogaps)),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i],
                      x = init_x)
  }
  
  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           inits = init,
                           n.chains = 3)
  
  jags.out   <- coda.samples(model = j.model,variable.names = c("tau_add","tau_obs"), n.iter = 10000)
  
  m   <- coda.samples(model = j.model,
                      variable.names = c("x","tau_add","tau_obs", "x_obs"),
                      n.iter = 10000,
                      thin = 5)
  
  model_output <- m %>%
    spread_draws(x_obs[day]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = full_time$time[day]) %>%
    ungroup() %>%
    select(time, x_obs, ensemble)
  
  obs <- tibble(time = full_time$time,
                obs = y_wgaps)
  
  model_output %>% 
    group_by(time) %>% 
    summarise(mean = mean(x_obs),
              upper = quantile(x_obs, 0.975),
              lower = quantile(x_obs, 0.025),.groups = "drop") %>% 
    ggplot(aes(x = time, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
    geom_point(data = obs, aes(x = time, y = obs), color = "red") +
    labs(x = "Date", y = "nee")

    forecast_saved_tmp <- model_output %>%
      filter(time > start_forecast) %>%
      rename(nee = x_obs) %>% 
      mutate(data_assimilation = 0,
             forecast = 1,
             obs_flag = 2,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = team_name)
    
    forecast_saved_nee <- rbind(forecast_saved_nee, forecast_saved_tmp)

}

# Latent heat model

forecast_saved_le <- NULL

for(s in 1:length(site_names)){
  
  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max_time
  # This is key here - I added 16 days on the end of the data for the forecast period
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "1 day"))
  
  site_data_var <- left_join(full_time, site_data_var)
  
  # NEE
  
  #Full time series with gaps
  y_wgaps <- site_data_var$le
  time <- c(site_data_var$time)
  #Remove gaps
  y_nogaps <- y_wgaps[!is.na(y_wgaps)]
  #Indexes of full time series with gaps
  y_wgaps_index <- 1:length(y_wgaps)
  #keep indexes to reference the gappy time series
  y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
  
  init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
  
  data <- list(y = y_nogaps,
               y_wgaps_index = y_wgaps_index,
               nobs = length(y_wgaps_index),
               n = length(y_wgaps),
               x_ic = 0.0)
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(tau_add = 1/var(diff(y_nogaps)),
                      tau_obs = 10 * 1/var(diff(y_nogaps)),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i],
                      x = init_x)
  }
  
  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           inits = init,
                           n.chains = 3)
  
  jags.out   <- coda.samples(model = j.model,variable.names = c("tau_add","tau_obs"), n.iter = 10000)
  
  m   <- coda.samples(model = j.model,
                      variable.names = c("x","tau_add","tau_obs", "x_obs"),
                      n.iter = 10000,
                      thin = 5)
  
  model_output <- m %>%
    spread_draws(x_obs[day]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = full_time$time[day]) %>%
    ungroup() %>%
    select(time, x_obs, ensemble)
  
  obs <- tibble(time = full_time$time,
                obs = y_wgaps)
  
  model_output %>% 
    group_by(time) %>% 
    summarise(mean = mean(x_obs),
              upper = quantile(x_obs, 0.975),
              lower = quantile(x_obs, 0.025),.groups = "drop") %>% 
    ggplot(aes(x = time, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
    geom_point(data = obs, aes(x = time, y = obs), color = "red") +
    labs(x = "Date", y = "nee")
  
    forecast_saved_tmp <- model_output %>%
      filter(time > start_forecast) %>%
      rename(le = x_obs) %>% 
      mutate(data_assimilation = 0,
             forecast = 1,
             obs_flag = 2,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = team_name)
    
    forecast_saved_le <- rbind(forecast_saved_le, forecast_saved_tmp)
}


forecast_saved <- cbind(forecast_saved_nee, forecast_saved_le$le) %>% 
  rename(le = `forecast_saved_le$le`) %>% 
  select(time, ensemble, siteID, obs_flag, nee, le, forecast, data_assimilation)


forecast_file_name_base <- paste0("terrestrial-",as_date(start_forecast),"-",team_name)
forecast_file <- paste0(forecast_file_name_base, ".csv.gz")
write_csv(forecast_saved, forecast_file)

# Generate metadata

curr_time <- with_tz(Sys.time(), "UTC")
#forecast_issue_time <- format(curr_time,format = "%Y-%m-%d %H:%M:%SZ", usetz = F)
forecast_issue_time <- as_date(curr_time)
forecast_iteration_id <- start_forecast
forecast_model_id <- team_name

## define variable names, units, etc
## in practice, this might be kept in a spreadsheet
attributes <- tibble::tribble(
  ~attributeName,     ~attributeDefinition,                          ~unit,                  ~formatString,  ~definition, ~numberType,
  "time",              "[dimension]{time}",                          "year",                 "YYYY-MM-DD",   NA,          "datetime",
  "ensemble",          "[dimension]{index of ensemble member}",      "dimensionless",         NA,            NA,          "integer",
  "siteID",             "[dimension]{neon site}",                     NA,                     NA,           "NEON site ID",  "character",
  "obs_flag",          "[flag]{observation error}",                  "dimensionless",         NA,           NA,           "integer",
  "nee",               "[variable]{net ecosystem exchange}",         "numberPerMeterSquared", NA,           NA,           "real",
  "le",                "[variable]{latent heat}",                    "numberPerMeterSquared", NA,           NA,          "real",
  "forecast",          "[flag]{whether represents forecast}",        "dimensionless",         NA,           NA,          "integer",
  "data_assimilation", "[flag]{whether time step assimilated data}", "dimensionless",         NA,           NA,          "integer"
) 

## note: EML uses a different unit standard than UDUNITS. For now use EML. EFI needs to provide a custom unitList.
attrList <- EML::set_attributes(attributes, 
                                col_classes = c("Date", "numeric", "character","numeric", 
                                                "numeric","numeric", "numeric","numeric"))

physical <- set_physical(forecast_file)

dataTable <- eml$dataTable(
  entityName = "forecast",  ## this is a standard name to allow us to distinguish this entity from 
  entityDescription = "Forecast of NEE and LE for four NEON sites",
  physical = physical,
  attributeList = attrList)

#meta <- neonstore::neon_index(ext="xml", product = "DP4.00200.001")
#all <- lapply(meta$path, emld::as_emld)
#geo <- lapply(all, function(x) x$dataset$coverage$geographicCoverage)
#sites_ids <- lapply(geo, function(x) x$id) %>% unlist() 

#first_name <- rep(NA, length(site_names))
#for(i in 1:length(site_names)){
#  first_name[i] <- min(which(sites_ids == site_names[i]))
#  geo[[first_name[i]]]$boundingCoordinates$boundingAltitudes$altitudeMinimum <- round(as.numeric(geo[[first_name[i]]]$boundingCoordinates$boundingAltitudes$altitudeMinimum), 4)
#  geo[[first_name[i]]]$boundingCoordinates$boundingAltitudes$altitudeMaximum <- round(as.numeric(geo[[first_name[i]]]$boundingCoordinates$boundingAltitudes$altitudeMaximum), 4)
#}
#geo[first_name] %>% toJSON() %>% fromJSON() %>% distinct() %>% write_json("meta/terrestrial_geo.json", auto_unbox=TRUE)

temporalCoverage <- list(rangeOfDates =
                           list(beginDate = list(calendarDate = min(forecast_saved$time)),
                                endDate = list(calendarDate = max(forecast_saved$time))))

geographicCoverage = jsonlite::read_json("meta/terrestrial_geo.json")

coverage <- list(geographicCoverage = geographicCoverage,
                 temporalCoverage = temporalCoverage)

dataset = eml$dataset(
  title = "Daily persistence null forecast for nee and lee",
  creator = team_list,
  contact = list(references=team_list[[1]]$id),
  pubDate = as_date(forecast_issue_time),
  intellectualRights = "https://creativecommons.org/licenses/by/4.0/",
  dataTable = dataTable,
  coverage = coverage
)

#Minimal metdata
additionalMetadata <- eml$additionalMetadata(
  metadata = list(
    forecast = list(
      ## Basic elements
      timestep = "1 day", ## should be udunits parsable; already in coverage -> temporalCoverage?
      forecast_horizon = "35 days",
      forecast_issue_time = forecast_issue_time,
      forecast_iteration_id = forecast_iteration_id,
      forecast_project_id = forecast_project_id,
      metadata_standard_version = "0.3",
      model_description = list(
        forecast_model_id = forecast_model_id,
        name = "state-space Bayesian null",
        type = "empirical",
        repository = "https://github.com/eco4cast/neon4cast-terrestrial"
      ),
      ## MODEL STRUCTURE & UNCERTAINTY CLASSES
      initial_conditions = list(
        # Possible values: absent, present, data_driven, propagates, assimilates
        status = "assimilates",
        complexity = 2,
        propagation = list(
          type = "ensemble",
          size = max(forecast_saved$ensemble)),
        assimilation = list(
          type = "refit",
          reference = "NA",
          complexity = 4)
      ),
      drivers = list(
        status = "absent"
      ),
      parameters = list(
        status = "assimilates",
        complexity = 2,
        propagation = list(
          type = "ensemble",
          size = max(forecast_saved$ensemble)),
        assimilation = list(
          type = "refit",
          reference = "NA",
          complexity = 4)
      ),
      random_effects = list(
        status = "absent"
      ),
      process_error = list(
        status = "assimilates",
        complexity = 2,
        propagation = list(
          type = "ensemble",
          size = max(forecast_saved$ensemble)),
        assimilation = list(
          type = "refit",
          reference = "NA",
          complexity = 4),
        covariance = FALSE
      ),
      obs_error = list(
        status = "present",
        complexity = 2
      )
    ) # forecast
  ) # metadata
) # eml$additionalMetadata


my_eml <- eml$eml(dataset = dataset,
                  additionalMetadata = additionalMetadata,
                  packageId = forecast_iteration_id , 
                  system = "datetime"  ## system used to generate packageId
)

## check base EML
EML::eml_validate(my_eml)

EFIstandards::forecast_validator(my_eml)

meta_data_filename <-  paste0(forecast_file_name_base,".xml")

write_eml(my_eml, meta_data_filename)

## Publish the forecast automatically. (EFI-only)

source("../neon4cast-shared-utilities/publish.R")
publish(code = "03_terrestrial_flux_daily_null.R",
        data_in = "terrestrial-daily-targets.csv.gz",
        data_out = forecast_file,
        meta = meta_data_filename,
        prefix = "terrestrial/",
        bucket = "forecasts")
