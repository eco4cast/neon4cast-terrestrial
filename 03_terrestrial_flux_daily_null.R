#'# Ecological Forecasting Initiative Null Model 

#'## Set-up

print(paste0("Running Creating Daily Terrestrial Forecasts at ", Sys.time()))

#'Load renv.lock file that includes the versions of all the packages used
#'You can generate using the command renv::snapshot()
renv::restore()

#' Required packages.  
#' EFIstandards is at remotes::install_github("eco4cast/EFIstandards")
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

#' set the random number for reproducible MCMC runs
set.seed(329)

#'Generate plot to visualized forecast
generate_plots <- TRUE
#'Is the forecast run on the Ecological Forecasting Initiative Server?
#'Setting to TRUE published the forecast on the server.
efi_server <- TRUE

#' List of team members. Used in the generation of the metadata
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

#'Team name code
team_name <- "pers_null_daily"

#'Download target file from the server
download.file("https://data.ecoforecast.org/targets/terrestrial/terrestrial-daily-targets.csv.gz",
              "terrestrial-daily-targets.csv.gz")

#'Read in target file.  The guess_max is specified because there could be a lot of
#'NA values at the beginning of the file
terrestrial_targets <- read_csv("terrestrial-daily-targets.csv.gz", guess_max = 10000)

#'Focal sites
site_names <- c("BART","KONZ","OSBS","SRER")

#'Generic random walk state-space model is JAGS format.  We use this model for 
#'both the NEE and LE null forecasts
RandomWalk = "
model{

  # Priors
  x[1] ~ dnorm(x_ic,tau_obs)
  tau_add ~ dgamma(0.1,0.1)
  tau_obs ~ dgamma(0.1,0.1)

  # Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
    x_obs[t] ~ dnorm(x[t],tau_obs)
  }

  # Data Model
  for(i in 1:nobs){
    y[i] ~ dnorm(x[y_wgaps_index[i]], tau_obs)
  }

}
"

#'## NEE Model

#'Create variable for combined forecasts across sites
forecast_saved_nee <- NULL
nee_figures <- list()
#+ message = FALSE
#' Loop through sites
for(s in 1:length(site_names)){
  
  # Select site
  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  # Find the last day in the observed data and add one day for the start of the 
  # forecast
  start_forecast <- max(site_data_var$time) + days(1)
  
  # This is key here - I added 35 days on the end of the data for the forecast period
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "1 day"))
  
  # Join the full time with the site_data_var so there aren't gaps in the time column
  site_data_var <- left_join(full_time, site_data_var)
  
  #observed NEE: Full time series with gaps
  y_wgaps <- site_data_var$nee
  time <- c(site_data_var$time)
  #observed NEE: time series without gaps
  y_nogaps <- y_wgaps[!is.na(y_wgaps)]
  #Index: time series with gaps
  y_wgaps_index <- 1:length(y_wgaps)
  #Index: the index of the non-NA values in time series with gaps
  y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
  
  #Generate starting initial conditions for latent states
  init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
  
  #Create a list of the data for use in JAGS.  Include vector lengths (nobs, n)
  data <- list(y = y_nogaps,
               y_wgaps_index = y_wgaps_index,
               nobs = length(y_wgaps_index),
               n = length(y_wgaps),
               x_ic = 0.0)
  
  #Initialize parameters 
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
  
  #Initialize JAGS model
  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           inits = init,
                           n.chains = 3)
  
  #Run JAGS model as the burn-in
  jags.out   <- coda.samples(model = j.model,variable.names = c("tau_add","tau_obs"), n.iter = 10000)
  
  #Run JAGS model again and sample from the posteriors
  m   <- coda.samples(model = j.model,
                      variable.names = c("x","tau_add","tau_obs", "x_obs"),
                      n.iter = 10000,
                      thin = 5)
  
  #Use TidyBayes package to clean up the JAGS output
  model_output <- m %>%
    spread_draws(x_obs[day]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = full_time$time[day]) %>%
    ungroup() %>%
    select(time, x_obs, ensemble)
  
  if(generate_plots){
    #Pull in the observed data for plotting
    obs <- tibble(time = full_time$time,
                  obs = y_wgaps)
    
    
    #Post past and future
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
    
    ggsave(paste0("nee_daily_",site_names[s],"_figure.pdf"), device = "pdf")
  }
  
  #Filter only the forecasted dates and add columns for required variable
  forecast_saved_tmp <- model_output %>%
    filter(time > start_forecast) %>%
    rename(nee = x_obs) %>% 
    mutate(data_assimilation = 0,
           forecast = 1,
           obs_flag = 2,
           siteID = site_names[s]) %>%
    mutate(forecast_iteration_id = start_forecast) %>%
    mutate(forecast_project_id = team_name)
  
  # Combined with the previous sites
  forecast_saved_nee <- rbind(forecast_saved_nee, forecast_saved_tmp)
  
}

#'## Latent heat model
#' 
#' See notes from the NEE section above
#+ message = FALSE
forecast_saved_le <- NULL
le_figures <- list()
for(s in 1:length(site_names)){
  
  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max_time
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "1 day"))
  
  site_data_var <- left_join(full_time, site_data_var)
  
  y_wgaps <- site_data_var$le
  time <- c(site_data_var$time)
  
  y_nogaps <- y_wgaps[!is.na(y_wgaps)]
  
  y_wgaps_index <- 1:length(y_wgaps)
  
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
                      tau_obs = 1/var(diff(y_nogaps)),
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
  
  if(generate_plots){
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
      labs(x = "Date", y = "le")
    
    ggsave(paste0("le_daily_",site_names[s],"_figure.pdf"), device = "pdf")
  }
  
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

#'## Soil moisture
#' 
#' See notes from the NEE section above
#+ message = FALSE
forecast_saved_soil_moisture <- NULL
soil_moisture_figures <- list()
for(s in 1:length(site_names)){
  
  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max_time
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "1 day"))
  
  site_data_var <- left_join(full_time, site_data_var)
  
  y_wgaps <- site_data_var$vswc
  time <- c(site_data_var$time)
  
  y_nogaps <- y_wgaps[!is.na(y_wgaps)]
  
  y_wgaps_index <- 1:length(y_wgaps)
  
  y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
  
  init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
  
  data <- list(y = y_nogaps,
               y_wgaps_index = y_wgaps_index,
               nobs = length(y_wgaps_index),
               n = length(y_wgaps),
               x_ic = 0.3)
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(tau_add = 1/var(diff(y_nogaps)),
                      tau_obs = 1/var(diff(y_nogaps)),
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
  
  if(generate_plots){
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
      labs(x = "Date", y = "vswc", title = site_names[s])
    
    ggsave(paste0("sm_daily_",site_names[s],"_figure.pdf"), device = "pdf")
  }
  
  forecast_saved_tmp <- model_output %>%
    filter(time > start_forecast) %>%
    rename(vswc = x_obs) %>% 
    mutate(data_assimilation = 0,
           forecast = 1,
           obs_flag = 2,
           siteID = site_names[s]) %>%
    mutate(forecast_iteration_id = start_forecast) %>%
    mutate(forecast_project_id = team_name)
  
  forecast_saved_soil_moisture <- rbind(forecast_saved_soil_moisture, forecast_saved_tmp)
}

#'Combined the NEE and LE forecasts together and re-order column
forecast_saved <- cbind(forecast_saved_nee, forecast_saved_le$le, forecast_saved_soil_moisture$vswc) %>% 
  rename(le = `forecast_saved_le$le`,
         vswc = `forecast_saved_soil_moisture$vswc`) %>% 
  select(time, ensemble, siteID, obs_flag, nee, le, vswc, forecast, data_assimilation)

#'Save file as CSV in the
#'[theme_name]-[year]-[month]-[date]-[team_name].csv
forecast_file_name_base <- paste0("terrestrial-",as_date(start_forecast),"-",team_name)
forecast_file <- paste0(forecast_file_name_base, ".csv.gz")
write_csv(forecast_saved, forecast_file)

#'#Generate metadata

#'Get system time for setting the issue time of the forecast
curr_time <- with_tz(Sys.time(), "UTC")
#forecast_issue_time <- format(curr_time,format = "%Y-%m-%d %H:%M:%SZ", usetz = F)
forecast_issue_time <- as_date(curr_time)
forecast_iteration_id <- start_forecast

#' The team name is the `forecast_model_id`
forecast_model_id <- team_name

#'Build attribute table. Other models won't likely change this for the challenge
#'note: need to fix the units for nee and le because the units do not pass EML unit checks
attributes <- tibble::tribble(
  ~attributeName,     ~attributeDefinition,                          ~unit,                  ~formatString,  ~definition, ~numberType,
  "time",              "[dimension]{time}",                          "year",                 "YYYY-MM-DD",   NA,          "datetime",
  "ensemble",          "[dimension]{index of ensemble member}",      "dimensionless",         NA,            NA,          "integer",
  "siteID",             "[dimension]{neon site}",                     NA,                     NA,           "NEON site ID",  "character",
  "obs_flag",          "[flag]{observation error}",                  "dimensionless",         NA,           NA,           "integer",
  "nee",               "[variable]{net ecosystem exchange}",         "dimensionless",         NA,           NA,           "real",
  "vswc",              "[variable]{volumetric soil water content}",  "dimensionless",         NA,           NA,           "real",
  "le",                "[variable]{latent heat}",                    "dimensionless",         NA,           NA,          "real",
  "forecast",          "[flag]{whether represents forecast}",        "dimensionless",         NA,           NA,          "integer",
  "data_assimilation", "[flag]{whether time step assimilated data}", "dimensionless",         NA,           NA,          "integer"
) 

#' use `EML` package to build the attribute list
attrList <- EML::set_attributes(attributes, 
                                col_classes = c("Date", "numeric", "character","numeric","numeric", 
                                                "numeric","numeric", "numeric","numeric"))
#' use `EML` package to build the physical list
physical <- EML::set_physical(forecast_file)

#' use `EML` package to dataTable
dataTable <- eml$dataTable(
  entityName = "forecast",  ## this is a standard name to allow us to distinguish this entity from 
  entityDescription = "Forecast of NEE and LE for four NEON sites",
  physical = physical,
  attributeList = attrList)

#'This code is for generating the geographicCoverage from NEON supplied EML
#'It uses the `neonstore` package. The for-loop extracts only one set of 
#'geographicCoverage for each site.  We have already extracted and saved as a
#'JSON file.
#'
#+ eval=FALSE
meta <- neonstore::neon_index(ext="xml", product = "DP4.00200.001")
all <- lapply(meta$path, emld::as_emld)
geo <- lapply(all, function(x) x$dataset$coverage$geographicCoverage)
sites_ids <- lapply(geo, function(x) x$id) %>% unlist()

first_name <- rep(NA, length(site_names))
for(i in 1:length(site_names)){
  first_name[i] <- min(which(sites_ids == site_names[i]))
  geo[[first_name[i]]]$boundingCoordinates$boundingAltitudes$altitudeMinimum <- round(as.numeric(geo[[first_name[i]]]$boundingCoordinates$boundingAltitudes$altitudeMinimum), 4)
  geo[[first_name[i]]]$boundingCoordinates$boundingAltitudes$altitudeMaximum <- round(as.numeric(geo[[first_name[i]]]$boundingCoordinates$boundingAltitudes$altitudeMaximum), 4)
}
geo[first_name] %>% toJSON() %>% fromJSON() %>% distinct() %>% write_json("meta/terrestrial_geo.json", auto_unbox=TRUE)

#'Read JSON file
geographicCoverage = jsonlite::read_json("meta/terrestrial_geo.json")

#'Get start and end dates
temporalCoverage <- list(rangeOfDates =
                           list(beginDate = list(calendarDate = min(forecast_saved$time)),
                                endDate = list(calendarDate = max(forecast_saved$time))))

#'Create the coverage EML
coverage <- list(geographicCoverage = geographicCoverage,
                 temporalCoverage = temporalCoverage)

#'Create the dataset EML
dataset <- eml$dataset(
  title = "Daily persistence null forecast for nee and lee",
  creator = team_list,
  contact = list(references=team_list[[1]]$id),
  pubDate = as_date(forecast_issue_time),
  intellectualRights = "https://creativecommons.org/licenses/by/4.0/",
  dataTable = dataTable,
  coverage = coverage
)

#'Create extra metadata required for submissions to the challenge
#'The metadata follows the EFI Forecasting Standards
additionalMetadata <- eml$additionalMetadata(
  metadata = list(
    forecast = list(
      # Basic elements
      timestep = "1 day", # should be udunits parsable; already in coverage -> temporalCoverage?
      forecast_horizon = "35 days",
      forecast_issue_time = forecast_issue_time,
      forecast_iteration_id = forecast_iteration_id,
      forecast_project_id = team_name,
      metadata_standard_version = "0.3",
      model_description = list(
        forecast_model_id = forecast_model_id,
        name = "state-space Bayesian null",
        type = "empirical",
        repository = "https://github.com/eco4cast/neon4cast-terrestrial"
      ),
      # MODEL STRUCTURE & UNCERTAINTY CLASSES
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

#'Build full EML
my_eml <- eml$eml(dataset = dataset,
                  additionalMetadata = additionalMetadata,
                  packageId = forecast_iteration_id , 
                  system = "datetime"  ## system used to generate packageId
)

#'Check base EML
EML::eml_validate(my_eml)

#'Check that EML matches EFI Standards
EFIstandards::forecast_validator(my_eml)

#'Write metadata
meta_data_filename <-  paste0(forecast_file_name_base,".xml")
write_eml(my_eml, meta_data_filename)

#'Publish the forecast automatically.  Run only on EFI Challenge server
if(efi_server){
  source("../neon4cast-shared-utilities/publish.R")
  publish(code = "03_terrestrial_flux_daily_null.R",
          data_in = "terrestrial-daily-targets.csv.gz",
          data_out = forecast_file,
          meta = meta_data_filename,
          prefix = "terrestrial/",
          bucket = "forecasts")
}
