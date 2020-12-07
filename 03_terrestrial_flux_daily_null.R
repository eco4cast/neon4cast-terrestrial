#'# Ecological Forecasting Initiative Null Model 

#'## Set-up

print(paste0("Running Creating Daily Terrestrial Forecasts at ", Sys.time()))

#'Load renv.lock file that includes the versions of all the packages used
#'You can generate using the command renv::snapshot()

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
team_name <- "EFInulldaily"

#'Download target file from the server
download.file("https://data.ecoforecast.org/targets/terrestrial/terrestrial_daily-targets.csv.gz",
              "terrestrial_daily-targets.csv.gz")

#'Read in target file.  The guess_max is specified because there could be a lot of
#'NA values at the beginning of the file
terrestrial_targets <- read_csv("terrestrial_daily-targets.csv.gz", guess_max = 10000)

terrestrial_targets <- terrestrial_targets %>% 
  filter(time < as_date("2020-09-01"))

download.file("https://data.ecoforecast.org/targets/terrestrial/terrestrial_30min-targets.csv.gz",
              "terrestrial_30min-targets.csv.gz")

terrestrial_targets_30min <- read_csv("terrestrial_30min-targets.csv.gz", guess_max = 10000)

nee_sd <- (sqrt(2) * terrestrial_targets_30min$nee_sd_intercept) * ((12 / 1000000) * (60 * 60 * 24)) / sqrt(48)
le_sd <- (sqrt(2) * terrestrial_targets_30min$nee_sd_intercept) * ((12 / 1000000) * (60 * 60 * 24)) / sqrt(48)

#'Focal sites
site_names <- c("BART","KONZ","OSBS","SRER")

#'Generic random walk state-space model is JAGS format.  We use this model for 
#'both the NEE and LE null forecasts
RandomWalk = "
model{

  # Priors
  x[1] ~ dnorm(x_ic,tau_init)
  tau_add ~ dgamma(0.1,0.1)
  tau_init ~ dgamma(0.1,0.1)
  
  # Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
    x_obs[t] ~ dnorm(x[t],tau_obs[t])
  }

  # Data Model
  for(i in 1:nobs){
    y[i] ~ dnorm(x[y_wgaps_index[i]], tau_obs[y_wgaps_index[i]])
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
               tau_obs = 1/(nee_sd ^ 2),
               n = length(y_wgaps),
               x_ic = 0.0)
  
  #Initialize parameters 
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(tau_add = 1/var(diff(y_nogaps)),
                      tau_init = mean( 1/var(diff(y_nogaps)), na.rm = TRUE),
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
  jags.out   <- coda.samples(model = j.model,variable.names = c("tau_add","tau_init"), n.iter = 10000)
  
  #Run JAGS model again and sample from the posteriors
  m   <- coda.samples(model = j.model,
                      variable.names = c("x","tau_add","tau_init", "x_obs"),
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
               tau_obs = 1/(le_sd ^ 2),
               n = length(y_wgaps),
               x_ic = 0.0)
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(tau_add = 1/var(diff(y_nogaps)),
                      tau_init = mean(1/var(diff(y_nogaps)), na.rm = TRUE),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i],
                      x = init_x)
  }
  
  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           inits = init,
                           n.chains = 3)
  
  jags.out   <- coda.samples(model = j.model,variable.names = c("tau_add","tau_init"), n.iter = 10000)
  
  m   <- coda.samples(model = j.model,
                      variable.names = c("x","tau_add","tau_init", "x_obs"),
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

vswc_sd <-  rep(mean(terrestrial_targets$vswc_sd, na.rm = TRUE), length(terrestrial_targets$vswc_sd))
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
               tau_obs = 1/(vswc_sd ^ 2),
               n = length(y_wgaps),
               x_ic = 0.3)
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(tau_add = 1/var(diff(y_nogaps)),
                      tau_init = mean(1/var(diff(y_nogaps)), na.rm = TRUE),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i],
                      x = init_x)
  }
  
  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           inits = init,
                           n.chains = 3)
  
  jags.out   <- coda.samples(model = j.model,variable.names = c("tau_add","tau_init"), n.iter = 10000)
  
  m   <- coda.samples(model = j.model,
                      variable.names = c("x","tau_add","tau_init", "x_obs"),
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

ggplot(forecast_saved, aes(x = time, y = vswc, group = ensemble)) + 
  geom_line() +
  facet_wrap(~siteID)

#'#Generate metadata

#'Get system time for setting the issue time of the forecast
curr_time <- with_tz(Sys.time(), "UTC")
#forecast_issue_time <- format(curr_time,format = "%Y-%m-%d %H:%M:%SZ", usetz = F)
forecast_issue_time <- as_date(curr_time)
forecast_iteration_id <- start_forecast

#' The team name is the `forecast_model_id`
forecast_model_id <- team_name

source("generate_metadata.R")

meta_data_filename <- generate_metadata(forecast_file =  forecast_file,
                                        metadata_yaml = "metadata.yml",
                                        forecast_issue_time = as_date(with_tz(Sys.time(), "UTC")),
                                        forecast_iteration_id = start_forecast,
                                        forecast_file_name_base = forecast_file_name_base)


#'Publish the forecast automatically.  Run only on EFI Challenge server
if(efi_server){
  source("../neon4cast-shared-utilities/publish.R")
  publish(code = "03_terrestrial_flux_daily_null.R",
          data_in = "terrestrial_daily-targets.csv.gz",
          data_out = forecast_file,
          meta = meta_data_filename,
          prefix = "terrestrial/",
          bucket = "forecasts")
}
