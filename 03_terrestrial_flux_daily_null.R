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
library(EFIstandards)
library(EML)
library(jsonlite)

#' set the random number for reproducible MCMC runs
set.seed(329)

#'Generate plot to visualized forecast
generate_plots <- FALSE

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
team_name <- "persistence"

#'Download target file from the server
download.file("https://data.ecoforecast.org/neon4cast-targets/terrestrial_daily/terrestrial_daily-targets.csv.gz",
              "terrestrial_daily-targets.csv.gz")

#'Read in target file.  The guess_max is specified because there could be a lot of
#'NA values at the beginning of the file
terrestrial_targets <- read_csv("terrestrial_daily-targets.csv.gz", guess_max = 10000)

terrestrial_targets |> 
  group_by(variable) |> 
  summarize(mean = quantile(observed, 0.75, na.rm = TRUE))

terrestrial_targets <- terrestrial_targets #%>% 
  #filter(time < as_date("2020-12-01"))

#download.file("https://data.ecoforecast.org/neon4cast-targets/terrestrial_30min/terrestrial_30min-targets.csv.gz",
#              "terrestrial_30min-targets.csv.gz")

#terrestrial_targets_30min <- read_csv("terrestrial_30min-targets.csv.gz", guess_max = 10000)

#variable_sd <- terrestrial_targets_30min |> 
#  mutate(sd = (sqrt(2) * sd_intercept) * ((12 / 1000000) * (60 * 60 * 24)) / sqrt(48)) |> 
#  group_by(site_id, variable) |> 
#  summarize(sd = mean(sd, na.rm = TRUE), .groups = "drop")

#'Focal sites
sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(terrestrial == 1)

site_names <- sites$field_site_id

#'Generic random walk state-space model is JAGS format.  We use this model for 
#'both the NEE and LE null forecasts
RandomWalk = "
model{

  # Priors
  x[1] ~ dnorm(x_ic,tau_add)
  
  sd_add  ~ dunif(0.0000001, 100)
  tau_add <- 1/ pow(sd_add, 2)

  
  # Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
    #Data Model
    y[t] ~ dnorm(x[t],tau_obs)
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
  
  message(paste0("NEE: ", site_names[s]))
  
  #site_sd <- variable_sd |> 
  #  filter(site_id == site_names[s],
  #         variable == "nee") |> 
  #  pull(sd)
    
  
  # Select site
  site_data_var <- terrestrial_targets %>%
    filter(variable == "nee") |> 
    filter(site_id == site_names[s], 
           time >= lubridate::as_date("2020-01-01")) 
  
  # Find the last day in the observed data and add one day for the start of the 
  # forecast
  start_forecast <- max(site_data_var$time) + days(1)
  
  # This is key here - I added 35 days on the end of the data for the forecast period
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "1 day"))
  
  # Join the full time with the site_data_var so there aren't gaps in the time column
  site_data_var <- left_join(full_time, site_data_var)
  
  #observed NEE: Full time series with gaps
  y_wgaps <- site_data_var$observed
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
  data <- list(y = y_wgaps,
               tau_obs = 1/(0.05 ^ 2),
               n = length(y_wgaps),
               x_ic = 0.0)
  
  #Initialize parameters 
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd_add = sd(diff(y_nogaps)),
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
  jags.out   <- coda.samples(model = j.model,variable.names = c("sd_add"), n.iter = 10000)
  
  #Run JAGS model again and sample from the posteriors
  m   <- coda.samples(model = j.model,
                      variable.names = c("y","sd_add"),
                      n.iter = 10000,
                      thin = 5)
  
  #Use TidyBayes package to clean up the JAGS output
  model_output <- m %>%
    spread_draws(y[day]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = full_time$time[day]) %>%
    ungroup() %>%
    select(time, y, ensemble)
  
  if(generate_plots){
    #Pull in the observed data for plotting
    obs <- tibble(time = full_time$time,
                  obs = y_wgaps)
    
    
    #Post past and future
    model_output %>% 
      group_by(time) %>% 
      summarise(mean = mean(y),
                upper = quantile(y, 0.975),
                lower = quantile(y, 0.025),.groups = "drop") %>% 
      ggplot(aes(x = time, y = mean)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
      geom_point(data = obs, aes(x = time, y = obs), color = "red") +
      labs(x = "Date", y = "nee")
    
    ggsave(paste0("nee_daily_",site_names[s],"_figure.pdf"), device = "pdf")
  }
  
  #Filter only the forecasted dates and add columns for required variable
  forecast_saved_tmp <- model_output %>%
    filter(time >= start_forecast) %>%
    rename(predicted = y) %>% 
    mutate(variable = "nee",
           site_id = site_names[s]) %>%
    mutate(forecast_iteration_id = start_forecast) %>%
    mutate(forecast_project_id = team_name)
  
  # Combined with the previous sites
  forecast_saved_nee <- bind_rows(forecast_saved_nee, forecast_saved_tmp)
  
}

#'## Latent heat model
#' 
#' See notes from the NEE section above
#+ message = FALSE


forecast_saved_le <- NULL
le_figures <- list()
for(s in 1:length(site_names)){
  
  message(paste0("LE: ", site_names[s]))
  
  site_data_var <- terrestrial_targets %>%
    filter(variable == "le") |> 
    filter(site_id == site_names[s], 
           time >= lubridate::as_date("2020-01-01"))
    
    #site_sd <- variable_sd |> 
    #  filter(site_id == site_names[s],
    #         variable == "le") |> 
    #  pull(sd)
    
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max_time
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "1 day"))
  
  site_data_var <- left_join(full_time, site_data_var)
  
  y_wgaps <- site_data_var$observed
  time <- c(site_data_var$time)
  
  y_nogaps <- y_wgaps[!is.na(y_wgaps)]
  
  y_wgaps_index <- 1:length(y_wgaps)
  
  y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
  
  init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
  
  data <- list(y = y_wgaps,
               tau_obs = 1/(0.1 ^ 2),
               n = length(y_wgaps),
               x_ic = 0.0)
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd_add = sd(diff(y_nogaps)),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i],
                      x = init_x)
  }
  
  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           inits = init,
                           n.chains = 3)
  
  jags.out   <- coda.samples(model = j.model,variable.names = c("sd_add"), n.iter = 10000)
  
  m   <- coda.samples(model = j.model,
                      variable.names = c("y","sd_add"),
                      n.iter = 10000,
                      thin = 5)
  
  model_output <- m %>%
    spread_draws(y[day]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = full_time$time[day]) %>%
    ungroup() %>%
    select(time, y, ensemble)
  
  if(generate_plots){
    obs <- tibble(time = full_time$time,
                  obs = y_wgaps)
    
    model_output %>% 
      group_by(time) %>% 
      summarise(mean = mean(y),
                upper = quantile(y, 0.975),
                lower = quantile(y, 0.025),.groups = "drop") %>% 
      ggplot(aes(x = time, y = mean)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
      geom_point(data = obs, aes(x = time, y = obs), color = "red") +
      labs(x = "Date", y = "le")
    
    ggsave(paste0("le_daily_",site_names[s],"_figure.pdf"), device = "pdf")
  }
  
  forecast_saved_tmp <- model_output %>%
    filter(time >= start_forecast) %>%
    rename(predicted = y) %>% 
    mutate(variable = "le",
           site_id = site_names[s]) %>%
    mutate(forecast_iteration_id = start_forecast) %>%
    mutate(forecast_project_id = team_name)
  
  forecast_saved_le <- bind_rows(forecast_saved_le, forecast_saved_tmp)
}


#'Combined the NEE and LE forecasts together and re-order column
forecast_saved <- bind_rows(forecast_saved_nee, forecast_saved_le) %>% 
  select(time, site_id, ensemble, variable, predicted)

#'Save file as CSV in the
#'[theme_name]-[year]-[month]-[date]-[team_name].csv
forecast_file_name_base <- paste0("terrestrial_daily-",as_date(Sys.Date()),"-",team_name)
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

#source("generate_metadata.R")

#meta_data_filename <- generate_metadata(forecast_file =  forecast_file,
#                                        metadata_yaml = "metadata.yml",
#                                        forecast_issue_time = as_date(with_tz(Sys.time(), "UTC")),
#                                        forecast_iteration_id = start_forecast,
#                                        forecast_file_name_base = forecast_file_name_base)


neon4cast::submit(forecast_file = forecast_file, 
                  metadata = NULL, 
                  ask = FALSE)

unlink(forecast_file)
#unlink(meta_data_filename)