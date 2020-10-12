library(tidyverse)
library(lubridate)
library(rjags)
library(tidybayes)
library(modelr)
library(nimble)
library(coda)

set.seed(329)

download.file("https://data.ecoforecast.org/targets/terrestrial_fluxes/terrestrial-30min-targets.csv.gz",
              "terrestrial-30min-targets.csv.gz")

terrestrial_targets <- read_csv("terrestrial-30min-targets.csv.gz", guess_max = 10000)

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

nee_figures <- list()
for(s in 1:length(site_names)){
  
  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max_time
  # This is key here - I added 16 days on the end of the data for the forecast period
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "30 min"))
  
  site_data_var <- left_join(full_time, site_data_var)
  
  # NEE
  
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
  
  nee_figures[s] <- model_output %>% 
    group_by(time) %>% 
    summarise(mean = mean(x_obs),
              upper = quantile(x_obs, 0.975),
              lower = quantile(x_obs, 0.025),.groups = "drop") %>% 
    ggplot(aes(x = time, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
    geom_point(data = obs, aes(x = time, y = obs), color = "red") +
    labs(x = "Date", y = "nee", title = site_names[s])
  
  if(s == 1){
    forecast_saved_nee <- model_output %>%
      filter(time > start_forecast) %>%
      rename(nee = x_obs) %>% 
      mutate(data_assimilation = 0,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = "EFInull")
  }else{
    
    forecast_saved_tmp <- model_output %>%
      filter(time > start_forecast) %>%
      rename(nee = x_obs) %>% 
      mutate(data_assimilation = 0,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = "EFInull")
    
    forecast_saved_nee <- rbind(forecast_saved_nee, forecast_saved_tmp)
  }
}

# Latent heat

le_figures <- list()

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
  
  le_figures[s] <- model_output %>% 
    group_by(time) %>% 
    summarise(mean = mean(x_obs),
              upper = quantile(x_obs, 0.975),
              lower = quantile(x_obs, 0.025),.groups = "drop") %>% 
    ggplot(aes(x = time, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
    geom_point(data = obs, aes(x = time, y = obs), color = "red") +
    labs(x = "Date", y = "le", title = site_names[s])
  
  if(s == 1){
    forecast_saved_le <- model_output %>%
      filter(time > start_forecast) %>%
      rename(le = x_obs) %>% 
      mutate(data_assimilation = 0,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = "EFInull")
  }else{
    
    forecast_saved_tmp <- model_output %>%
      filter(time > start_forecast) %>%
      rename(le = x_obs) %>% 
      mutate(data_assimilation = 0,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = "EFInull")
    
    forecast_saved_le <- rbind(forecast_saved_le, forecast_saved_tmp)
  }
}

forecast_saved <- cbind(forecast_saved_nee, forecast_saved_le$le) %>% 
  rename(le = `forecast_saved_le$le`) %>% 
  select(time, siteID, ensemble, nee, le,data_assimilation, forecast_iteration_id, forecast_project_id)

forecast_saved %>% 
  select(time, nee, )


forecast_file_name <- paste0("terrestrial-EFInull30min-",as_date(start_forecast),".csv.gz")
write_csv(forecast_saved, forecast_file_name)

source("../EFI_terrestrial/R/publish.R")
publish(code = "03_terrestrial_flux_daily_null",
        data_in = "terrestrial-daily-targets.csv.gz",
        data_out = forecast_file_name,
        prefix = "terrestrial_fluxes/",
        bucket = "forecasts")