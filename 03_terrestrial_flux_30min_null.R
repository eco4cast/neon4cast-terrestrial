print(paste0("Running Creating 30-minute Terrestrial Forecasts at ", Sys.time()))

pecan_flux_uncertainty <- "../pecan/modules/uncertainty/R/flux_uncertainty.R"

source(pecan_flux_uncertainty)

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

team_name <- "pnull30min"
forecast_project_id <- "efi_null"

download.file("https://data.ecoforecast.org/targets/terrestrial/terrestrial_30min-targets.csv.gz",
              "terrestrial_30min-targets.csv.gz")

terrestrial_targets <- read_csv("terrestrial_30min-targets.csv.gz", guess_max = 10000)

site_names <- c("BART","KONZ","OSBS","SRER")

RandomWalk = "
model{

  #### Priors
  for(t in 1:48){
    x[t] ~ dnorm(x_ic[t], tau_add)
    tau_obs[t] <- pow(ifelse(x[t] >= 0, (sqrt(2) * (obs_intercept + x[t] * obs_slopeP)), (sqrt(2) * (obs_intercept + x[t] * obs_slopeN))), -2)
    y[t] ~ dnorm(x[t], tau_obs[t])
  }
  

  sd_add  ~ dunif(0.000001, 100)
  tau_add <- 1/ pow(sd_add, 2)

  #### Process Model
  for(t in 49:n){
    x[t]~dnorm(x[t-48], tau_add)
    tau_obs[t] <- pow(ifelse(x[t] >= 0, (sqrt(2) * (obs_intercept + x[t] * obs_slopeP)), (sqrt(2) * (obs_intercept + x[t] * obs_slopeN))), -2)
    y[t] ~ dnorm(x[t],tau_obs[t])
  }


}
"

nee_figures <- list()
forecast_saved_nee <- NULL
for(s in 1:length(site_names)){
  
  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  #unc <- flux.uncertainty(measurement = site_data_var$nee, 
  #                        QC = rep(0, length(site_data_var$nee)),
  #                        bin.num = 20)
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max(site_data_var$time) + days(1)
  
  min_time <- min(which(!is.na(site_data_var$nee)))
  # This is key here - I added 16 days on the end of the data for the forecast period
  full_time <- tibble(time = seq(site_data_var$time[min_time], max(site_data_var$time) + days(35), by = "30 min"))
  
  site_data_var <- left_join(full_time, site_data_var) %>% 
    filter(time >= max_time - months(5))
  
  
  
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
  
  data <- list(y = y_wgaps,
               n = length(y_wgaps),
               x_ic = rep(0.0, 48),
               obs_intercept = site_data_var$nee_sd_intercept[1],
               obs_slopeP = site_data_var$nee_sd_slopeP[1],
               obs_slopeN = site_data_var$nee_sd_slopeN[1]
               )
  
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
                      variable.names = c("sd_add","y"),
                      n.iter = 10000,
                      thin = 5)
  
  
  
  model_output <- m %>%
    spread_draws(y[day]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = site_data_var$time[day]) %>%
    ungroup() %>%
    select(time, y, ensemble)
  
  obs <- tibble(time = site_data_var$time,
                obs = y_wgaps)
  
  #ggplot(obs, aes(x = time, y = obs)) + geom_point()
  
  rm(m)
  gc()
  
  model_output %>% 
    group_by(time) %>% 
    filter(time > lubridate::as_date("2021-01-01")) %>% 
    summarise(mean = mean(y),
              upper = quantile(y, 0.975),
              lower = quantile(y, 0.025),.groups = "drop") %>% 
    ggplot(aes(x = time, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
    #geom_point(data = obs, aes(x = time, y = obs), color = "red") +
    labs(x = "Date", y = "nee", title = site_names[s])
  
  ggsave(paste0("nee_30min_",site_names[s],"_figure.pdf"), device = "pdf")
  
    forecast_saved_tmp <- model_output %>%
      filter(time > start_forecast) %>%
      rename(nee = y) %>% 
      mutate(data_assimilation = 0,
             forecast = 1,
             obs_flag = 2,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = team_name)
    
    forecast_saved_tmp %>% 
      group_by(time) %>% 
      summarise(mean = mean(nee),
                upper = quantile(nee, 0.975),
                lower = quantile(nee, 0.025),.groups = "drop") %>% 
      ggplot(aes(x = time, y = mean)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
      labs(x = "Date", y = "nee", title = site_names[s])
    
    forecast_saved_nee <- rbind(forecast_saved_nee, forecast_saved_tmp)
    
    rm(forecast_saved_tmp)
    gc()
}



# Latent heat

le_figures <- list()

for(s in 1:length(site_names)){
  
  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  unc <- flux.uncertainty(site_data_var$le, QC = rep(0, length(site_data_var$le)))
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max_time
  # This is key here - I added 16 days on the end of the data for the forecast period
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "30 min"))
  
  site_data_var <- left_join(full_time, site_data_var) %>% 
    filter(time >= max_time - months(5) & time < lubridate::as_date("2020-09-14"))
    #filter(time >= max_time - months(5))
  
  ggplot(site_data_var, aes(x = time, y = le)) + geom_point()
  
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
  
  data <- list(y = y_wgaps,
               n = length(y_wgaps),
               x_ic = rep(100, 48),
               obs_intercept = site_data_var$le_sd_intercept[1],
               obs_slopeP = site_data_var$le_sd_slopeP[1],
               obs_slopeN = site_data_var$le_sd_slopeN[1])
  
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
    mutate(time = site_data_var$time[day]) %>%
    ungroup() %>%
    select(time, y, ensemble)
  
  rm(m)
  gc()
  
  obs <- tibble(time = full_time$time,
                obs = y_wgaps)
  
  model_output %>% 
    group_by(time) %>% 
    filter(time >= max_time) %>% 
    summarise(mean = mean(y),
              upper = quantile(y, 0.975),
              lower = quantile(y, 0.025),.groups = "drop") %>% 
    ggplot(aes(x = time, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
    #geom_point(data = obs, aes(x = time, y = obs), color = "red") +
    labs(x = "Date", y = "le", title = site_names[s])
  
  ggsave(paste0("le_30min_",site_names[s],"_figure.pdf"), device = "pdf")
  
  if(s == 1){
    forecast_saved_le <- model_output %>%
      filter(time >= start_forecast) %>%
      rename(le = y) %>% 
      mutate(data_assimilation = 0,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = "EFInull")
  }else{
    
    forecast_saved_tmp <- model_output %>%
      filter(time >= start_forecast) %>%
      rename(le = y) %>% 
      mutate(data_assimilation = 0,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = "EFInull")
    
    forecast_saved_le <- rbind(forecast_saved_le, forecast_saved_tmp)
  }
}

ggsave("le_30min_figures.pdf", le_figures, device = "pdf")

# Soil moisture

soil_moisture_figures <- list()

for(s in 1:length(site_names)){
  
  site_data_var <- terrestrial_targets %>%
    filter(siteID == site_names[s])
  
  max_time <- max(site_data_var$time) + days(1)
  
  start_forecast <- max_time
  # This is key here - I added 16 days on the end of the data for the forecast period
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(35), by = "30 min"))
  
  site_data_var <- left_join(full_time, site_data_var) %>% 
    filter(time >= max_time - months(5))
  
  mean_sd <- mean(site_data_var$vswc_sd, na.rm = TRUE)
  
  site_data_var <- site_data_var %>% 
    mutate(vswc_sd = ifelse(!is.na(vswc_sd), vswc_sd, mean_sd))
  
  # NEE
  
  #Full time series with gaps
  y_wgaps <- site_data_var$vswc
  time <- c(site_data_var$time)
  #Remove gaps
  y_nogaps <- y_wgaps[!is.na(y_wgaps)]
  #Indexes of full time series with gaps
  y_wgaps_index <- 1:length(y_wgaps)
  #keep indexes to reference the gappy time series
  y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
  
  init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
  
  data <- list(y = y_wgaps,
               n = length(y_wgaps),
               x_ic = rep(0.0, 48),
               sd_obs = site_data_var$vswc_sd)
  
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
  
  rm(m)
  gc()
  
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
    labs(x = "Date", y = "vswc", title = site_names[s])
  
  ggsave(paste0("soil_moisture_30min_",site_names[s],"_figure.pdf"), device = "pdf")
  
  if(s == 1){
    forecast_saved_soil_moisture <- model_output %>%
      filter(time >= start_forecast) %>%
      rename(le = y) %>% 
      mutate(data_assimilation = 0,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = "EFInull")
  }else{
    
    forecast_saved_tmp <- model_output %>%
      filter(time >= start_forecast) %>%
      rename(vswc = y) %>% 
      mutate(data_assimilation = 0,
             siteID = site_names[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = "EFInull")
    
    forecast_saved_soil_moisture <- rbind(forecast_saved_soil_moisture, forecast_saved_tmp)
  }
}

forecast_saved <- cbind(forecast_saved_nee, forecast_saved_le$le, forecast_saved_soil_moisture$vswc) %>% 
  rename(le = `forecast_saved_le$le`,
         vswc = `forecast_saved_soil_moisture$vswc`) %>% 
  select(time, ensemble, siteID, obs_flag, nee, le, vswc, forecast, data_assimilation)

forecast_file_name_base <- paste0("terrestrial-",as_date(start_forecast),"-",team_name)
forecast_file <- paste0(forecast_file_name_base, ".csv.gz")
write_csv(forecast_saved, forecast_file)

# Generate metadata

curr_time <- with_tz(Sys.time(), "UTC")
#forecast_issue_time <- format(curr_time,format = "%Y-%m-%d %H:%M:%SZ", usetz = F)
forecast_issue_time <- as_date(curr_time)
forecast_iteration_id <- start_forecast
forecast_model_id <- team_name

source("generate_metadata.R")
meta_data_filename <- generate_metadata(forecast_file =  forecast_file,
                                        metadata_yaml = "metadata_30min.yml",
                                        forecast_issue_time = as_date(with_tz(Sys.time(), "UTC")),
                                        forecast_iteration_id = start_forecast,
                                        forecast_file_name_base = forecast_file_name_base)
## Publish the forecast automatically. (EFI-only)

source("../neon4cast-shared-utilities/publish.R")
publish(code = c("03_terrestrial_flux_30min_null.R", pecan_flux_uncertainty),
        data_in = "terrestrial-30min-targets.csv.gz",
        data_out = forecast_file,
        meta = meta_data_filename,
        prefix = "terrestrial/",
        bucket = "forecasts")



