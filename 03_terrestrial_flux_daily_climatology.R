#'# Ecological Forecasting Initiative Null Model 

#'## Set-up

print(paste0("Running Creating Daily Terrestrial Forecasts at ", Sys.time()))

#'Load renv.lock file that includes the versions of all the packages used
#'You can generate using the command renv::snapshot()

#' Required packages.  
#' EFIstandards is at remotes::install_github("eco4cast/EFIstandards")
library(tidyverse)
library(lubridate)
library(aws.s3)
library(EFIstandards)
library(EML)
library(jsonlite)
library(imputeTS)


#' set the random number for reproducible MCMC runs
set.seed(329)

#'Generate plot to visualized forecast
generate_plots <- TRUE
#'Is the forecast run on the Ecological Forecasting Initiative Server?
#'Setting to TRUE published the forecast on the server.
efi_server <- TRUE

#' List of team members. Used in the generation of the metadata
#team_list <- list(list(individualName = list(givenName = "Quinn", surName = "Thomas"), 
#                       id = "https://orcid.org/0000-0003-1282-7825"),
#                  list(individualName = list(givenName = "Others",  surName ="Pending")),
#)

#'Team name code
team_name <- "climatology"

#'Read in target file.  The guess_max is specified because there could be a lot of
#'NA values at the beginning of the file
targets <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/terrestrial_daily/terrestrial_daily-targets.csv.gz", guess_max = 10000)

sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(terrestrial == 1)

site_names <- sites$field_site_id

target_clim <- targets %>%  
  mutate(doy = yday(time)) %>% 
  group_by(doy, site_id, variable) %>% 
  summarise(clim_mean = mean(observed, na.rm = TRUE),
            clim_sd = sd(observed, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(clim_mean = ifelse(is.nan(clim_mean), NA, clim_mean))

#curr_month <- month(Sys.Date())
curr_month <- month(Sys.Date())
if(curr_month < 10){
  curr_month <- paste0("0", curr_month)
}

curr_year <- year(Sys.Date())
start_date <- Sys.Date() + days(1)

forecast_dates <- seq(start_date, as_date(start_date + days(35)), "1 day")
forecast_doy <- yday(forecast_dates)

forecast <- target_clim %>%
  mutate(doy = as.integer(doy)) %>% 
  filter(doy %in% forecast_doy) %>% 
  mutate(time = as_date(ifelse(doy > last(doy),
                               as_date((doy-1), origin = paste(year(Sys.Date())+1, "01", "01", sep = "-")),
                               as_date((doy-1), origin = paste(year(Sys.Date()), "01", "01", sep = "-")))))

subseted_site_names <- unique(forecast$site_id)
site_vector <- NULL
for(i in 1:length(subseted_site_names)){
  site_vector <- c(site_vector, rep(subseted_site_names[i], length(forecast_dates)))
}

forecast_tibble1 <- tibble(time = rep(forecast_dates, length(subseted_site_names)),
                          site_id = site_vector,
                          variable = "nee")

forecast_tibble2 <- tibble(time = rep(forecast_dates, length(subseted_site_names)),
                          site_id = site_vector,
                          variable = "nee")

forecast_tibble <- bind_rows(forecast_tibble1, forecast_tibble2)

foreast <- right_join(forecast, forecast_tibble)

combined <- forecast %>% 
  select(time, site_id, variable, clim_mean, clim_sd) %>% 
  rename(mean = clim_mean,
         sd = clim_sd) %>% 
  group_by(site_id, variable) %>% 
  mutate(mu = imputeTS::na_interpolation(x = mean),
         sigma = median(sd, na.rm = TRUE))
combined <- combined %>%
  pivot_longer(c("mu", "sigma"),names_to = "parameter", values_to = "predicted") |> 
  mutate(family = "normal") |> 
  mutate(start_time = min(combined$time) - lubridate::days(1)) |> 
  select(time, start_time, site_id, variable, family, parameter, predicted)

combined %>% 
  filter(variable == "nee") |> 
  pivot_wider(names_from = parameter, values_from = predicted) %>% 
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin=mu - sigma*1.96, ymax=mu + sigma*1.96), alpha = 0.1) + 
  geom_point(aes(y = mu)) +
  facet_wrap(~site_id)

forecast_file <- paste("terrestrial_daily", min(combined$time), "climatology.csv.gz", sep = "-")

write_csv(combined, forecast_file)

neon4cast::submit(forecast_file = forecast_file, 
                  metadata = NULL, 
                  ask = FALSE)

unlink(forecast_file)



