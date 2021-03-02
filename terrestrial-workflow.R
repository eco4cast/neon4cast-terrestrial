print(paste0("Running Creating Terrestrial Targets at ", Sys.time()))

Sys.setenv("NEONSTORE_HOME" = "/efi_neon_challenge/neonstore")

remotes::install_deps()

library(neonUtilities)
library(neonstore)
library(tidyverse)
library(lubridate)
library(contentid)

run_full_workflow <- TRUE
generate_null_daily <- TRUE
generate_null_30min <- TRUE

# Terrestrial
#DP4.00200.001 & DP1.00094.001
sites <- c("BART", "KONZ", "SRER", "OSBS")
start_date <- NA

print("Downloading: DP4.00200.001")
new_data1 <- neonstore::neon_download(product = "DP4.00200.001", site = sites, type = "basic", start_date = start_date, .token = Sys.getenv("NEON_TOKEN"))
print("Downloading: DP1.00094.001")
new_data2 <- neonstore::neon_download(product = "DP1.00094.001", site = sites, type = "basic", start_date = start_date, .token = Sys.getenv("NEON_TOKEN"))

if(!is.null(new_data1) | !is.null(new_data2) | run_full_workflow){
  
  source("02_terrestrial_targets.R")
  
  print(paste0("Completed Target at ", Sys.time()))
  
  if(generate_null_daily){
    
    print(paste0("Running daily Null at ", Sys.time()))
    source("03_terrestrial_flux_daily_null.R")
    print(paste0("Completed daily Null at ", Sys.time()))
  }
  
  if(generate_null_30min){
    print(paste0("Running 30 min Null at ", Sys.time()))
    source("04_terrestrial_flux_30min_clim.R")
    print(paste0("Completed 30 min Null at ", Sys.time()))  
  }
}
