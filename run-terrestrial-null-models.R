print(paste0("Running Creating Terrestrial Targets at ", Sys.time()))

generate_null_daily <- TRUE
generate_null_30min <- TRUE

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
