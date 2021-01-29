

generate_metadata <- function(forecast_file, metadata_yaml, forecast_issue_time, forecast_iteration_id, forecast_file_name_base,start_time=NULL,stop_time=NULL){
  
  metadata <- yaml::read_yaml(metadata_yaml)
  
  if(tools::file_ext(forecast_file) == "csv"){
    forecast <- readr::read_csv(forecast_file)
  } else {
    forecast <- NULL
  }
  
  theme <- unlist(stringr::str_split(stringr::str_split(forecast_file, "-")[[1]][1], "_")[[1]][1])
  team_name <- unlist(stringr::str_split(unlist(stringr::str_split(forecast_file, "-"))[5], ".csv"))[1]
  
  if(theme == "terrestrial"){
    #'Build attribute table. Other models won't likely change this for the challenge
    #'note: need to fix the units for nee and le because the units do not pass EML unit checks
    attributes <- tibble::tribble(
      ~attributeName,     ~attributeDefinition,                          ~unit,                  ~formatString,  ~definition, ~numberType,
      "time",              "[dimension]{time}",                          "year",                 "YYYY-MM-DD",   NA,          "datetime",
      "ensemble",          "[dimension]{index of ensemble member}",      "dimensionless",         NA,            NA,          "integer",
      "siteID",             "[dimension]{neon site}",                     NA,                     NA,           "NEON site ID",  "character",
      "obs_flag",          "[flag]{observation error}",                  "dimensionless",         NA,           NA,           "integer",
      "nee",               "[variable]{net ecosystem exchange}",         "dimensionless",         NA,           NA,           "real",
      "le",                "[variable]{latent heat}",                    "dimensionless",         NA,           NA,          "real",
      "vswc",              "[variable]{volumetric soil water content}",  "dimensionless",         NA,           NA,           "real",
      "forecast",          "[flag]{whether represents forecast}",        "dimensionless",         NA,           NA,          "integer",
      "data_assimilation", "[flag]{whether time step assimilated data}", "dimensionless",         NA,           NA,          "integer"
    )
    
    if(!is.null(forecast)){
      if(all.equal(names(forecast), attributes$attributeName) != TRUE){
        message("Column names in file do not match required names for complete metadate")
        message(paste0("File names are: ",names(forecast)))
        message(paste0("Required names are: ",attributes$attributeName))
        stop()
      }
    }
    
    #' use `EML` package to build the attribute list
    attrList <- EML::set_attributes(attributes, 
                                    col_classes = c("Date", "numeric", "character","numeric","numeric", 
                                                    "numeric","numeric", "numeric","numeric"))
    
    entityDescription_text = "Forecast of NEE and LE for four NEON sites"
  }
  
  #' use `EML` package to build the physical list
  physical <- EML::set_physical(forecast_file)
  
  #' use `EML` package to dataTable
  dataTable <- eml$dataTable(
    entityName = "forecast",  ## this is a standard name to allow us to distinguish this entity from 
    entityDescription = entityDescription_text,
    physical = physical,
    attributeList = attrList)
  
  fullgeographicCoverage <- jsonlite::read_json("meta/terrestrial_geo.json")
  
  site_id_index <- NULL
  for(i in 1:length(fullgeographicCoverage)){
    if(fullgeographicCoverage[[i]]$id %in% metadata$sites)
      site_id_index <- c(site_id_index, i)
  }
  
  geographicCoverage <- fullgeographicCoverage[site_id_index]
  
  if(is.null(start_time)) start_time = min(forecast$time)
  if(is.null(stop_time)) stop_time = max(forecast$time)
  
  temporalCoverage <- list(rangeOfDates =
                             list(beginDate = list(calendarDate = start_time),
                                  endDate = list(calendarDate = stop_time)))
  #'Create the coverage EML
  coverage <- list(geographicCoverage = geographicCoverage,
                   temporalCoverage = temporalCoverage)
  
  
  
  
  #'Create the dataset EML
  dataset <- eml$dataset(
    title = "Daily persistence null forecast for nee and lee",
    creator = metadata$team_list,
    contact = metadata$team_list[[1]],
    pubDate = as_date(forecast_issue_time),
    intellectualRights = "https://creativecommons.org/licenses/by/4.0/",
    dataTable = dataTable,
    coverage = coverage
  )
  
  metadata$metadata$forecast$forecast_issue_time <- as_date(forecast_issue_time)
  metadata$metadata$forecast$forecast_iteration_id <- forecast_iteration_id
  metadata$metadata$forecast$forecast_project_id <- team_name
  
  
  my_eml <- eml$eml(dataset = dataset,
                    additionalMetadata = eml$additionalMetadata(metadata = metadata$metadata),
                    packageId = forecast_iteration_id , 
                    system = "datetime"  ## system used to generate packageId
  )
  
  
  #'Check base EML
  if(!EML::eml_validate(my_eml)){
    message("Error in EML metadata")
  }
  
  #'Check that EML matches EFI Standards
  if(EFIstandards::forecast_validator(my_eml)){
    #'Write metadata
    meta_data_filename <-  paste0(forecast_file_name_base,".xml")
    EML::write_eml(my_eml, meta_data_filename)
  }else{
    message("Error in EFI metadata")
  }
  return(meta_data_filename)
}

