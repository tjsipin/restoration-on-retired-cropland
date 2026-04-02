#' Predict suitability by month
#'
#' This function reads in a folder of monthly averages for environmental data, 
#' then maps predicted suitability from a Maxent model by month. 
#' The outputs are 12 rasters, one for each month,
#' matching the extent, resolution, and crs of the input files. 
#' 
#' @author Nick McManus, TJ Sipin
#' @param model the Maxent model used for predicting suitability
#' @param spp the species being modeled
#' @param model_years the model and time range of supplied env variables. This is only used for file naming convention, so best to include future model with years for clarity (e.g. "2000_2022" or "MIROC45_2070_2099"). 
#' @param bcmPath path to directory of bcm monthly average rasters
#' @param soilPath path to directory of soil rasters
#' @param pathOut path to output directory where new rasters will be saved
#' @return 12 TIFs with predicted habitat suitability by month


get_raster_safe <- function(path, pattern) {
    files <- list.files(path = path, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
        stop(paste("No raster file found matching pattern:", pattern, "in", path))
    }
    raster::raster(files[1])
}
## Monthly pred fxn --------------------------------------------------------
pred_fxn <- function(model, spp, months, model_years, bcmPath, soilPath, salinityPath, pathOut) {
    message(paste0(spp, months, model_years, collapse = " | "))
    ## Read each variable
    aet <- raster::raster(list.files(path = bcmPath,
                                     pattern = paste0("aet", ".+", month),
                                     full = TRUE))  
    tmx <- raster::raster(list.files(path = bcmPath,
                                     pattern = paste0("tmx",".+", month),
                                     full = TRUE))
    tdiff <- raster::raster(list.files(path = bcmPath,
                                       pattern = paste0("tdiff",".+", month),
                                       full = TRUE)) 
    tmxSummer <- raster::raster(list.files(path = bcmPath, 
                                           pattern = "tmx_summer",
                                           full = TRUE))
    pptWinter <- raster::raster(list.files(path = bcmPath, 
                                           pattern = "ppt_winter",
                                           full = TRUE))
    om <- raster::raster(list.files(path = soilPath,
                                    pattern = "natsgo_om_270m_CA_2023",
                                    full = TRUE)) %>% 
        resample(aet)
    ph <- raster::raster(list.files(path = soilPath,
                                    pattern = "natsgo_ph_270m_CA_2023",
                                    full = TRUE)) %>% 
        resample(aet)
    drclass <- raster::raster(list.files(path = soilPath,
                                         pattern = "natsgo_drclass_270m_CA_2023",
                                         full = TRUE)) %>% 
        resample(aet, method='mode')
    cec <- raster::raster(list.files(path = soilPath,
                                     pattern = "natsgo_cec_270m_CA_2023",
                                     full = TRUE)) %>% 
        resample(aet)
    salinity <- raster::raster(list.files(path = salinityPath,
                                          pattern = 'salinity_full_CA_res.tif',
                                          full = TRUE)) %>% 
        resample(aet, method='mode')
    
    ## Stack and match lyr names to model variables
    stack <- raster::stack(aet, tdiff, tmx, tmxSummer, pptWinter, om, ph, drclass, cec)
    names(stack) <- c(
        'aet', 'tdiff', 'tmx', 'tmx_summer_mean', 'ppt_winter_sum', 
        "om", "ph", "drclass", "cec", 'salinity'
    )
    
    ## Predict suitability with model, then save raster
    pred <- dismo::predict(model, stack) 
    
    writeRaster(pred, paste0(pathOut, month, '_', spp, '_', model_years, '.tif'), overwrite=TRUE)
} ### END `pred_fxn()`

pred_month <- function(model, spp, model_years, bcmPath, soilPath, salinityPath = salinityPath, pathOut) {
    
    ## List of months for map
    month <- c('jan', 'feb', 'mar', 'apr', 'may', 'jun',
               'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
    
    ## Iterate fxn over list of months 
    purrr::map(
        .x=month,
        ~ pred_fxn(
            model = model, spp=spp, model_years = model_years, 
            month = .x, bcmPath = bcmPath, soilPath = soilPath, salinityPath = salinityPath, pathOut = pathOut
        )
    )
    
} 



# Aggregated pred function ------------------------------------------------


## Monthly pred fxn --------------------------------------------------------
agg_pred_fxn <- function(model, spp, months, model_years, bcmPath, soilPath, salinityPath, pathOut) {
    message(paste0(spp, months, model_years, collapse = " | "))
    ## Read each variable
    aet <- terra::rast(
        list.files(
            path = bcmPath,
            pattern = paste0("aet", ".+", months, collapse="|"),
            full = TRUE
        )
    ) %>% 
        mean()
    tmx <- terra::rast(
        list.files(
            path = bcmPath,
            pattern = paste0("tmx", ".+", months, collapse="|"),
            full = TRUE
        )
    ) %>% 
        mean()
    
    tmn <- terra::rast(
        list.files(
            path = bcmPath,
            pattern = paste0("tmn", ".+", months, collapse="|"),
            full = TRUE
        )
    ) %>% 
        mean()
    tdiff <- tmx - tmn
    tmxSummer <- terra::rast(
        list.files(
            path = bcmPath, 
            pattern = "tmx_summer",
            full = TRUE
        )
    ) 
    pptWinter <- terra::rast(
        list.files(
            path = bcmPath, 
            pattern = "ppt_winter",
            full = TRUE
        )
    )
    om <- terra::rast(
        list.files(
            path = soilPath,
            pattern = "natsgo_om_270m_CA_2023",
            full = TRUE
        )
    ) %>% 
        resample(aet)
    ph <- terra::rast(
        list.files(
            path = soilPath,
            pattern = "natsgo_ph_270m_CA_2023",
            full = TRUE
        )
    ) %>% 
        resample(aet)
    drclass <- terra::rast(
        list.files(
            path = soilPath,
            pattern = "natsgo_drclass_270m_CA_2023",
            full = TRUE
        )
    ) %>% 
        resample(aet, method='mode')
    cec <- terra::rast(
        list.files(
            path = soilPath,
            pattern = "natsgo_cec_270m_CA_2023",
            full = TRUE
        )
    ) %>% 
        resample(aet)
    salinity <- terra::rast(
        list.files(
            path = salinityPath,
            pattern = "salinity_full_CA_res.tif",
            full = TRUE
        )
    ) %>% 
        resample(aet, method='mode')
    
    ## Stack and match lyr names to model variables
    stack <- rast(list(aet, tdiff, tmx, tmxSummer, pptWinter, om, ph, drclass, cec, salinity)) %>% 
        project("EPSG:3310") %>% 
        stack() 
    names(stack) <- c(
        'aet', 'tdiff', 'tmx', 'tmx_summer_mean', 'ppt_winter_sum', 
        "om", "ph", "drclass", "cec", "salinity"
    )
    
    ## Predict suitability with model, then save raster
    pred <- predict(model, stack) 
    
    writeRaster(pred, paste0(pathOut, spp, '_', model_years, '_agg.tif'), overwrite=TRUE)
} ### END `pred_fxn()`

pred_agg <- function(model, months_vec, spp, model_years, bcmPath, soilPath, salinityPath, pathOut) {
    
    ## List of months for map
    months = tibble(
        mon = c('jan', 'feb', 'mar', 'apr', 'may', 'jun',
                'jul', 'aug', 'sep', 'oct', 'nov', 'dec'),
        mm = 1:12
    ) %>% 
        filter(mm %in% months_vec) %>% 
        pull(mon) #%>% 
    # paste(collapse="|")
    
    agg_pred_fxn(
        model= model, spp = spp, model_years = model_years, months = months, pathOut = pathOut, bcmPath = bcmPath, soilPath = soilPath, salinityPath = salinityPath
    )
    
} ### END SCRIPT
