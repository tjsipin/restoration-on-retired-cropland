#' Predict suitability by month
#'
#' This function reads in a folder of monthly averages for environmental data, 
#' then maps predicted suitability from a Maxent model by month. 
#' The outputs are 12 rasters, one for each month,
#' matching the extent, resolution, and crs of the input files. 
#' 
#' @author Nick McManus
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
pred_fxn <- function(model, spp, month, model_years, bcmPath, soilPath, pathOut) {
    out.filename = paste0(pathOut, month, '_', spp, '_', model_years, '.tif')
    # if(file.exists(out.filename)) return(NULL)
    
    message(paste0(spp, month, model_years, collapse = " | "))
    aet_filename = list.files(path = bcmPath,
                              pattern = paste0("aet", ".+", month),
                              full = TRUE)
    tmx_filename = list.files(path = bcmPath,
                              pattern = paste0("tmx",".+", month),
                              full = TRUE)
    tdiff_filename = list.files(path = bcmPath,
                                pattern = paste0("tdiff",".+", month),
                                full = TRUE)
    tmxSummer_filename = list.files(path = bcmPath, 
                                    pattern = "tmx_summer",
                                    full = TRUE)
    pptWinter_filename = list.files(path = bcmPath, 
                                    pattern = "ppt_winter",
                                    full = TRUE)
    
    om_filename = list.files(path = soilPath,
                             pattern = "natsgo_om_270m_CA_2023",
                             full = TRUE)
    ph_filename = list.files(path = soilPath,
                             pattern = "natsgo_ph_270m_CA_2023",
                             full = TRUE)
    drclass_filename = list.files(path = soilPath,
                                  pattern = "natsgo_drclass_270m_CA_2023",
                                  full = TRUE)
    cec_filename = list.files(path = soilPath,
                              pattern = "natsgo_cec_270m_CA_2023",
                              full = TRUE)
    
    salinity_filename = 'data/0_env/salinity/salinity_full_CA_res.tif'
    
    ## Read each variable
    aet <- terra::rast(aet_filename)  
    tmx <- terra::rast(tmx_filename)
    tdiff <- terra::rast(tdiff_filename) 
    tmxSummer <- terra::rast(tmxSummer_filename)
    pptWinter <- terra::rast(pptWinter_filename)
    
    om <- terra::rast(om_filename) %>% 
        resample(aet)
    ph <- terra::rast(ph_filename) %>% 
        resample(aet)
    drclass <- terra::rast(drclass_filename) %>% 
        resample(aet, method='mode')
    cec <- terra::rast(cec_filename) %>% 
        resample(aet)
    
    salinity <- terra::rast(salinity_filename) %>% 
        resample(aet, method='mode')
    
    ## Stack and match lyr names to model variables
    stack <- raster::stack(rast(list(aet, tdiff, tmx, tmxSummer, pptWinter, om, ph, drclass, cec, salinity)))
    names(stack) <- c(
        'aet', 'tdiff', 'tmx', 'tmx_summer_mean', 'ppt_winter_sum', 
        'om', 'ph', 'drclass', 'cec', 'salinity'
    )
    ## Predict suitability with model, then save raster
    pred <- dismo::predict(model, stack)
    # pred_new <- dismo::predict(model, stack_new)
    writeRaster(pred, out.filename, overwrite=TRUE)
} ### END `pred_fxn()`

pred_month <- function(model, spp, model_years, bcmPath, soilPath, pathOut) {
    
    ## List of months for map
    month <- c('jan', 'feb', 'mar', 'apr', 'may', 'jun',
               'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
    
    ## Iterate fxn over list of months 
    purrr::map(
        .x=month,
        ~ pred_fxn(
            model = model, spp=spp, model_years = model_years, 
            month = .x, bcmPath = bcmPath, soilPath = soilPath, pathOut = pathOut
        )
    )
    
} ### END SCRIPT