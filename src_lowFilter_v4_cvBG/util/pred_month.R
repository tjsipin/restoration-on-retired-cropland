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
    bd_POLARIS_filename = 'data/0_env/polaris/rasters/bd_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif'
    clay_POLARIS_filename =  'data/0_env/polaris/rasters/clay_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif'
    om_POLARIS_filename = 'data/0_env/polaris/rasters/om_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif'
    ph_POLARIS_filename = 'data/0_env/polaris/rasters/ph_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif'
    silt_POLARIS_filename = 'data/0_env/polaris/rasters/silt_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif'
    
    
    ## Read each variable
    aet <- raster::raster(aet_filename)  
    tmx <- raster::raster(tmx_filename)
    tdiff <- raster::raster(tdiff_filename) 
    tmxSummer <- raster::raster(tmxSummer_filename)
    pptWinter <- raster::raster(pptWinter_filename)
    
    bd_POLARIS <- raster::raster(bd_POLARIS_filename) %>% 
        resample(aet)
    clay_POLARIS <- raster::raster(clay_POLARIS_filename) %>% 
        resample(aet)
    om_POLARIS <- raster::raster(om_POLARIS_filename) %>% 
        resample(aet)
    ph_POLARIS <- raster::raster(ph_POLARIS_filename) %>% 
        resample(aet)
    silt_POLARIS <- raster::raster(silt_POLARIS_filename) %>% 
        resample(aet)
    
    ## Stack and match lyr names to model variables
    stack <- raster::stack(aet, tdiff, tmx, tmxSummer, pptWinter, bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS)
    names(stack) <- c(
        'aet', 'tdiff', 'tmx', 'tmx_summer_mean', 'ppt_winter_sum', 
        'bd_POLARIS', 'clay_POLARIS', 'om_POLARIS', 'ph_POLARIS', 'silt_POLARIS',
    )
    ## Predict suitability with model, then save raster
    pred <- dismo::predict(model, stack)
    # pred_new <- dismo::predict(model, stack_new)
    writeRaster(pred, paste0(pathOut, month, '_', spp, '_', model_years, '.tif'), overwrite=TRUE)
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