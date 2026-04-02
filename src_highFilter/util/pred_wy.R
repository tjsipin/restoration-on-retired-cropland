#' Predict suitability by month
#'
#' This function reads in a folder of monthly averages for environmental data, 
#' then maps predicted suitability from a Maxent model by month. 
#' The outputs are 12 rasters, one for each month,
#' matching the extent, resolution, and crs of the input files. 
#' 
#' @author Nick McManus, TJ Sipin
#' @param model the Maxent model used for predicting suitability
#' @param sp the species being modeled
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


## WY pred fxn --------------------------------------------------------
wy_pred_fxn <- function(model, sp, model_years, bcmPath, soilPath, salinityPath, pathOut) {
    message(paste0(sp, model_years, collapse = " | "))
    out.filename = paste0(pathOut, sp, '_', model_years, '_wy.tif')
    if(file.exists(out.filename)) return(NULL)
    ## Read each variable
    aet <- terra::rast(
        list.files(
            path = bcmPath,
            pattern = paste0("aet"),
            full = TRUE
        )
    ) %>% 
        mean()
    if(length(
        list.files(
            path = bcmPath,
            pattern = paste0("tmx_[a-z]{3}_avg.tif|tmx_[a-z]{3}_2070_2099_MIROC"),
            full = TRUE
        )
    ) != 12) return(NULL)
    tmx <- terra::rast(
        list.files(
            path = bcmPath,
            pattern = paste0("tmx_[a-z]{3}_avg.tif|tmx_[a-z]{3}_2070_2099_MIROC"),
            full = TRUE
        )
    ) %>% 
        mean()
    
    tmn <- terra::rast(
        list.files(
            path = bcmPath,
            pattern = paste0("tmn"),
            full = TRUE
        )
    ) %>% 
        mean()
    tdiff <- tmx - tmn
    tmx10 = rast(
        list.files(
            bcmPath,
            pattern='tmx_oct',
            full.names=T
        )
    )
    tmx11 = rast(
        list.files(
            bcmPath,
            pattern='tmx_nov',
            full.names=T
        )
    )
    tmx12 = rast(
        list.files(
            bcmPath,
            pattern='tmx_dec',
            full.names=T
        )
    )
    tmx1 = rast(
        list.files(
            bcmPath,
            pattern='tmx_jan',
            full.names=T
        )
    )
    tmx2 = rast(
        list.files(
            bcmPath,
            pattern='tmx_feb',
            full.names=T
        )
    )
    tmx3 = rast(
        list.files(
            bcmPath,
            pattern='tmx_mar',
            full.names=T
        )
    )
    tmx4 = rast(
        list.files(
            bcmPath,
            pattern='tmx_apr',
            full.names=T
        )
    )
    tmx5 = rast(
        list.files(
            bcmPath,
            pattern='tmx_may',
            full.names=T
        )
    )
    tmx6 = rast(
        list.files(
            bcmPath,
            pattern='tmx_jun',
            full.names=T
        )
    )
    tmx7 = rast(
        list.files(
            bcmPath,
            pattern='tmx_jul',
            full.names=T
        )
    )
    tmx8 = rast(
        list.files(
            bcmPath,
            pattern='tmx_aug',
            full.names=T
        )
    )
    tmx9 = rast(
        list.files(
            bcmPath,
            pattern='tmx_sep',
            full.names=T
        )
    )
    tmx_summer_mean = rast(list(tmx6, tmx7, tmx8, tmx9)) %>% 
        mean()
    ppt10 = rast(
        list.files(
            bcmPath,
            pattern='ppt_oct',
            full.names=T
        )
    )
    ppt11 = rast(
        list.files(
            bcmPath,
            pattern='ppt_nov',
            full.names=T
        )
    )
    ppt12 = rast(
        list.files(
            bcmPath,
            pattern='ppt_dec',
            full.names=T
        )
    )
    ppt1 = rast(
        list.files(
            bcmPath,
            pattern='ppt_jan',
            full.names=T
        )
    )
    ppt2 = rast(
        list.files(
            bcmPath,
            pattern='ppt_feb',
            full.names=T
        )
    )
    ppt3 = rast(
        list.files(
            bcmPath,
            pattern='ppt_mar',
            full.names=T
        )
    )
    ppt4 = rast(
        list.files(
            bcmPath,
            pattern='ppt_apr',
            full.names=T
        )
    )
    ppt5 = rast(
        list.files(
            bcmPath,
            pattern='ppt_may',
            full.names=T
        )
    )
    ppt6 = rast(
        list.files(
            bcmPath,
            pattern='ppt_jun',
            full.names=T
        )
    )
    ppt7 = rast(
        list.files(
            bcmPath,
            pattern='ppt_jul',
            full.names=T
        )
    )
    ppt8 = rast(
        list.files(
            bcmPath,
            pattern='ppt_aug',
            full.names=T
        )
    )
    ppt9 = rast(
        list.files(
            bcmPath,
            pattern='ppt_sep',
            full.names=T
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
    stack <- rast(list(
        aet, tmx, tdiff, 
        tmx10, tmx11, tmx12,
        tmx1, tmx2, tmx3,
        tmx4, tmx5, tmx6,
        tmx7, tmx8, tmx9,
        tmx_summer_mean,
        ppt10, ppt11, ppt12,
        ppt1, ppt2, ppt3,
        ppt4, ppt5, ppt6,
        ppt7, ppt8, ppt9,
        om, ph, drclass, cec, salinity
    )) %>% 
        project("EPSG:3310")
    names(stack) <- c(
        'aet', 'tmx', 'tdiff', 
        'tmx10', 'tmx11', 'tmx12',
        'tmx1', 'tmx2', 'tmx3',
        'tmx4', 'tmx5', 'tmx6',
        'tmx7', 'tmx8', 'tmx9',
        'tmx_summer_mean',
        'ppt10', 'ppt11', 'ppt12',
        'ppt1', 'ppt2', 'ppt3',
        'ppt4', 'ppt5', 'ppt6',
        'ppt7', 'ppt8', 'ppt9',
        'om', 'ph', 'drclass', 'cec', 'salinity'
    )
    
    #Get the superlative variables ({min/max}_month_{tmx/ppt})
    tmx_stack = rast(list(
        tmx1, tmx2, tmx3, tmx4, tmx5, tmx6,
        tmx7, tmx8, tmx9, tmx10, tmx11, tmx12
    ))
    names(tmx_stack) = c('tmx1', 'tmx2', 'tmx3',
                         'tmx4', 'tmx5', 'tmx6',
                         'tmx7', 'tmx8', 'tmx9',
                         'tmx10', 'tmx11', 'tmx12')
    tmx_stack_super = tmx_stack
    
    tmx_stack_super$tmx_min_month = which.min(tmx_stack)
    tmx_stack_super$tmx_max_month = which.max(tmx_stack)
    
    ppt_stack = rast(list(
        ppt1, ppt2, ppt3, ppt4, ppt5, ppt6,
        ppt7, ppt8, ppt9, ppt10, ppt11, ppt12
    ))
    names(ppt_stack) = c('ppt1', 'ppt2', 'ppt3',
                         'ppt4', 'ppt5', 'ppt6',
                         'ppt7', 'ppt8', 'ppt9',
                         'ppt10', 'ppt11', 'ppt12')
    ppt_stack_super = ppt_stack
    
    ppt_stack_super$ppt_min_month = which.min(ppt_stack)
    ppt_stack_super$ppt_max_month = which.max(ppt_stack)
    
    stack_super = rast(list(
        tmx_min_month = tmx_stack_super$tmx_min_month,
        tmx_max_month = tmx_stack_super$tmx_max_month,
        ppt_min_month = ppt_stack_super$ppt_min_month,
        ppt_max_month = ppt_stack_super$ppt_max_month
    ))
    
    stack = rast(
        list(
            stack, stack_super
        )
    ) %>% 
        tidyterra::select(-c(tmx10:tmx9)) %>% 
        stack() 
    
    ## Predict suitability with model, then save raster
    pred <- predict(model, stack) 
    
    writeRaster(pred, out.filename, overwrite=TRUE)
} ### END `pred_fxn()`

pred_wy <- function(model, sp, model_years, bcmPath, soilPath, salinityPath, pathOut) {
    
    ## List of months for map
    months = tibble(
        mon = c('jan', 'feb', 'mar', 'apr', 'may', 'jun',
                'jul', 'aug', 'sep', 'oct', 'nov', 'dec'),
        mm = 1:12
    ) %>% 
        pull(mon) #%>% 
    # paste(collapse="|")
    
    wy_pred_fxn(
        model= model, sp = sp, model_years = model_years, pathOut = pathOut, bcmPath = bcmPath, soilPath = soilPath, salinityPath = salinityPath
    )
    
} ### END SCRIPT
