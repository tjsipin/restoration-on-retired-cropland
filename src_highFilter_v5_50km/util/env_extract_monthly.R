#' Environmental data extraction for points of species occurrence
#'
#' This function reads in a folder of monthly environmental data (in .tif format)
#' from the Basin Characterization Model version 8 (BCMv8) and NATSGO soil data,
#' then filters species occurrence data for that specific month, and then extracts 
#' raster values for each point. The resulting data frame is in a
#' Samples with Data (SWD) format for performing a species distribution model (SDM)
#' with Maxent. 
#' Note: This code was written to be compatible with the specific naming 
#' convention of BCMv8 data. 
#' 
#' @author Nick McManus
#' @param startYear character or numeric. The first water year in dataset.
#' @param endYear character or numeric. The last water year in dataset.
#' @param pathMonth character. File path to directory with monthly rasters.
#' @param pathQuarter character. File path to directory with quarterly data (generated from `quarterly_rast()` function).
#' @param pathSoil character. File path to directory with NATSGO soil rasters.
#' @param occ data frame. Contains species occurrence or background points data for extraction.
#' @param lon character. Variable name for longitude values in occ df (default = "lon")
#' @param lat character. Variable name for the latitude values in occ df (default = "lat")
#' @param crs character. The coordinate reference system of the spp occurrence data (default = "WGS84")
#' @return data frame with environmental values at point of species observation


extractEnvMonthly <- function(startYear, endYear, pathMonth, pathQuarter, pathNatsgo, pathPolaris, pathSalinity,
                              occ, lon = "lon", lat = "lat", crs = "WGS84") {
    
    # CA <- urbnmapr::get_urbn_map(
    #     map = 'counties',
    #     sf = T
    # ) %>% 
    #     terra::vect() %>% 
    #     terra::project("EPSG:4326") %>% 
    #     filter(state_name == "California")
    CA <- tigris::states() %>% 
        vect() %>% 
        filter(STUSPS=="CA") %>% 
        project("EPSG:3310")
    
    ## Warnings
    if (startYear < 2000)
        warning("This function was built with a dataset starting in water year 2000. \nEnsure your start year matches with available data.")
    if (endYear >2023)
        warning("This function was built with a dataset ending at water year 2023.\nEnsure your end year matches with available data.")
    
    
    ## Create df for range of dates -----------------------------------------
    dates_df <- data.frame(wy = rep(startYear:endYear, each = 12), 
                           ### Repeat each yr 12 times to create a row 
                           ### for each mo/yr
                           mon = c('oct', 'nov', 'dec', 'jan', 'feb', 'mar', 
                                   'apr', 'may', 'jun', 'jul', 'aug', 'sep'),
                           mon_num = rep(c(10:12, 1:9), each=1)) %>% 
        ## generate calendar year from wy 
        mutate(year = case_when(mon_num %in% c(10, 11, 12) ~(.$wy-1),
                                .default = .$wy))
    t0=Sys.time()
    ## Extract loop --------------------------------------------------
    ## empty df to store loop results
    extract_df <- data.frame()
    
    ## LOOP START 
    ## run through every month of time period
    for (i in 1:nrow(dates_df)) {
        
        ## Filter obs to mo/yr
        occ_filter <- occ %>% 
            filter(year == dates_df$year[i],
                   month == dates_df$mon_num[i])
        
        ## If filtered df has obs for that month, 
        ## then read in env data, vectorize occ data, and extract.
        if(nrow(occ_filter) > 0) {
            print(i)
            ## List of monthly raster files
            filesMonth <- list.files(
                path = pathMonth,
                ## only list those with matching yr/mo in name
                pattern = paste0(dates_df$year[i], dates_df$mon[i]),
                full.names = TRUE
            )
            ## List of quarterly raster files
            filesQuarter <- list.files(
                path = pathQuarter,
                pattern = paste0(dates_df$wy[i]),
                full.names = TRUE
            )
            ## Soil properties
            filesNatsgo <- list.files(
                path = pathNatsgo,
                pattern = paste0("natsgo_", ".+", ".tif"),
                full.names = TRUE
            )
            filesPolaris <- list.files(
                path = pathPolaris,
                full.names = TRUE
            )
            filesSalinity <- list.files(
                path = pathSalinity,
                pattern = 'CA_res',
                full.names = TRUE
            )
            month_rast = rast(filesMonth) %>% 
                crop(CA)
            quarter_rast = rast(filesQuarter) %>% 
                resample(month_rast)
            natsgo_rast = rast(filesNatsgo) %>% 
                resample(month_rast, method='mode')
            polaris_rast = rast(filesPolaris) %>% 
                resample(month_rast)
            salinity_rast = rast(filesSalinity) %>% 
                resample(month_rast, method='mode')
            names(salinity_rast) = "salinity"
            
            ## Stack all rasters
            env_stack <- terra::rast(list(month_rast, quarter_rast, natsgo_rast, polaris_rast, salinity_rast))
            
            ## vectorize and reproject to env data crs
            occ_vect <- occ_filter %>%
                terra::vect(geom = c(lon, lat), crs = "EPSG:4326") %>%
                terra::project(crs(env_stack))
            
            ## extract and tidy df
            sppExtract <- terra::extract(env_stack, occ_vect, method = "simple", ID = FALSE) %>% 
                ## only keep first 3 chars of monthly columns
                ## (e.g. "cwd2021jan" becomes "cwd")
                rename_with(.fn= ~substr(., 1, 3), 
                            ## only for columns of monthly variables
                            .cols = 1:length(filesMonth)) %>% 
                ## rename quarterly rasts; replace wy with _
                ## (e.g. "ppt2021winter_mean" becomes "ppt_winter_mean")
                rename_with(.fn= ~gsub(dates_df$wy[i], "_", x=.), 
                            ## only for columns of quarterly variables
                            .cols = (length(filesMonth)+1):(ncol(.)-length(filesPolaris))) %>% 
                ## merge occ data w/extract data
                cbind(occ_filter, .) %>% 
                rename(
                    bd_POLARIS=bd_0_200_mean,
                    clay_POLARIS=clay_0_200_mean,
                    om_POLARIS=om_0_200_mean,
                    ph_POLARIS=ph_0_200_mean,
                    sand_POLARIS=sand_0_200_mean,
                    silt_POLARIS=silt_0_200_mean
                )
            
            ## append results to df
            extract_df <- rbind(extract_df, sppExtract)
            
            ## If no obs for a yr/mo, skip  
        } else {
            next
        } ### END if/else statement
        
        
    } ### END LOOP
    Sys.time()-t0
    
    return(extract_df)
    
} ### end fxn

trainingTestingFunc = function(
        spp, occ.filename, back.filename,
        training.filename, testing.filename){
    occ = read_csv(occ.filename) %>% 
        select(-source) %>% 
        mutate(presence = 1)
    back = read_csv(back.filename) %>% 
        mutate(
            date = my(paste0(month, "_", year)),
            id = paste0('back_', row_number()),
            presence = 0
        ) %>% 
        select(any_of(names(occ)))
    
    combined = occ %>% 
        rbind(back)
    
    #Create training and testing sets
    
    set.seed(123)
    split = initial_split(combined, prop = 0.7, strata = "presence")
    training = training(split)
    testing = testing(split)
    write_csv(training, training.filename)
    write_csv(testing, testing.filename)
}
