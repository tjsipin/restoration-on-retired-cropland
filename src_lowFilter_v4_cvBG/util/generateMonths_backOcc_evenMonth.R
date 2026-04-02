#' Background (pseudo-absence) point  generation
#'
#' This function creates background points within a buffered range of supplied 
#' occurrence data (`sppOcc`). The temporal distribution of background points 
#' by year matches that of the occurrence data.
#' 
#' @author Nick McManus
#' @param sppOcc data frame. Contains species occurrence points, with x and y as separate variables and a dates variable as a "Date" class
#' @param back_n integer. Target number of background points to be generated (default = 10000). Function will likely return slightly larger number
#' @param raster SpatRaster (terra raster format). Reference raster for the spatial resolution of background points. 
#' As written, no more than one point will be generated per raster cell. Recommended that raster matches the spatial resolution of occurrence data.
#' @param buffer integer. Distance (in meters) from occurrence points to generate background points.
#' @param lon character. Variable name for longitude values in sppOcc df (default = "lon")
#' @param lat character. Variable name for the latitude values in sppOcc df (default = "lat")
#' @param crs character. The coordinate reference system of the sppOcc occurrence data (default = "WGS84")
#' @return data frame with coordinates, month, and year of background points
set.seed(123)

generateEvenMonths <- function(sp) {
    background.filename = paste0("data_cvBG/2_background/back_", sp, "_5km_lowFilter.csv")
    background = read_csv(background.filename)
    
    sppOcc.filename = paste0("data/1_occ/combined_spp_occ/", sp, "_lowFilter.csv")
    sppOcc = read_csv(sppOcc.filename) 
    
    N = nrow(background)
    
    ## Find occ ratio to apply to bkg pts
    wy_freq = sppOcc %>% 
        ## only keep occs in 2000-2022 wy
        mutate(wy = lfstat::water_year(.$date, origin = 10),
               ##convert factor to numeric
               wy = as.numeric(levels(wy))[wy]) %>% 
        dplyr::filter(wy >= 2000 & wy <= 2023) %>% 
        ## find counts and % counts for each water year and month
        count(wy) %>% 
        mutate(freq = n/sum(n)) %>% 
        mutate(amount = ceiling(freq*N))
    
    background_month_wy = background %>% 
        group_split(wy) %>% 
        map(
            function(df){
                
                this_wy = df %>% 
                    pull(wy) %>% 
                    unique()
                
                occ_filt <- wy_freq %>% 
                    dplyr::filter(wy == this_wy)
                
                ## If there are occs in the wy, 
                ## then generate number of bkg pts based on occ ratio
                if(nrow(occ_filt) != 0) {
                    n = nrow(df)
                    
                    ## Assign mo/yr to new bkg pts
                    dates <- data.frame(month = rep(c(10:12, 1:9), 
                                                    each = (n/12)))
                    backOcc_dates <- cbind(df, dates) %>% 
                        ## we want calendar yr for later env extraction
                        mutate(year = case_when(month %in% c(10, 11, 12) ~(wy-1),
                                                .default = wy))
                }
                
                return(backOcc_dates)
            }
        ) %>% 
        bind_rows()
    return(background_month_wy)
    
}## End fxn
