#' Background (pseudo-absence) point  generation
#'
#' This function creates background points within a buffered range of supplied 
#' occurrence data (`sppOcc`). The temporal distribution of background points 
#' by year matches that of the occurrence data.
#' TJ: will be temporally distributing background points by month-year.
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

generateMatchingMonths <- function(sp) {
    background.filename = paste0("data/2_background/back_", sp, "_5km_lowFilter.csv")
    background = read_csv(background.filename)
    
    sppOcc.filename = paste0("data/1_occ/combined_spp_occ/", sp, "_lowFilter.csv")
    sppOcc = read_csv(sppOcc.filename) 
    
    ## Find occ ratio to apply to bkg pts
    month_wy_freq = sppOcc %>% 
        ## only keep occs in 2000-2022 wy
        mutate(wy = lfstat::water_year(.$date, origin = 10),
               ##convert factor to numeric
               wy = as.numeric(levels(wy))[wy]) %>% 
        dplyr::filter(wy >= 2000 & wy <= 2023) %>% 
        ## find counts and % counts for each water year and month
        count(wy, month) %>%
        group_by(wy) %>%
        mutate(freq = n/sum(n)) %>%
        ungroup()
    
    background_month_wy = background %>% 
        group_split(wy) %>% 
        map(
            function(df){
                this_wy = df %>% 
                    pull(wy) %>% 
                    unique()
                N = nrow(df)
                month_wy_amount = month_wy_freq %>% 
                    filter(wy == this_wy) %>% 
                    mutate(amount = ceiling(freq*N))
                months_to_sample = rep(month_wy_amount$month, month_wy_amount$amount)
                if(length(months_to_sample)==1){
                    df$month = months_to_sample
                } else{
                    set.seed(123)
                    month_vec = sample(
                        months_to_sample,
                        size = N
                    )
                    out = df %>% 
                        mutate(
                            month = month_vec,
                            ## we want calendar yr for later env extraction
                            year = case_when(month %in% c(10, 11, 12) ~(wy-1),
                                             .default = wy)
                        )
                }
                return(out)
            }
        ) %>% 
        bind_rows()
    return(background_month_wy)
}## End fxn
