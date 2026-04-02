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

backOcc <- function(sppOcc, referenceRaster, cvReferenceRaster, exclusionBuffer, outerBuffer, back_n = 10000,
                    lon = "lon", lat = "lat", crs = "epsg:4326") {
    
    sppOcc = sppOcc %>% 
        select(lon, lat, month, year, everything())
    
    sppOcc_vect = terra::vect(sppOcc, 
                              geom = c(lon, lat),
                              crs = crs)
    sppExclusion = sppOcc_vect %>% 
        terra::buffer(width=exclusionBuffer)
    ## Vectorize sppOcc, find convex hull of pts, then add buffer
    sppZone <- terra::convHull(sppOcc_vect) %>% 
        terra::buffer(., width = outerBuffer)
    cvZone = cvReferenceRaster %>% 
        as.polygons() %>% 
        terra::buffer(width=outerBuffer)
    
    ## Keep only area of reference raster w/in sppZone
    convHullBuffer <- referenceRaster %>% 
        ## match crs of raster to spp polygon
        terra::project(y=crs(sppZone)) %>% 
        ## crop and mask to make all cells outside polygon NA
        # terra::crop(y=sppZone, mask = TRUE)
        terra::mask(sppZone)
    cvBuffer = referenceRaster %>% 
        # terra::crop(y=cvZone, mask = TRUE)
        terra::mask(cvZone)
    convHullBufferExclusion = convHullBuffer %>% 
        terra::mask(sppExclusion, inverse = T)
    cvBufferExclusion = cvBuffer %>% 
        terra::mask(sppExclusion, inverse = T)
    combineSampleRaster = rast(list(
        convHull=convHullBufferExclusion, cv=cvBufferExclusion
    ))
    combineSampleRaster$convHull = ifel(is.na(combineSampleRaster$convHull), 0, 1)
    combineSampleRaster$cv = ifel(is.na(combineSampleRaster$cv), 0, 1)
    combineSampleRaster$combine = combineSampleRaster$convHull + combineSampleRaster$cv
    
        
    ## Raster pkg format required for dismo 
    r <- raster::raster(combineSampleRaster$combine)
    # r = sppZone %>% raster()
    
    ## Find occ ratio to apply to bkg pts
    occ_count <- sppOcc %>% 
        ## only keep occs in 2000-2022 wy
        mutate(wy = lfstat::water_year(.$date, origin = 10),
               ##convert factor to numeric
               wy = as.numeric(levels(wy))[wy]) %>% 
        dplyr::filter(wy >= 2000 & wy <= 2023) %>% 
        ## find counts and % counts for each water year
        group_by(wy) %>% 
        summarize(occ_n = n()) %>% 
        mutate(occ_n_pct = occ_n/sum(.$occ_n))
    
    ## list of water years
    w_years <- rep(2000:2023, each = 1)
    
    ## Start empty df for loop below
    ## (will remove these NA values at end)
    backOcc_total <- data.frame("lon" = NA, "lat" = NA, "month"=NA, "year" = NA, "wy" = NA)
    
    ## Loop through every wy to generate bkg pts ----------------------------
    backOcc_total = map(
        w_years, 
        .f = function(yyyy){
            print(yyyy)
            occ_filt <- occ_count %>% 
                dplyr::filter(wy == yyyy)
            
            ## If there are occs in the wy, 
            ## then generate number of bkg pts based on occ ratio
            if(nrow(occ_filt) != 0) {
                ## number pts to be generated rounded up to be equal per month
                roundUp <- function(x) {12*ceiling(x/12)}
                x <- ceiling(occ_filt$occ_n_pct * back_n)
                n <- roundUp(x)
                
                ## Generate random pts
                backOcc <- as.data.frame(
                    dismo::randomPoints(mask = r, 
                                        n = n,  
                                        prob = T,
                                        ## ensures bkg pt not in same cell as spp occ
                                        p = sppOcc_vect %>% 
                                            st_as_sf() %>% 
                                            as_Spatial(),
                                        excludep = TRUE)
                ) %>% 
                    rename(lon=x, lat=y) %>% 
                    mutate(wy = yyyy)
                
                return(backOcc)
                
                ## if no occs in a wy, skip to the next wy in list
            } else {
                return(NULL)
            }##End if/else
        }
    ) %>% 
        bind_rows()
    return(backOcc_total)
}## End fxn
