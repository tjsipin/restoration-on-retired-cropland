#' @author Nick McManus, Will Dean, Madi Calbert, TJ Sipin
#' @date 10/27/2025
#' @description
#' This script generates and preps the environmental data used in species distribution modeling. 
#' Outputs are read in for extraction in the `spp_occ_background.Rmd`.
#' Aggregation.

library(tidyverse)    ## always
library(here)         ## reading/writing data
library(purrr)        ## run iterative fxns faster
library(sf)           ## vector objects
library(terra)        ## Better/faster GIS
library(raster)       ## GIS format required by Dismo
library(tidyterra)

# BCMv8 Data

### Convert ASC to TIF
#BCMv8 data can be directly downloaded from USGS at: https://www.sciencebase.gov/catalog/item/5f29c62d82cef313ed9edb39

## Read in fxn
source('src_lowFilter_v2/util/asc_to_tif.R')

## Assign file path (selects all .asc files in directory)
filepath = 'data/0_env/bcm/bcmv8_historic/2000_2023_monthly//'

## Run fxn
asc_to_tif(filepath, remove = TRUE)


### Quarterly rasters

#### Avg by water year
## Read in fxn
source("src_lowFilter_v2//util/quarterly_rast.R")

## Define fxn variables
pathIn = "data/0_env/bcm/bcmv8_historic/2000_2023_monthly//"
pathOut = "data/0_env/bcm/bcmv8_historic/quarterly_avgs//"
dir.create(pathOut)
startYear = 2000
endYear = 2023

## Run for winter ppt
quarter_rast(var="ppt", quarter="winter", method="sum",
            startYear, endYear,
            pathIn, pathOut)

## Run for summer tmx
quarter_rast(var="tmx", quarter="summer", method = "mean",
            startYear, endYear,
            pathIn, pathOut)

#### Avg for entire time period
## Avg of winter precip ----------------------------------------------------
dir.create('data/0_env/bcm/bcmv8_historic/monthly_avgs')

pptFiles <- list.files(
  ## read in the qtr avgs created in previous code chunk
  path = "data/0_env/bcm/bcmv8_historic/quarterly_avgs//",
  ## only select winter ppt
  pattern = paste0("ppt", ".+", "winter_sum"),
  full.names = TRUE
)

## read in files as raster "stack"
pptStack <- terra::rast(c(pptFiles))
## find mean of all rasts
pptStack_avg <- terra::app(pptStack, fun = 'mean')
## export mean rast
writeRaster(pptStack_avg,
            "data/0_env/bcm/bcmv8_historic/monthly_avgs/ppt_winter_avg.tif", overwrite=TRUE)


## Avg of summer tmx --------------------------------------------------------
tmxFiles <- list.files(
  path = "data/0_env/bcm/bcmv8_historic/quarterly_avgs//",
  ## only select tmx files
  pattern = "tmx",
  full.names = TRUE
)

## read in as stack, find mean, then export
tmxStack <- terra::rast(c(tmxFiles))
tmxStack_avg <- terra::app(tmxStack, fun = "mean")
writeRaster(tmxStack_avg,
            "data/0_env/bcm/bcmv8_historic/monthly_avgs/tmx_summer_avg.tif", overwrite=TRUE)


### Avg monthly rasters:

## df of variables to avg
vars <- data.frame(variable = rep(c("aet", "ppt", "tmn", "tmx", "cwd"), each = 12),
                   month = c('jan', 'feb', 'mar', 'apr', 'may', 'jun',
                             'jul', 'aug', 'sep', 'oct', 'nov', 'dec'))

## Fxn to avg -------------------------------------------------
var_avg <- function(variable, month, pathIn, pathOut) {
    ## Read in all files for that var/mo
    files <- list.files(path = pathIn, 
                        ## only list those with matching yr/mo in name
                        pattern = paste0(variable, ".+", month),
                        full=TRUE)
    ## Stack and avg
    stack <- terra::rast(c(files))
    stack_avg <- terra::app(stack, fun = 'mean')
    
    ## Save
    writeRaster(stack_avg, 
                paste0(pathOut, variable, "_", month, "_avg.tif"), 
                overwrite = TRUE)
}

## Run fxn w/pmap
pathIn = "data/0_env/bcm/bcmv8_historic/2000_2023_monthly//"
pathOut = "data/0_env/bcm/bcmv8_historic/monthly_avgs//"
purrr::pmap(.l=vars, .f=var_avg, pathIn, pathOut, .progress = TRUE)

vars <- data.frame(month = c('oct', 'nov', 'dec', 'jan', 'feb', 'mar', 
                             'apr', 'may', 'jun', 'jul', 'aug', 'sep'))

tdiff_avg <- function(month, pathIn, pathOut) {
    ## stack all tmx and tmn by month
    tmx.list <- list.files(path = pathIn,
                           pattern = paste0("tmx", ".+", month),
                           full.names = TRUE)
    tmn.list <- list.files(path = pathIn,
                           pattern = paste0("tmn", ".+", month),
                           full.names = TRUE)
    
    tmx <- terra::rast(c(tmx.list))
    tmn <- terra::rast(c(tmn.list))
    
    ## Find tdiff by month, then avg 
    tdiff <- tmx - tmn
    avg <- terra::app(tdiff, fun = "mean")
    
    ## Save
    writeRaster(avg, 
                paste0(pathOut, "tdiff_", month, "_avg.tif"), 
                overwrite = TRUE)
}

## Run fxn w/pmap
pathIn = "data/0_env/bcm/bcmv8_historic/2000_2023_monthly//"
pathOut = "data/0_env/bcm/bcmv8_historic/monthly_avgs//"
purrr::pmap(.l=vars, .f=tdiff_avg, pathIn, pathOut, .progress = TRUE)



### Future BCM data
#Finally, to generate future species distribution probability maps, we'll have to wrangle projected environmental data. 

### Function to reproject, resample, and save rasters ----------------------
bcm_future_reproj <- function(model, month, ref_r) {
  ## path to "raw" files
  pathIn <- paste0("data/0_env/bcm/bcm_future//",
                   model,
                   "/raw")
  ## path to save new files
  pathOut <- paste0("data/0_env/bcm/bcm_future//",
                    model,
                    "/resampled//")
  
  # Read in all the files for 30yr monthly data
  raw_files <- list.files(
    path = pathIn,
    pattern = paste0(month, ".+", ".tif$"),
    ## look through every folder in directory
    recursive = TRUE,
    full.names = TRUE
  )
  ## Both .rgb.tif and .tif files present, so keep only the .tif
  raw_files <- grep(raw_files,
                    pattern = ".rgb",
                    invert = TRUE,
                    value = TRUE)
  
  ## Loop through each file for that month and save
  purrr::map(.x=raw_files, function(x) {
    r <- rast(x)
    ## manually assign crs
    ## epsg:9001 cannot do terra::project() transformation
    crs(r) <- crs(ref_r)
    ## resample to ensure matching extent
    r_res <- resample(r, ref_r, method = "bilinear")
    
    ## only keep var from the long file name
    ## (needs to match var name from model)
    old_name <- names(r_res)
    var <- substr(x = old_name, start = 1, stop = 3)
    names(r_res) <- var
    
    ## Save
    writeRaster(r_res,
                paste0(pathOut, var, "_", month, "_2070_2099_", model, ".tif"),
                overwrite = TRUE)
  })##End map
  
  ## Generate TDIFF raster for month
  tmx <- rast(paste0(pathOut, "tmx_", month, "_2070_2099_", model, ".tif"))
  tmn <- rast(paste0(pathOut, "tmn_", month, "_2070_2099_", model, ".tif"))
  tdiff <- tmx - tmn
  names(tdiff) <- "tdiff"
  writeRaster(tdiff,
              paste0(pathOut, "tdiff_", month, "_2070_2099_", model, ".tif"),
              overwrite = TRUE)
  
}##END fxn


### Define fxn variables and iterate over all months -----------------------
## reference raster for crs and resampling
ref_r <- rast("data/0_env/bcm/bcmv8_historic/monthly_avgs/aet_apr_avg.tif")

## df of months
months <- data.frame(month = c("jan", "feb", "mar", "apr", "may", "jun",
                               "jul", "aug", "sep", "oct", "nov", "dec"))

## iterate fxn over months for specific model
purrr::pmap(.l = months, 
            .f = bcm_future_reproj, 
            model = "MIROC85",
            ref_r = ref_r,
            .progress = TRUE)



## Winter precip ----------------------------------------------------
pptFiles <- list.files(
    ## read in the data just generated
    path = "data/0_env/bcm/bcm_future/MIROC85/resampled//",
    ## only select winter ppt
    pattern = paste0("ppt", ".+", "dec|ppt", ".+", "jan|ppt", ".+", "feb"),
    full.names = TRUE
)

## read in files as raster "stack"
pptStack <- terra::rast(c(pptFiles))
## find mean of all rasts
pptStack_sum <- terra::app(pptStack, fun = 'sum')
## export mean rast
writeRaster(pptStack_sum,
            "data/0_env/bcm/bcm_future/MIROC85/resampled/ppt_winter_2070_2099_MIROC85.tif")

## Summer tmx --------------------------------------------------------
tmxFiles <- list.files(
    ## read in the data just generated
    path = "data/0_env/bcm/bcm_future/MIROC85/resampled/",
    ## only select winter ppt
    pattern = paste0("tmx", ".+", "jun|tmx", ".+", "jul|tmx", ".+", "aug"),
    full.names = TRUE
)

## read in files as raster "stack"
tmxStack <- terra::rast(c(tmxFiles))
## find mean of all rasts
tmxStack_avg <- terra::app(tmxStack, fun = 'mean')
## export mean rast
writeRaster(tmxStack_avg,
            "data/0_env/bcm/bcm_future/MIROC85/resampled/tmx_summer_2070_2099_MIROC85.tif")




# gNATSGO Data


### Aggregate soil data

## database pathway
gdb = "data/0_env/natsgo/gNATSGO_CA/gNATSGO_CA.gdb"

## Check the layers in the .gdb
# st_layers(gdb)

## Read in soil data using sf package to specify layer in database
## Comes in as df
horizon = st_read(gdb, layer = "chorizon")
component = st_read(gdb, layer = "component")
mapunit = st_read(gdb, layer = "mapunit")
## provides pre-aggregated data on drainage class,
## but not for other soil properties of interest
muagg = st_read(gdb, layer = "muaggatt")

## Read in and run fxn
source('src_lowFilter_v2/util/natsgo_agg.R')

natsgo_agg(horizon, component, mapunit, muagg,
           depth = 200, pathOut = "data/0_env/natsgo//")


### Rasterize soil data

## Read in soil data and raster
mu_r <- rast("data/0_env/natsgo/rasters/MapunitRaster_10m_CA_2023.tif")

## Change resolution to 270m
mu270_r <- aggregate(mu_r, fact=27, 
                     ## mapunits are categorical
                     fun="modal")
## Save
writeRaster(mu270_r, "data/0_env/natsgo/rasters/MapunitRaster_270m_CA_2023.tif")


#Now we'll create a new raster for each soil variable by reclassifying the mapunit raster and matching the extent of BCM data:
## Read in data ----------------------------------------------------------
## mapunit raster
mu270_r <- rast("data/0_env/natsgo/rasters/MapunitRaster_270m_CA_2023.tif")
## sample bcm raster
bcm_r<- rast("data/0_env/bcmv8/monthly_avgs/tmx_jan_avg.tif")
## soil data
soil_df <- read_csv("data/0_env/natsgo/horizon_200cm_CA.csv") %>% 
  ## turn drainage into factor, set levels
  mutate(drclass = factor(drclass, 
                             levels = c("Excessively drained",
                                        "Somewhat excessively drained",
                                        "Well drained",
                                        "Moderately well drained",
                                        "Somewhat poorly drained",
                                        "Poorly drained",
                                        "Very poorly drained"))) %>% 
  ## then reclass as numeric
  mutate(drclass = as.numeric(drclass))


## Create fxn to reclassify and match BCM raster extent -------------------
soil_rast <- function(var, type) {
  ## Select one variable at a time
  rcl <- soil_df %>% 
    dplyr::select(mukey, var) %>% 
    rbind(c(0, NA)) #make outside areas NA
  
  ## reclassify
  rcl_r <- classify(mu270_r, rcl)
  
  ## if continuous variable, resamp w/bilinear method
  if (type == "cont") {
    rcl_resamp_r <- rcl_r %>%
      project(y = crs(bcm_r), method = "bilinear") %>%
      resample(y = bcm_r, method = "bilinear")
  ## if categorical, use near method
  } else {
    rcl_resamp_r <- rcl_r %>%
      project(y = crs(bcm_r), method = "near") %>%
      resample(y = bcm_r, method = "near")
  }
  
  ## proper raster file name
  names(rcl_resamp_r) <- var
  
  ## Save
  writeRaster(rcl_resamp_r,
              paste0(
                "data/0_env/natsgo/rasters/natsgo_",
                var,
                "_270m_CA_2023.tif"
              ),
              overwrite = TRUE)
} ##End fxn

soil_rast(var = "drclass", type = "cat")

## Iterate fxn over all soil variables ------------------------------------
## var names must match columns in `soil_df`
soil_vars = data.frame(var = c("cec", "ph", "om", "drclass"),
                       type = c("cont", "cont", "cont", "cat"))
## run w/pmap
purrr::pmap(.l = soil_vars, 
            .f = soil_rast, 
            .progress = TRUE)

##Salinity from 
##https://data.isric.org/geonetwork/srv/eng/catalog.search;jsessionid=9251411A3E92851C12FAA0C06EB6745F#/metadata/c59d0162-a258-4210-af80-777d7929c512
#Read in salinity raster that is within CA
sal <- rast('data/0_env/salinity/salMap2016-0000000000-0000000000.tif')
# CA = tigris::states() %>%
#     filter(STUSPS=="CA") %>%
#     vect()
## reference raster for crs and resampling
ref_r <- rast("data/0_env/bcm/bcmv8_historic/monthly_avgs/aet_apr_avg.tif") %>% 
    project('epsg:4326')

# ref_r_sal = ref_r %>% 
#     project(sal)
t0 = Sys.time()
sal_CA = sal %>% 
    crop(ref_r, mask = T) %>% 
    project("EPSG:3310", method = 'mode')
Sys.time() - t0
t0 = Sys.time()
x = sal_CA
res(x) = 270
sal_CA_res = sal_CA %>% 
    resample(x, method = 'mode')
Sys.time() - t0
writeRaster(sal, 'data/0_env/salinity/salinity_full.tif')
writeRaster(sal_CA_res, 'data/0_env/salinity/salinity_full_CA_res.tif')



# To align BCM and Polaris Data:

# polaris_silt <- rast(here("data/0_env/polaris/rasters/polaris_silt_0200_270m_lowFilter.tif"))
# polaris_ph <- rast(here("data/0_env/polaris/rasters/polaris_ph_0200_270m_lowFilter.tif"))
# polaris_om <- rast(here("data/0_env/polaris/rasters/polaris_om_0200_270m_lowFilter.tif"))
# polaris_sand <- rast(here("data/0_env/polaris/rasters/polaris_sand_0200_270m_lowFilter.tif"))
# polaris_clay <-rast(here("data/0_env/polaris/rasters/polaris_clay_0200_270m_lowFilter.tif"))
# polaris_bd <-rast(here("data/0_env/polaris/rasters/polaris_bd_0200_270m_lowFilter.tif"))
# 
# #Read in sample BCM raster to match extent and resolution
# bcm_r <- rast(here("data/0_env/bcm/bcmv8_historic/monthly_avgs/tmx_jan_avg.tif"))
# 
# #Reproject and resample POLARIS rasters to match BCM resolution (270m)
# polaris_silt_resampled <- resample(polaris_silt, bcm_r, method = "bilinear")
# polaris_ph_resampled <- resample(polaris_ph, bcm_r, method = "bilinear")
# polaris_om_resampled <- resample(polaris_om, bcm_r, method = "bilinear")
# polaris_sand_resampled <- resample(polaris_sand, bcm_r, method = "bilinear")
# polaris_clay_resampled <- resample(polaris_clay, bcm_r, method = "bilinear")
# polaris_bd_resampled <- resample(polaris_bd, bcm_r, method = "bilinear")
# 
# # Create a function to update the POLARIS rasters with the associated BCM data
# update_raster <- function(polaris_raster, bcm_r, var_name) {
#   # Resample the POLARIS raster to match BCM resolution and extent
#   resampled_raster <- resample(polaris_raster, bcm_r, method = "bilinear")
#   
#   # Save the updated raster with the same name but prefixed with 'updated_'
#   writeRaster(resampled_raster, 
#               here(paste0("data/0_env/polaris/rasters/updated_", var_name, "_270m_lowFilter.tif")), 
#               overwrite = TRUE)
# }
# 
# # Update and save each POLARIS raster
# update_raster(polaris_silt, bcm_r, "polaris_silt")
# update_raster(polaris_ph, bcm_r, "polaris_ph")
# update_raster(polaris_om, bcm_r, "polaris_om")
# update_raster(polaris_sand, bcm_r, "polaris_sand")
# update_raster(polaris_clay, bcm_r, "polaris_clay")
# update_raster(polaris_bd, bcm_r, "polaris_bd")
