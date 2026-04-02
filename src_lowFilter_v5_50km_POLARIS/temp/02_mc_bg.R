#' @author Nick McManus, Will Dean, Madi Calbert, TJ Sipin
#' @date 10/27/2025
#' @description
#' This script generates and preps occurrence and background data for running the SDMs.
#' First, occurrence records from GBIF and CalFlora are filtered, merged, and spatially thinned.
#' Second, a set of random background points are created for each species.
#' Finally, environmental data are extracted for all occurrence and background points.


library(tidyverse)            ## always
library(here)                 ## reading/writing data
library(rgbif)                ## download GBIF data
library(CoordinateCleaner)    ## filter/clean GBIF data
library(terra)                ## rast pkg for quicker reprojecting
library(raster)               ## rast format that plays w/dismo
library(enmSdmX)              ## spatially thinning data
library(dismo)                ## generating background points
library(lfstat)               ## water year fxn
library(tidymodels)
library(tidyterra)
library(sf)
tidymodels_prefer()

names <- c(
    "a_polycarpa",
    "p_arborea",
    "c_pungens",
    "l_pentachaeta",
    "p_ciliata",
    "a_menziesii",
    "a_intermedia",
    "c_lasiophyllus",
    'l_californica',
    'l_gracilis'
)

#####Monthly#####
## Read in fxn
spp_occ_back = tibble(
    spp = rep(names, 2),
    occ_back = c(rep('occ', length(names)), rep('back', length(names)))
)

## Fxn variables
startYear = 2000
endYear = 2023
pathMonth = "data/0_env/bcm/bcmv8_historic/2000_2023_monthly/"
pathQuarter = "data/0_env/bcm/bcmv8_historic/quarterly_avgs/"
pathNatsgo = "data/0_env/natsgo/rasters/"
pathPolaris = "data/0_env/polaris/rasters/"
pathSalinity = "data/0_env/salinity/"

dir.create('data_v5_50km/3_swd/monthly/selective', recursive=T)
source('src_lowFilter_v5_50km/util/env_extract_monthly.R')

##Extract for even months background set
dir.create('data_v5_50km/3_swd/monthly/comprehensive', recursive=T)
set.seed(123)
map(
    1:nrow(spp_occ_back),
    function(i){
        spp = spp_occ_back$spp[i]
        occ_back = spp_occ_back$occ_back[i]
        out.filename = paste0('data_v5_50km/3_swd/monthly/comprehensive/swd_', spp, '_', occ_back, '_soil200cm_lowFilter_monthly_comprehensive.csv')
        if(file.exists(out.filename)) return(NULL)
        print(i)
        if(occ_back=='occ'){
            points.filename = paste0('data/1_occ/combined_spp_occ/', spp, '_lowFilter.csv')
        } else{
            points.filename = paste0("data_v5_50km/2_background/back_", spp, "_5km_lowFilter_even.csv")
        }
        points = read_csv(points.filename)
        
        print(spp)
        out = extractEnvMonthly(
            startYear = startYear,
            endYear = endYear,
            pathMonth = pathMonth,
            pathQuarter = pathQuarter,
            pathNatsgo = pathNatsgo,
            pathPolaris = pathPolaris,
            pathSalinity = pathSalinity,
            occ = points,
            lon = 'lon',
            lat = 'lat'
        )
        
        write_csv(out, out.filename)
        
    },
    .progress = T
)

map(
    names,
    function(sp){
        print(sp)
        out.dir = paste0('data_v5_50km/3_swd/monthly/comprehensive')
        #Get occ and back extracted sets
        occ.filename = paste0(out.dir, '/swd_', sp, '_occ_soil200cm_lowFilter_monthly_comprehensive.csv')
        back.filename = paste0(out.dir, '/swd_', sp, '_back_soil200cm_lowFilter_monthly_comprehensive.csv')
        
        training.filename = paste0(out.dir, '/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')
        testing.filename = paste0(out.dir, '/testing_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')
        trainingTestingFunc(
            sp=sp,
            occ.filename=occ.filename,
            back.filename=back.filename,
            training.filename=training.filename,
            testing.filename=testing.filename
        )
    }
)