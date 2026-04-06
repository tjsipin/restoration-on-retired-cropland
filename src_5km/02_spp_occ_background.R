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

#Call in occurrences from data_MAIN

# Background points

##Generate background points without month, which will be assigned later
source('src_5km/util/generate_backOcc.R')
dir.create('data_5km/2_background')
### 5km buffer
buffer = 5000
CA = tigris::states() %>%
    vect() %>%
    project("EPSG:3310") %>%
    filter(STUSPS=="CA")
### reference raster
raster = rast("data/0_env/natsgo/rasters/natsgo_ph_270m_CA_2023.tif") %>%
    crop(CA, mask=T)
purrr::map(
    .x = 1:length(names),
    .f = function(i) {
        print(i)
        set.seed(123)
        ## read in spp occurrence points
        sppOcc = read_csv(paste0("data/1_occ/combined_spp_occ//",
                                 names[i],
                                 "_lowFilter.csv"))
        
        ## Generate pts w/fxn
        backOcc_pts <- backOcc(sppOcc, raster=raster, buffer=buffer)
        
        ## Save
        write_csv(backOcc_pts, paste0("data_5km/2_background/back_",
                                      names[i],
                                      "_5km_lowFilter.csv"))
    }
)


## Read in fxn and set parameters --------------------------
#Use this for selective monthly model
source("src_5km/util/generateMonths_backOcc_matchMonth.R")
#Use this for agg and comprehensive monthly model
source("src_5km/util/generateMonths_backOcc_evenMonth.R")


## Generate backOccs for each spp in list -------------------
purrr::map(
    .x = 1:length(names),
    .f = function(i) {
        set.seed(123)
        print(i)
        matching_backOcc_pts.filename = paste0("data_5km/2_background/back_", names[i], "_5km_lowFilter_match.csv")
        matching_backOcc_pts = generateMatchingMonths(sp = names[i])
        write_csv(matching_backOcc_pts, matching_backOcc_pts.filename)
        
        even_backOcc_pts.filename = paste0("data_5km/2_background/back_", names[i], "_5km_lowFilter_even.csv")
        even_backOcc_pts = generateEvenMonths(sp = names[i])
        write_csv(even_backOcc_pts, even_backOcc_pts.filename)
    }
)


# Extract environmental data


# Extract environmental data AKA SWD ----------------------------------------------

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

dir.create('data_5km/3_swd/monthly/selective', recursive=T)
source('src_5km//util/env_extract_monthly.R')

##Extract for matched months background set
map(
    1:nrow(spp_occ_back),
    function(i){
        set.seed(123)
        spp = spp_occ_back$spp[i]
        occ_back = spp_occ_back$occ_back[i]
        out.filename = paste0('data_5km/3_swd/monthly/selective/swd_', spp, '_', occ_back, '_soil200cm_lowFilter_monthly_selective.csv')
        
        # if(file.exists(out.filename)) return(NULL)
        if(occ_back=='occ'){
            points.filename = paste0('data/1_occ/combined_spp_occ/', spp, '_lowFilter.csv')
        } else{
            points.filename = paste0("data_5km/2_background/back_", spp, "_5km_lowFilter_match.csv")
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
        out.dir = paste0('data_5km/3_swd/monthly/selective')
        #Get occ and back extracted sets
        occ.filename = paste0(out.dir, '/swd_', sp, '_occ_soil200cm_lowFilter_monthly_selective.csv')
        back.filename = paste0(out.dir, '/swd_', sp, '_back_soil200cm_lowFilter_monthly_selective.csv')
        
        training.filename = paste0(out.dir, '/training_', sp, '_soil200cm_lowFilter_monthly_selective.csv')
        testing.filename = paste0(out.dir, '/testing_', sp, '_soil200cm_lowFilter_monthly_selective.csv')
        trainingTestingFunc(
            sp=sp,
            occ.filename=occ.filename,
            back.filename=back.filename,
            training.filename=training.filename,
            testing.filename=testing.filename
        )
    }
)

##Extract for even months background set
dir.create('data_5km/3_swd/monthly/', recursive=T)
set.seed(123)
map(
    1:nrow(spp_occ_back),
    function(i){
        spp = spp_occ_back$spp[i]
        occ_back = spp_occ_back$occ_back[i]
        out.filename = paste0('data_5km/3_swd/monthly/swd_', spp, '_', occ_back, '_soil200cm_lowFilter_monthly.csv')
        # if(file.exists(out.filename)) return(NULL)
        if(occ_back=='occ'){
            points.filename = paste0('data/1_occ/combined_spp_occ/', spp, '_lowFilter.csv')
        } else{
            points.filename = paste0("data_5km/2_background/back_", spp, "_5km_lowFilter_even.csv")
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
        out.dir = paste0('data_5km/3_swd/monthly/comprehensive')
        #Get occ and back extracted sets
        occ.filename = paste0(out.dir, '/swd_', sp, '_occ_soil200cm_lowFilter_monthly.csv')
        back.filename = paste0(out.dir, '/swd_', sp, '_back_soil200cm_lowFilter_monthly.csv')
        
        training.filename = paste0(out.dir, '/training_', sp, '_soil200cm_lowFilter_monthly.csv')
        testing.filename = paste0(out.dir, '/testing_', sp, '_soil200cm_lowFilter_monthly.csv')
        trainingTestingFunc(
            sp=sp,
            occ.filename=occ.filename,
            back.filename=back.filename,
            training.filename=training.filename,
            testing.filename=testing.filename
        )
    }
)

###Water year####
# Read in fxn
source('src_5km//util/env_extract_wy.R')
dir.create('data_5km/3_swd/wy', recursive=T)

set.seed(123)
map(
    names,
    function(sp){
        print(sp)
        extractEnvWY(sp = sp)
    },
    .progress = T
)

set.seed(123)
map(
    names,
    function(sp){
        out.dir = paste0('data_5km/3_swd/wy')
        #Get occ and back extracted sets
        in.filename = paste0(out.dir, '/swd_', sp, '_soil200cm_lowFilter_wy.csv')
        
        training.filename = paste0(out.dir, '/training_', sp, '_soil200cm_lowFilter_wy.csv')
        testing.filename = paste0(out.dir, '/testing_', sp, '_soil200cm_lowFilter_wy.csv')
        
        trainingTestingFunc(
            sp = sp,
            in.filename = in.filename,
            training.filename = training.filename,
            testing.filename = testing.filename
        )
    }
)

####Agg####
## Read in fxn
source('src_5km//util/env_extract_agg.R')
dir.create('data_5km/3_swd/agg', recursive=T)
spp_relevant_months = tibble(
    sp = c(
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
    ),
    mons = c(
        list(2:8), #a_polycarpa
        list(c(10:12, 1:7)), #p_arborea
        list(2:9), #c_pungens
        list(c(10:12, 1:5)), #l_pentachaeta
        list(1:5), #p_ciliata
        list(c(11:12, 1:5)), #a_menziesii
        list(c(11:12, 1:6)), #a_intermedia
        list(c(11:12, 1:6)), #c_lasiophyllus
        list(c(10:6)), #l_californica
        list(c(10:6)) #l_gracilis
    )
)

map(
    1:nrow(spp_relevant_months),
    function(i){
        sp = spp_relevant_months$sp[i]
        relevant_months = spp_relevant_months$mons[i] %>% unlist()
        print(sp)
        extractEnvAgg(sp = sp, relevant_months = relevant_months)
    },
    .progress = T
)

map(
    names,
    function(sp){
        out.dir = paste0('data_5km/3_swd/agg')
        #Get occ and back extracted sets
        in.filename = paste0(out.dir, '/swd_', sp, '_soil200cm_lowFilter_agg.csv')
        
        training.filename = paste0(out.dir, '/training_', sp, '_soil200cm_lowFilter_agg.csv')
        testing.filename = paste0(out.dir, '/testing_', sp, '_soil200cm_lowFilter_agg.csv')
        
        trainingTestingFunc(
            sp = sp,
            in.filename = in.filename,
            training.filename = training.filename,
            testing.filename = testing.filename
        )
    }
)

