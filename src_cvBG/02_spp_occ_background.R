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


# Background points

##Generate background points without month, which will be assigned later
source('src_cvBG/util/generate_backOcc.R')
dir.create('data_cvBG/2_background', recursive=T)
### 5km exclusion buffer to reduce false negatives
exclusionBuffer = 5000
### 20km around convex hull and CV to expand environmental state space and to inform model of unknown conditions for extrapolation (https://www.sciencedirect.com/science/article/pii/S030438001500215X)
outerBuffer = 20000
CA = tigris::states() %>%
    vect() %>%
    project("EPSG:3310") %>%
    filter(STUSPS=="CA")
CV = vect('data/central_valley/ds2632.gdb/')
### reference raster
bcmRaster = rast('data/0_env/bcm/bcmv8_historic/2000_2023_monthly/aet1999dec.tif')
referenceRaster = rast("data/0_env/natsgo/rasters/natsgo_drclass_270m_CA_2023.tif") %>%
    crop(CA, mask=T) %>%
    project('epsg:4326')
cvReferenceRaster = rast("data/0_env/natsgo/rasters/natsgo_drclass_270m_CA_2023.tif") %>%
    crop(CV, mask=T) %>%
    project('epsg:4326')
cvReferenceRaster = ifel(is.na(cvReferenceRaster), NA, 1)
referenceRaster = ifel(is.na(referenceRaster), NA, 1)
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
        backOcc_pts <- backOcc(sppOcc=sppOcc, referenceRaster=referenceRaster, cvReferenceRaster=cvReferenceRaster, exclusionBuffer=exclusionBuffer, outerBuffer=outerBuffer)

        ## Save
        write_csv(backOcc_pts, paste0("data_cvBG/2_background/back_",
                                      names[i],
                                      "_5km_lowFilter.csv"))
    }
)


## Read in fxn and set parameters --------------------------
#Use this for agg and comprehensive monthly model
source("src_cvBG/util/generateMonths_backOcc_evenMonth.R")


## Generate backOccs for each spp in list -------------------
purrr::map(
    .x = 1:length(names),
    .f = function(i) {
        set.seed(123)
        print(i)

        even_backOcc_pts.filename = paste0("data_cvBG/2_background/back_", names[i], "_5km_lowFilter_even.csv")
        even_backOcc_pts = generateEvenMonths(sp = names[i])
        write_csv(even_backOcc_pts, even_backOcc_pts.filename)
    }
)


# Extract environmental data


# Extract environmental data AKA SWD ----------------------------------------------

###Water year####
# Read in fxn
source('src_cvBG/util/env_extract_wy.R')
dir.create('data_cvBG/3_swd/wy', recursive=T)

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
        out.dir = paste0('data_cvBG/3_swd/wy')
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