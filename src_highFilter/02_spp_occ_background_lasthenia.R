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

# Occurrence Records

## Download/import

### Phenology dictionary
repro_values <- read_csv("data/1_occ/gbif_codes.csv")

### CalFlora
calflora_path = "data/1_occ/calflora_downloads_030626/"

## List of "accepted" sources
sources <- paste(c("BLM",
                   "Bureau",
                   "USDA",
                   "DFW",
                   "USGS",
                   "Nature Conservancy",
                   "TNC",
                   "CNPS",
                   "Taylor",
                   "Hrusa"), collapse = "|")

## Loop through each spp to read in, filter, and export
for (i in 1:length(names)) {
    print(names[i])
    ## Read in data
    input.filename = list.files(
        paste0(calflora_path),
        pattern = paste0(names[i], '_calflora_download'),
        full.names = T
    )
    df <- read_csv(input.filename) %>%
        janitor::clean_names()
    print(df$date %>% head())

    ## filter and keep coords
    df_filter <- df %>%
        ## "LOW" filter criteria
        filter(location_quality %in% c("high", "medium"),
               accuracy_square_meters <= 72900,
               date >= "1999-10-01") %>% 
        filter((str_detect(source, sources) | dataset == "cch2"))
        # filter(str_detect(source, sources))

    df_repro_filter = df_filter %>%
        filter(phenology_code %in% c("F", "R") | is.na(phenology_code)) %>%
        mutate(phenology_code = as.character(phenology_code))

    ## export
    write.csv(df_filter, paste0(calflora_path, names[i], "_calflora_highFilter_noReproFilter.csv"))
    write.csv(df_repro_filter, paste0(calflora_path, names[i], "_calflora_highFilter.csv"))
}


### GBIF
## Pull taxon keys from list of spp
taxon_keys <- name_backbone_checklist(c(
    ## shrubs
    "Atriplex polycarpa",
    "Cleomella arborea var. globosa (Coville) J.C.Hall & Roalson",
    ## forbs
    "Centromadia pungens subsp. pungens",
    "Layia pentachaeta subsp. albida D.D.Keck",
    "Phacelia ciliata Benth.",
    "Amsinckia menziesii (Lehm.) A.Nelson & J.F.Macbr.",
    "Amsinckia intermedia Fisch. & C.A.Mey",
    "Caulanthus lasiophyllus (Hook. & Arn.) Payson",
    "Lasthenia californica DC. ex Lindl",
    "Lasthenia gracilis (DC.) Greene"
)) %>%
    pull(usageKey)


dir.create('data/1_occ/gbif_downloads_122625/extracts', recursive = T)
## Loop through spp list and save each
for (i in 1:length(taxon_keys)) {
    #CAUTION: Do not run as this may overwrite our gbif downloads!
    ### download filtered data for A. polycarpa
    gbif_download <- occ_download(
        pred_in("taxonKey", taxon_keys[i]),
        ## remove geospatial issues
        pred("hasGeospatialIssue", FALSE),
        ## ensure coords
        pred("hasCoordinate", TRUE),
        ## remove "absent" occurences
        pred("occurrenceStatus", "PRESENT"),
        ## within US
        pred("country", "US"),
        ## within CA
        pred("stateProvince", "California"),
        ## output as CSV
        format = "DWCA",####################################################################################################
        ## enter your GBIF credentials below:
        user = "username",
        pwd = "pw",
        email = "email"
    )

    ### check on download status
    occ_download_wait(gbif_download)

    ### import GBIF data into env and filter using CoordinateCleaner pkg
    species <- occ_download_get(
        gbif_download,
        path = "data/1_occ/gbif_downloads_122625/extracts/",
        overwrite = F
    ) %>%
        occ_download_import() %>%
        ## set lowercase column names to work with CC
        setNames(tolower(names(.))) %>%
        ## filter out duplicate points
        distinct(decimallongitude, decimallatitude,
                 specieskey, datasetkey, .keep_all = TRUE) %>%
        ## filter known uncertainty below 270 and keep NAs
        filter(coordinateuncertaintyinmeters < 270 |
                   is.na(coordinateuncertaintyinmeters)) %>%
        ## known inaccurate default values
        filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>%
        ## remove herbaria/zoo locations
        cc_inst(lon = "decimallongitude", lat = "decimallatitude",
                buffer = 270, value = "clean", verbose = TRUE) %>%
        ## remove ocean values
        cc_sea(lon = "decimallongitude", lat = "decimallatitude") %>%
        ## remove points before 2000 wy
        filter(eventdate >= "1999-10-01") %>%
        mutate(reproductivecondition_clean = tolower(str_trim(reproductivecondition))) %>%
        left_join(repro_values)

    species_repro_filter = species %>%
        filter(phenology_code %in% c("F", "R") | is.na(phenology_code))

    ## export file
    dir.create(paste0("data/1_occ/gbif_downloads_122625/",names[i]))
    write_csv(species, paste0("data/1_occ/gbif_downloads_122625/",names[i],"/", names[i], "_gbif_lowFilter_noReproFilter.csv"))
    write_csv(species_repro_filter, paste0("data/1_occ/gbif_downloads_122625/", names[i],"/", names[i], "_gbif_lowFilter.csv"))

} ## END LOOP

#See lifestage counts per species
map(
    1:length(names),
    function(i){
        read_csv(paste0("data/1_occ/calflora_downloads_030626/", names[i], "_calflora_highFilter.csv")) %>%
            #Fix the phenology_code column type
            mutate(phenology_code = ifelse(phenology_code==F, "F", phenology_code) %>% as.character()) %>%
            mutate(species = names[i]) %>%
            group_by(species, phenology_code) %>%
            count() %>%
            ungroup()
    }
) %>%
    bind_rows() %>%
    pivot_wider(names_from = phenology_code, values_from = n)

map(
    1:length(names),
    function(i){
        read.csv(paste0("data/1_occ/gbif_downloads_122625/", names[i], "/", names[i], "_gbif_lowFilter.csv")) %>%
            #Fix the phenology_code column type
            mutate(phenology_code = ifelse(phenology_code==F, "F", phenology_code) %>% as.character()) %>%
            mutate(species = names[i]) %>%
            group_by(species, phenology_code) %>%
            count() %>%
            ungroup()
    }
) %>%
    bind_rows() %>%
    pivot_wider(names_from = phenology_code, values_from = n)


## Merge and Spatially Thin
set.seed(123)
## file paths
parent_path_gbif <- "data/1_occ/gbif_downloads_122625/"
parent_path_calflora <- "data/1_occ/calflora_downloads_030626/"

## reference raster for thinning
rast <- rast('/home/tjsipin/will_tj/data/0_env/bcm/bcmv8_historic/2000_2023_monthly/aet2020dec.tif') %>%
  ## match crs to spp occ data
  project(y = "WGS84")

dir.create("data/1_occ/combined_spp_occ", recursive=T)

# Read in, merge, thin, and export for each spp on list
for (i in 1:length(names)) {
    print(i)
    ## GBIF ---------------------
    path_gbif = paste0(parent_path_gbif, names[i], "/")
    ### If no data for that species, doesn't contribute to merged file
    if (length(list.files(path_gbif, pattern = paste0(names[i], "_gbif_lowFilter.csv"))) == 0) {
        gbif = NULL
    } else {
        gbif <- read_csv(paste0(path_gbif, names[i], "_gbif_lowFilter.csv")) %>%
            ## only select vars of interest
            dplyr::select(c(gbifid, decimallatitude, decimallongitude, eventdate)) %>%
            ## consistent var names
            rename(id = gbifid,
                   lat = decimallatitude,
                   lon = decimallongitude,
                   date = eventdate) %>%
            ## add source
            mutate(source = "gbif",
                   ## remove time from date
                   date = as.Date(date)) %>%
            ## remove points after 2023 wy
            dplyr::filter(date < "2023-11-01")
    }

    ## CalFlora -----------------
    if (length(list.files(parent_path_calflora, pattern = names[i])) == 0) {
        calflora = NULL
    } else {
        calflora_date = read_csv(paste0(parent_path_calflora, names[i], "_calflora_highFilter.csv")) %>%
            pull(date) %>%
            is.Date()
        if(calflora_date){
            calflora <- read_csv(paste0(parent_path_calflora, names[i], "_calflora_highFilter.csv")) %>%
                dplyr::select(c(id, latitude, longitude, date)) %>%
                rename(lat = latitude,
                       lon = longitude) %>%
                mutate(source = "calflora") %>%
                dplyr::filter(date < "2023-11-01")
        } else{
            calflora <- read_csv(paste0(parent_path_calflora, names[i], "_calflora_highFilter.csv")) %>%
                dplyr::select(c(id, latitude, longitude, date)) %>%
                mutate(date = mdy(date)) %>%
                rename(lat = latitude,
                       lon = longitude) %>%
                mutate(source = "calflora") %>%
                dplyr::filter(date < "2023-11-01")
        }
    }

    ## Merge and thin -----------
    combo <- rbind(gbif, calflora) %>%
        ## Create sep vars for yr and mo
        mutate(year = lubridate::year(date),
               month = lubridate::month(date),
               .before = source)


    ###
    ## Only keep 1pt per rast cell
    comboThin <- elimCellDuplicates(combo, rast, longLat = c("lon", "lat"))

    ## Export
    write_csv(comboThin, paste0("data/1_occ/combined_spp_occ//",
                                names[i],
                                "_highFilter.csv"))
} ## END LOOP


# Background points

##Generate background points without month, which will be assigned later
source('src_highFilter/util/generate_backOcc.R')
dir.create('data_highFilter/2_background', recursive=T)
### 50km buffer
buffer = 50000
CA = tigris::states() %>%
    vect() %>%
    project("EPSG:3310") %>%
    filter(STUSPS=="CA")
### reference raster
referenceRaster = rast("data/0_env/natsgo/rasters/natsgo_drclass_270m_CA_2023.tif") %>%
    crop(CA, mask=T) %>%
    project('epsg:4326')
purrr::map(
    .x = 1:length(names),
    .f = function(i) {
        print(i)
        set.seed(123)
        ## read in spp occurrence points
        sppOcc = read_csv(paste0("data/1_occ/combined_spp_occ//",
                                 names[i],
                                 "_highFilter.csv"))

        ## Generate pts w/fxn
        backOcc_pts <- backOcc(sppOcc, raster=referenceRaster, buffer=buffer)

        ## Save
        write_csv(backOcc_pts, paste0("data_highFilter/2_background/back_",
                                      names[i],
                                      "_5km_highFilter.csv"))
    }
)


## Read in fxn and set parameters --------------------------
#Use this for wy model
source("src_highFilter/util/generateMonths_backOcc_evenMonth.R")


## Generate backOccs for each spp in list -------------------
purrr::map(
    .x = 1:length(names),
    .f = function(i) {
        set.seed(123)
        print(i)

        even_backOcc_pts.filename = paste0("data_highFilter/2_background/back_", names[i], "_5km_highFilter_even.csv")
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


###Water year####
# Read in fxn
source('src_highFilter/util/env_extract_wy.R')
dir.create('data_highFilter/3_swd/wy', recursive=T)

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
        out.dir = paste0('data_highFilter/3_swd/wy')
        #Get occ and back extracted sets
        in.filename = paste0(out.dir, '/swd_', sp, '_soil200cm_highFilter_wy.csv')
        
        training.filename = paste0(out.dir, '/training_', sp, '_soil200cm_highFilter_wy.csv')
        testing.filename = paste0(out.dir, '/testing_', sp, '_soil200cm_highFilter_wy.csv')
        
        trainingTestingFunc(
            sp = sp,
            in.filename = in.filename,
            training.filename = training.filename,
            testing.filename = testing.filename
        )
    }
)