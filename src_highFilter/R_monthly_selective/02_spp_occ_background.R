library(tidyverse)            ## always
# library(here)                 ## reading/writing data
library(rgbif)                ## download GBIF data
library(CoordinateCleaner)    ## filter/clean GBIF data
library(terra)                ## rast pkg for quicker reprojecting
library(raster)               ## rast format that plays w/dismo
library(enmSdmX)              ## spatially thinning data
library(dismo)                ## generating background points
library(lfstat)               ## water year fxn
library(tidymodels)
library(tidyterra)
tidymodels_prefer()
setwd('/home/tjsipin/will_tj')

#' This script generates and preps occurrence and background data for running the SDMs.
#' First, occurrence records from GBIF and CalFlora are filtered, merged, and spatially thinned.
#' Second, a set of random background points are created for each species.
#' Finally, environmental data are extracted for all occurrence and background points.
#' These "samples with data" files are used for the models in the `kern_sdm` markdown.


# Occurrence Records

## Download/import

### Phenology dictionary

repro_values <- read_csv("data/occ/gbif_codes.csv")



### CalFlora
path = "data/occ/calflora_downloads_072425/"
names <- c(
    "a_polycarpa",
    "p_arborea",
    "c_pungens",
    "l_pentachaeta",
    "p_ciliata",
    "a_menziesii",
    "c_lasiophyllus"
)

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
    ## Read in data
    df <- read_csv(paste0(path, names[i], "_calflora_071525.csv")) %>%
        janitor::clean_names()

    ## filter and keep coords
    df_filter <- df %>%
        ## "LOW" filter criteria
        filter(location_quality %in% c("high", "medium"),
               accuracy_square_meters <= 72900,
               date >= "1999-10-01") #%>%
    ## "HIGH" filter criteria
    ## Either from CCH or one of the other sources
    # filter(str_detect(.$source, sources) | dataset == "cch2")

    df_repro_filter = df_filter %>%
        filter(phenology_code %in% c("F", "R") | is.na(phenology_code)) %>%
        mutate(phenology_code = as.character(phenology_code))

    ## export
    write.csv(df_filter, paste0(path, names[i], "_calflora_lowFilter_noReproFilter.csv"))
    write.csv(df_repro_filter, paste0(path, names[i], "_calflora_lowFilter.csv"))
}


### GBIF
## Pull taxon keys from list of spp
taxon_keys <- name_backbone_checklist(c(
    ## shrubs
    "Atriplex polycarpa",
    "Peritoma arborea",
    ## forbs
    "Centromadia pungens",
    "Layia pentachaeta subsp. albida D.D.Keck",
    "Phacelia ciliata Benth.",
    "Amsinckia menziesii (Lehm.) A.Nelson & J.F.Macbr.",
    "Caulanthus lasiophyllus (Hook. & Arn.) Payson"
)) %>%
    pull(usageKey)


names <- c(
    "a_polycarpa",
    "p_arborea",
    "c_pungens",
    "l_pentachaeta",
    "p_ciliata",
    "a_menziesii",
    "c_lasiophyllus"
)


## Loop through spp list and save each
for (i in 1:length(taxon_keys)) {
    #CAUTION: Do not run as this may overwrite our gbif downloads!
                # ### download filtered data for A. polycarpa
                # gbif_download <- occ_download(pred_in("taxonKey", taxon_keys[i]),
                #                        ## remove geospatial issues
                #                        pred("hasGeospatialIssue", FALSE),
                #                        ## ensure coords
                #                        pred("hasCoordinate", TRUE),
                #                        ## remove "absent" occurences
                #                        pred("occurrenceStatus", "PRESENT"),
                #                        ## within US
                #                        pred("country", "US"),
                #                        ## within CA
                #                        pred("stateProvince", "California"),
                #                        ## output as CSV
                #                        format = "SIMPLE_CSV",
                #                        ## enter your GBIF credentials below:
                #                         user = "tjsipin_ucsb",
                #                         pwd = "RKh4pXLx",
                #                         email = "tjsipin@ucsb.edu"
                #                        )
                #
                # ### check on download status
                # occ_download_wait(gbif_download)

    ### import GBIF data into env and filter using CoordinateCleaner pkg
    species <- list.files("data/occ/gbif_downloads_072425/zips/")[i] %>%
        str_sub(1, -5) %>%
        occ_download_get(
            path = ("/home/tjsipin/will_tj/data/occ/gbif_downloads_072425/extracts/"),
            overwrite = TRUE
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
        #TJ: This is changed from Nick/Madi/Will
        mutate(reproductivecondition_clean = tolower(str_trim(reproductivecondition))) %>%
        left_join(repro_values)

    #TJ: This is also changed from Nick/Madi/Will
    species_repro_filter = species %>%
        filter(phenology_code %in% c("F", "R") | is.na(phenology_code))

    ## export file
    write_csv(species, paste0(here("data/occ/gbif_downloads_072425/"),names[i],"/", names[i], "_gbif_lowFilter_noReproFilter.csv"))
    write_csv(species_repro_filter, paste0(here("data/occ/gbif_downloads_072425/"),names[i],"/", names[i], "_gbif_lowFilter.csv"))

} ## END LOOP

#See lifestage counts per species
map(
    1:7,
    function(i){
        read_csv(paste0(here("data/occ/calflora_downloads_072425/"),names[i], "_calflora_lowFilter.csv")) %>%
            #Fix the phenology_code column type
            mutate(phenology_code = ifelse(phenology_code==F, "F", phenology_code)) %>%
            mutate(species = names[i]) %>%
            group_by(species, phenology_code) %>%
            count() %>%
            ungroup()
    }
) %>%
    bind_rows() %>%
    pivot_wider(names_from = phenology_code, values_from = n)

map(
    1:7,
    function(i){
        read.csv(paste0(here("data/occ/gbif_downloads_072425/"),names[i],"/", names[i], "_gbif_lowFilter.csv")) %>%
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
parent_path_gbif <- "/home/tjsipin/will_tj/data/occ/gbif_downloads_072425/"
parent_path_calflora <- "/home/tjsipin/will_tj/data/occ/calflora_downloads_072425/"

## reference raster for thinning
rast <- rast(here('/home/tjsipin/madi/SDMplants/data/bcm/bcmv8_historic/2000_2023_monthly/aet2020dec.tif')) %>%
  ## match crs to spp occ data
  project(y = "WGS84")

dir.create("data/occ/combined_spp_occ", recursive=T)

## spp names
names <- c(
  "a_polycarpa",
  "p_arborea",
  "c_pungens",
  "l_pentachaeta",
  "p_ciliata",
  "a_menziesii",
  "c_lasiophyllus"
)

## Read in, merge, thin, and export for each spp on list
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
if (length(list.files(path_calflora, pattern = names[i])) == 0) {
    calflora = NULL
} else {
    if(i==1){
        calflora <- read_csv(paste0(path_calflora, names[i], "_calflora_lowFilter.csv")) %>%
            dplyr::select(c(id, latitude, longitude, date)) %>%
            mutate(date = mdy(date)) %>%
            rename(lat = latitude,
                   lon = longitude) %>%
            mutate(source = "calflora") %>%
            dplyr::filter(date < "2023-11-01")
    } else{
        calflora <- read_csv(paste0(path_calflora, names[i], "_calflora_lowFilter.csv")) %>%
            dplyr::select(c(id, latitude, longitude, date)) %>%
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
comboThin <- elimCellDuplicates(combo, rast, longLat = c("lon", "lat")) %>%
    mutate()

## Export
write_csv(comboThin, paste0(here("data/occ/combined_spp_occ//"),
                            names[i],
                            "_lowFilter.csv"))
} ## END LOOP




# Background points

## Read in fxn and set parameters --------------------------
source(here("R/generate_backOcc_tj.R"))
set.seed(123)
### 5km buffer
buffer = 5000
### reference raster
raster = rast(("/home/tjsipin/madi/SDMplants/data/natsgo/rasters/natsgo_ph_270m_CA_2023.tif")) #change reference raster from a gnatsgo one....'lowFilter' has no significance
### spp names
names <- c(
    "a_polycarpa",
    "p_arborea",
    "c_pungens",
    "l_pentachaeta",
    "p_ciliata",
    "a_menziesii",
    "c_lasiophyllus"
)

## Generate backOccs for each spp in list -------------------
purrr::map(
    .x = 1:7,
    .progress = TRUE,
    .f = function(i) {
        ## read in spp occurrence points
        sppOcc = read_csv(paste0(here("data/occ/combined_spp_occ//"),
                                 names[i],
                                 "_lowFilter.csv"))

        ## Generate pts w/fxn
        backOcc_pts <- backOcc(sppOcc, raster=raster, buffer=buffer)

        ## Save
        write_csv(backOcc_pts, paste0(here("data/background/back_"),
                                      names[i],
                                      "_5km_lowFilter.csv"))
    }
)


# Extract environmental data

## Each spp being modeled
names <- c(
  "a_polycarpa",
  "p_arborea",
  "c_pungens",
  "l_pentachaeta",
  "p_ciliata",
  "a_menziesii",
  "c_lasiophyllus"
)

## Read in fxn
source('/home/tjsipin/will_tj/R/env_extract_tj_3.R')
spp_occ_back = tibble(
    spp = rep(c("a_polycarpa", "p_arborea", "c_pungens", "l_pentachaeta", "p_ciliata", "a_menziesii", "c_lasiophyllus"), 2),
    occ_back = c(rep('occ', 7), rep('back', 7))
)

## Fxn variables
startYear = 2000
endYear = 2023
pathMonth = "/home/tjsipin/madi/SDMplants/data/bcm/bcmv8_historic/2000_2023_monthly/"
pathQuarter = "/home/tjsipin/madi/SDMplants/data/bcm/bcmv8_historic/quarterly_avgs/"
pathSoil = "/home/tjsipin/madi/SDMplants/data/natsgo/rasters/"

dir.create('data/swd_072925')

map(
    1:14,
    function(i){
        spp = spp_occ_back$spp[i]
        occ_back = spp_occ_back$occ_back[i]
        out.filename = paste0('data/swd/', spp, '/', occ_back, 'Extract_', spp, '_gNatsgo_lowFilter072925.csv')
        if(file.exists(out.filename)) return(NULL)
        if(occ_back=='occ'){
            points.filename = paste0('data/occ/combined_spp_occ/', spp, '_lowFilter.csv')
        } else{
            points.filename = paste0("data/background/back_", spp, "_5km_lowFilter.csv")
        }
        points = read_csv(points.filename)

        print(spp)
        out = env_extract(
            startYear = startYear,
            endYear = endYear,
            pathMonth = pathMonth,
            pathQuarter = pathQuarter,
            pathSoil = pathSoil,
            occ = points
        )

        write_csv(out, out.filename)

    },
    .progress = T
)

map(
    names,
    trainingTestingFunc
)


#####August 5, 2025#####
## Read in fxn
source('/home/tjsipin/will_tj/original/R/env_extract_tj_4.R')
spp_occ_back = tibble(
    spp = rep(c("a_polycarpa", "p_arborea", "c_pungens", "l_pentachaeta", "p_ciliata", "a_menziesii", "c_lasiophyllus"), 2),
    occ_back = c(rep('occ', 7), rep('back', 7))
)

## Fxn variables
startYear = 2000
endYear = 2023
pathMonth = "/home/tjsipin/madi/SDMplants/data/bcm/bcmv8_historic/2000_2023_monthly/"
pathQuarter = "/home/tjsipin/madi/SDMplants/data/bcm/bcmv8_historic/quarterly_avgs/"
pathSoil = "/home/tjsipin/madi/SDMplants/data/natsgo/rasters/"

dir.create('data/swd_072925')

map(
    1:14,
    function(i){
        spp = spp_occ_back$spp[i]
        occ_back = spp_occ_back$occ_back[i]
        out.filename = paste0('data/swd/', spp, '/', occ_back, 'Extract_', spp, '_gNatsgo_lowFilter080525.csv')
        if(file.exists(out.filename)) return(NULL)
        if(occ_back=='occ'){
            points.filename = paste0('data/occ/combined_spp_occ/', spp, '_lowFilter.csv')
        } else{
            points.filename = paste0("data/background/back_", spp, "_5km_lowFilter.csv")
        }
        points = read_csv(points.filename) 
        
        print(spp)
        out = env_extract(
            startYear = startYear, 
            endYear = endYear,
            pathMonth = pathMonth,
            pathQuarter = pathQuarter,
            pathSoil = pathSoil,
            occ = points
        )
        
        write_csv(out, out.filename)
        
    },
    .progress = T
)

map(
    names,
    function(spp){
        print(spp)
        out.dir = paste0('data/swd/', spp)
        #Get occ and back extracted sets
        occ.filename = paste0(out.dir, '/occExtract_', spp, '_gNatsgo_lowFilter080525.csv')
        back.filename = paste0(out.dir, '/backExtract_', spp, '_gNatsgo_lowFilter080525.csv')
        
        training.filename = paste0(out.dir, '/training_', spp, '_gNatsgo_lowFilter080525.csv')
        testing.filename = paste0(out.dir, '/testing_', spp, '_gNatsgo_lowFilter080525.csv')
        trainingTestingFunc(
            spp=spp, 
            occ.filename=occ.filename, 
            back.filename=back.filename,
            training.filename=training.filename, 
            testing.filename=testing.filename
        )
    }
)
