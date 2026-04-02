library(tidyverse)    ## always
library(here)         ## reading/writing data
library(purrr)        ## iterate fxns
library(sf)           ## vector objects
library(terra)        ## Better/faster GIS
library(raster)       ## GIS format required by Dismo
library(dismo)        ## Maxent pkg
library(rJava)        ## Needed for dismo
library(lubridate)    ## Dates and progress bar
library(corrplot)     ## Correlation matrix
library(pROC)
library(tidymodels)
library(ENMeval)
library(tidyterra)
tidymodels_prefer()

set.seed(123)

## Read in fxn
source('src_lowFilter_v2/util/pred_agg.R')

central_valley = vect('data/central_valley/ds2632.gdb') %>% 
    project('epsg:3310')
CA = tigris::states() %>% 
    vect() %>% 
    project('epsg:3310') %>% 
    filter(STUSPS=="CA")


spp_relevant_months = tibble(
    spp = c(
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

# ## Loop fxn for all spp
# ###Low filter
# map(
#     names[1:10],
#     function(sp){
#         print(paste0("Working on ", sp))
#         model = readRDS(paste0("data_v2/4_maxent_outputs/agg/", sp,
#                                "/lowFilter/model/",
#                                sp, "_final_sdm.rds"))
#         out.path = paste0("data_v2/4_maxent_outputs/agg/",
#                           sp,
#                           "/lowFilter/monthly_dist_hist/")
#         dir.create(out.path, recursive = T)
#         months_vec = spp_relevant_months %>%
#             filter(spp == sp) %>%
#             pull(mons) %>%
#             unlist()
#         ## Run fxn
#         pred_agg(
#             model = model,
#             spp = sp,
#             months_vec = months_vec,
#             model_years = "2000_2023",
#             bcmPath = "data/0_env/bcm/bcmv8_historic/monthly_avgs/",
#             soilPath = "data/0_env/natsgo/rasters/",
#             salinityPath = 'data/0_env/salinity/',
#             pathOut = out.path
#         )
#     }
# )

#RCP4.5
map(
    names,
    function(sp){
        ## Read in Maxent model
        print(paste0("Working on ", sp))
        model = readRDS(paste0("data_v2/4_maxent_outputs/agg/", sp,
                               "/lowFilter/model/",
                               sp, "_final_sdm.rds"))
        out.path = paste0("data_v2/4_maxent_outputs/agg//",
                          sp,
                          "/lowFilter/monthly_dist_MIROC45/")
        dir.create(out.path, recursive = T)
        months_vec = spp_relevant_months %>%
            filter(spp == sp) %>%
            pull(mons) %>%
            unlist()
        ## Run fxn
        pred_agg(
            model = model,
            spp = sp,
            months_vec = months_vec,
            model_years = "MIROC45_2070_2099",
            bcmPath = "data/0_env/bcm/bcm_future/MIROC45/resampled/",
            soilPath = "data/0_env/natsgo/rasters/",
            salinityPath = 'data/0_env/salinity/',
            pathOut = out.path
        )
    }
)

#RCP8.5
map(
    names,
    function(sp){
        ## Read in Maxent model
        print(paste0("Working on ", sp))
        model = readRDS(paste0("data_v2/4_maxent_outputs/agg/", sp,
                               "/lowFilter/model/",
                               sp, "_final_sdm.rds"))
        out.path = paste0("data_v2/4_maxent_outputs/agg/",
                          sp,
                          "/lowFilter/monthly_dist_MIROC85/")
        dir.create(out.path, recursive = T)
        months_vec = spp_relevant_months %>%
            filter(spp == sp) %>%
            pull(mons) %>%
            unlist()
        ## Run fxn
        pred_agg(
            model = model,
            spp = sp,
            months_vec = months_vec,
            model_years = "MIROC85_2070_2099",
            bcmPath = "data/0_env/bcm/bcm_future/MIROC85/resampled/",
            soilPath = "data/0_env/natsgo/rasters/",
            salinityPath = 'data/0_env/salinity/',
            pathOut = out.path
        )
    }
)


# Threshold rasters -------------------------------------------------------

## P10 threshold raster per species and model-years combination
getThresholdRasts = function(sp, model_years=c("2000_2023", "MIROC45_2070_2099", "MIROC85_2070_2099")){
    #Read in training data and model
    training.filename = paste0('data_v2/3_swd/agg/training_', sp, '_soil200cm_lowFilter_agg.csv')
    testing.filename = paste0('data_v2/3_swd/agg/testing_', sp, '_soil200cm_lowFilter_agg.csv')
    model.filename = paste0("data_v2/4_maxent_outputs/agg/", sp, "/lowFilter/model/", sp, "_final_sdm.rds")
    
    training = read_csv(training.filename)
    testing = read_csv(testing.filename)
    model = readRDS(model.filename)
    
    #Subset training data to occurrences only
    pred1 = training %>%
        rbind(testing) %>%
        filter(presence==1) %>%
        mutate(tdiff = tmx - tmn) %>%
        filter(complete.cases(.)) %>%
        mutate(pred = predict(model, .))
    
    #Get the 10th percentile of occurrence point predictions
    p10 = quantile(pred1$pred, 0.1)
    
    #Designate filepath string
    if(str_detect(model_years, "MIROC45")){
        monthly_dist_str = "monthly_dist_MIROC45/"
    } else if(str_detect(model_years, "MIROC85")){
        monthly_dist_str = "monthly_dist_MIROC85/"
    } else if(str_detect(model_years, "2023")){
        monthly_dist_str = "monthly_dist_hist/"
    }
    
    #List input prediction rast filenames
    input.filenames = list.files(
        paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/', monthly_dist_str),
        full.names = T,
        pattern='.tif'
    )
    
    #Create binary 10th percentile rasts
    output.rast = rast(input.filenames)
    output.rast = ifel(output.rast > p10, 1, 0)
    output.rast = output.rast %>% sum()
    
    dir.create(paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/p10/'), recursive=T, showWarnings=F)
    output.filename = paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/p10/', '/p10_', sp, '_', model_years, '_agg.tif')
    writeRaster(output.rast, output.filename, overwrite = T)
}

map(
    .x = names,
    .f = ~ getThresholdRasts(sp = .x, '2000_2023')
)
map(
    .x = names,
    .f = ~ getThresholdRasts(sp = .x, 'MIROC45_2070_2099')
)
map(
    .x = names,
    .f = ~ getThresholdRasts(sp = .x, 'MIROC85_2070_2099')
)


# To use in Figure 3 ------------------------------------------------------

sumThresholdRasts = function(model){
    out.dir = paste0("data_v2/4_maxent_outputs/agg/p10/")
    dir.create(out.dir, recursive=T)
    out.rast.filename.CA = paste0("data_v2/4_maxent_outputs/agg/p10/", model, "_lowFilter_sum_CA.tif")
    out.rast.filename.CV = paste0("data_v2/4_maxent_outputs/agg/p10/", model, "_lowFilter_sum_CV.tif")
    #Read in p10 rasters
    input_rast.filenames = list.files(
        paste0("data_v2/4_maxent_outputs/agg/"),
        pattern=paste0("p10_.*", model, '_agg'),
        recursive = T,
        full.names=T
    )
    #Sum across CA
    input_rast_CA = rast(input_rast.filenames) %>%
        sum() %>%
        mutate(sum = factor(sum, levels = as.integer(0:10))) %>%
        crop(CA, mask=T)
    #Sum across CV
    input_rast_CV = input_rast_CA %>%
        crop(central_valley, mask=T)
    
    writeRaster(input_rast_CA, out.rast.filename.CA, overwrite = T)
    writeRaster(input_rast_CV, out.rast.filename.CV, overwrite = T)
}

sumThresholdRasts(model = "2000_2023")
sumThresholdRasts(model = "MIROC45_2070_2099")
sumThresholdRasts(model = "MIROC85_2070_2099")
