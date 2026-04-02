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
library(tidyterra)
library(ENMeval)
tidymodels_prefer()

message('v5/mc/04')

## Read in fxn
source('src_lowFilter_v5_50km/util/pred_month.R')

central_valley = vect('data/central_valley/ds2632.gdb') %>% 
    project('epsg:3310')
CA = tigris::states() %>% 
    vect() %>% 
    project('epsg:3310') %>% 
    filter(STUSPS=="CA")

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

# Loop fxn for all spp
#Low filter
map(
    names,
    function(sp){
        ## Read in Maxent model
        print(paste0("Working on ", sp))
        model = readRDS(paste0("data_v5_50km/4_maxent_outputs/comprehensive/", sp,
                               "/lowFilter/model/",
                               sp, "_final_sdm.rds"))
        out.path = paste0("data_v5_50km/4_maxent_outputs/comprehensive/",
                          sp,
                          "/lowFilter/monthly_dist_hist/")
        dir.create(out.path, recursive = T)
        ## Run fxn
        pred_month(model = model,
                   spp = sp,
                   model_years = "2000_2023",
                   bcmPath = "data/0_env/bcm/bcmv8_historic/monthly_avgs/",
                   soilPath = "data/0_env/natsgo/rasters/",
                   pathOut = out.path
        )
    }
)

#RCP4.5
map(
    names,
    function(sp){
        ## Read in Maxent model
        print(paste0("Working on ", sp))
        model = readRDS(paste0("data_v5_50km/4_maxent_outputs/comprehensive/", sp, "/lowFilter/model/", sp, "_final_sdm.rds"))
        out.path = paste0("data_v5_50km/4_maxent_outputs/comprehensive//",
                          sp,
                          "/lowFilter/monthly_dist_MIROC45/")
        dir.create(out.path, recursive = T)
        ## Run fxn
        pred_month(model = model,
                   spp = sp,
                   model_years = "MIROC45_2070_2099",
                   bcmPath = "data/0_env/bcm/bcm_future/MIROC45/resampled/",
                   soilPath = "data/0_env/natsgo/rasters/",
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
        model = readRDS(paste0("data_v5_50km/4_maxent_outputs/comprehensive/", sp, "/lowFilter/model/", sp, "_final_sdm.rds"))
        out.path = paste0("data_v5_50km/4_maxent_outputs/comprehensive//",
                          sp,
                          "/lowFilter/monthly_dist_MIROC85/")
        dir.create(out.path, recursive = T)
        ## Run fxn
        pred_month(model = model,
                   spp = sp,
                   model_years = "MIROC85_2070_2099",
                   bcmPath = "data/0_env/bcm/bcm_future/MIROC85/resampled/",
                   soilPath = "data/0_env/natsgo/rasters/",
                   pathOut = out.path
        )
    }
)

# For figure 4 ------------------------------------------------------------


## P10 threshold per species, month, and model-years combination
getThresholdRasts = function(sp, model_years=c("2000_2023", "MIROC45_2070_2099", "MIROC85_2070_2099")){
    print(paste0(sp, "\n", model_years))
    dir.create(paste0('data_v5_50km/4_maxent_outputs/comprehensive/', sp, '/lowFilter/p10/'), recursive=T, showWarnings=F)
    output.filename.CV = paste0('data_v5_50km/4_maxent_outputs/comprehensive/', sp, '/lowFilter/p10/', '/p10_', sp, '_', model_years, '_monthly_CV.tif')
    output.filename.CA = paste0('data_v5_50km/4_maxent_outputs/comprehensive/', sp, '/lowFilter/p10/', '/p10_', sp, '_', model_years, '_monthly_CA.tif')
    
    #Read in training data and model
    training.filename = paste0('data_v5_50km/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')
    testing.filename = paste0('data_v5_50km/3_swd/monthly/comprehensive/testing_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')
    model.filename = paste0("data_v5_50km/4_maxent_outputs/comprehensive/", sp, "/lowFilter/model/", sp, "_final_sdm.rds")
    
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
    
    if(str_detect(model_years, "MIROC45")){
        monthly_dist_str = "monthly_dist_MIROC45/"
    } else if(str_detect(model_years, "MIROC85")){
        monthly_dist_str = "monthly_dist_MIROC85/"
    } else if(str_detect(model_years, "2023")){
        monthly_dist_str = "monthly_dist_hist/"
    }
    sp_theme = sp %>% 
        str_replace("_", ". ") %>% 
        str_to_sentence()
    
    input.filenames = list.files(
        paste0('data_v5_50km/4_maxent_outputs/comprehensive/', sp, '/lowFilter/', monthly_dist_str),
        full.names = T
    )
    map(
        input.filenames,
        function(fi){
            output.filename = str_replace(fi, 'monthly_dist', 'p10/p10')
            output.dir = str_split_i(output.filename, "//", 1)
            dir.create(output.dir)
            input.rast = rast(fi)
            thresh.rast = ifel(input.rast > p10, 1, 0)
            writeRaster(thresh.rast, output.filename, overwrite=T)
        },
        .progress = T
    )
    
    thresh.filenames = list.files(
        paste0("data_v5_50km/4_maxent_outputs/comprehensive/", sp, "/lowFilter/p10/"),
        pattern = paste0(model_years, ".tif"),
        recursive=T,
        full.names = T
    )
    output.rast = rast(thresh.filenames)
    output.rast = ifel(output.rast > p10, 1, 0) %>% 
        sum()
    output.rast.CA = output.rast %>% 
        crop(CA, mask=T)
    output.rast.CV = output.rast %>% 
        crop(central_valley, mask=T)
    
    writeRaster(output.rast.CV, output.filename.CV, overwrite = T)
    writeRaster(output.rast.CA, output.filename.CA, overwrite = T)
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

