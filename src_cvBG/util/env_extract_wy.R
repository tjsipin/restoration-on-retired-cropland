#This version, we will be extracting the previous 12 months (+ conccurrent) of environmental data

# setwd("/home/tjsipin/will_tj")

library(tidyverse)            ## always
library(rgbif)                ## download GBIF data
library(CoordinateCleaner)    ## filter/clean GBIF data
library(terra)                ## rast pkg for quicker reprojecting
library(raster)               ## rast format that plays w/dismo
library(enmSdmX)              ## spatially thinning data
library(dismo)                ## generating background points
library(lfstat)               ## water year fxn
library(tidyterra)            ## dplyr functions work on terra objects
library(rsample)
tidymodels::tidymodels_prefer()

CA <- tigris::states() %>% 
    vect() %>% 
    filter(STUSPS=="CA") %>% 
    project("EPSG:3310")

relevantWY = function(
        points,
        this_wy,
        bcm_filename_df,
        input_path='data/0_env/bcm/bcmv8_historic/2000_2023_monthly/', 
        target_crs, sp
){
    
    #Gather files
    bcm_filenames = list.files(
        "data/0_env/bcm/bcmv8_historic/2000_2023_monthly/",
        pattern = ".tif"
    ) 
    
    
    
    # Jan 2 -------------------------------------------------------------------
    
    bcm_tibble = bcm_filename_df %>%
        mutate(filename = paste0(input_path, filename)) %>% 
        arrange(var, year, month_num) %>%
        mutate(wy = lfstat::water_year(x = my(paste0(month_num, '-', year)), origin = 10)) %>% 
        group_split(var) %>%
        map(
            function(df){
                out = df %>% 
                    filter(wy == this_wy) 
                
                if(nrow(out) < 12) return(NULL)
                
                return(out)
            }
        ) %>%
        bind_rows() 
    if(nrow(bcm_tibble)==0) return(NULL)
    
    #For each variable, pull appropriate file names 
    aet_filenames = bcm_tibble %>% 
        filter(var == "aet") %>% 
        pull(filename)
    cwd_filenames = bcm_tibble %>% 
        filter(var == "cwd") %>% 
        pull(filename)
    
    ppt_filenames = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        pull(filename)
    ppt_10_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==10) %>% 
        pull(filename)
    ppt_11_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==11) %>% 
        pull(filename)
    ppt_12_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==12) %>% 
        pull(filename)
    ppt_1_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==1) %>% 
        pull(filename)
    ppt_2_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==2) %>% 
        pull(filename)
    ppt_3_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==3) %>% 
        pull(filename)
    ppt_4_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==4) %>% 
        pull(filename)
    ppt_5_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==5) %>% 
        pull(filename)
    ppt_6_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==6) %>% 
        pull(filename)
    ppt_7_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==7) %>% 
        pull(filename)
    ppt_8_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==8) %>% 
        pull(filename)
    ppt_9_filename = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        filter(month_num==9) %>% 
        pull(filename)
    
    tmx_filenames = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        pull(filename)
    tmx_10_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==10) %>% 
        pull(filename)
    tmx_11_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==11) %>% 
        pull(filename)
    tmx_12_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==12) %>% 
        pull(filename)
    tmx_1_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==1) %>% 
        pull(filename)
    tmx_2_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==2) %>% 
        pull(filename)
    tmx_3_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==3) %>% 
        pull(filename)
    tmx_4_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==4) %>% 
        pull(filename)
    tmx_5_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==5) %>% 
        pull(filename)
    tmx_6_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==6) %>% 
        pull(filename)
    tmx_7_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==7) %>% 
        pull(filename)
    tmx_8_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==8) %>% 
        pull(filename)
    tmx_9_filename = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(month_num==9) %>% 
        pull(filename)
    
    tmn_filenames = bcm_tibble %>% 
        filter(var == "tmn") %>% 
        pull(filename)
    
    #gNATSGO 
    ph_filename = paste0("data/0_env/natsgo/rasters/natsgo_ph_270m_CA_2023.tif")
    om_filename = paste0("data/0_env/natsgo/rasters/natsgo_om_270m_CA_2023.tif")
    cec_filename = paste0("data/0_env/natsgo/rasters/natsgo_cec_270m_CA_2023.tif")
    drclass_filename = paste0("data/0_env/natsgo/rasters/natsgo_drclass_270m_CA_2023.tif")
    
    #POLARIS
    bd_POLARIS_filename = paste0('data/0_env/polaris/rasters/bd_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif')
    clay_POLARIS_filename = paste0('data/0_env/polaris/rasters/clay_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif')
    om_POLARIS_filename = paste0('data/0_env/polaris/rasters/om_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif')
    ph_POLARIS_filename = paste0('data/0_env/polaris/rasters/ph_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif')
    sand_POLARIS_filename = paste0('data/0_env/polaris/rasters/sand_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif')
    silt_POLARIS_filename = paste0('data/0_env/polaris/rasters/silt_0_200_weightedAverages_scale270_bcmProjScale_reprojectFuncNoScale.tif')
    
    #Salinity
    salinity_filename = paste0("data/0_env/salinity/salinity_full_CA_res.tif")
    
    rast_filenames = list(
        aet_filenames,
        cwd_filenames, 
        ppt_filenames,
        tmx_filenames,
        tmn_filenames
    )
    
    for(f in rast_filenames){
        if(length(f)==0){
            message("Length of files for at least one BCM variable is 0.")
            return(NULL)
        } 
    }
    
    #Rasterize all sets of file names
    aet_annual = rast(aet_filenames) %>% 
        mean() %>% 
        crop(CA)
    cwd_annual = rast(cwd_filenames) %>% 
        sum() %>% 
        crop(CA)
    
    ppt_annual = rast(ppt_filenames) %>% 
        sum() %>% 
        crop(CA)
    ppt_10_annual = rast(ppt_10_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_11_annual = rast(ppt_11_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_12_annual = rast(ppt_12_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_1_annual = rast(ppt_1_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_2_annual = rast(ppt_2_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_3_annual = rast(ppt_3_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_4_annual = rast(ppt_4_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_5_annual = rast(ppt_5_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_6_annual = rast(ppt_6_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_7_annual = rast(ppt_7_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_8_annual = rast(ppt_8_filename) %>% 
        sum() %>% 
        crop(CA)
    ppt_9_annual = rast(ppt_9_filename) %>% 
        sum() %>% 
        crop(CA)
    
    
    tmx_annual = rast(tmx_filenames) %>% 
        mean() %>% 
        crop(CA)
    tmx_10_annual = rast(tmx_10_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_11_annual = rast(tmx_11_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_12_annual = rast(tmx_12_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_1_annual = rast(tmx_1_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_2_annual = rast(tmx_2_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_3_annual = rast(tmx_3_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_4_annual = rast(tmx_4_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_5_annual = rast(tmx_5_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_6_annual = rast(tmx_6_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_7_annual = rast(tmx_7_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_8_annual = rast(tmx_8_filename) %>% 
        mean() %>% 
        crop(CA)
    tmx_9_annual = rast(tmx_9_filename) %>% 
        mean() %>% 
        crop(CA)
    
    tmn_annual = rast(tmn_filenames) %>% 
        mean() %>% 
        crop(CA)
    
    tdiff_annual = tmx_annual - tmn_annual
    
    ph = rast(ph_filename) %>% 
        resample(aet_annual)
    om = rast(om_filename) %>% 
        resample(aet_annual)
    cec = rast(cec_filename) %>% 
        resample(aet_annual)
    drclass = rast(drclass_filename) %>% 
        resample(aet_annual, method='mode')
    
    bd_POLARIS = rast(bd_POLARIS_filename) %>% 
        resample(aet_annual)
    clay_POLARIS = rast(clay_POLARIS_filename) %>% 
        resample(aet_annual)
    om_POLARIS = rast(om_POLARIS_filename) %>% 
        resample(aet_annual)
    ph_POLARIS = rast(ph_POLARIS_filename) %>% 
        resample(aet_annual)
    sand_POLARIS = rast(sand_POLARIS_filename) %>% 
        resample(aet_annual)
    silt_POLARIS = rast(silt_POLARIS_filename) %>% 
        resample(aet_annual)
    
    salinity = rast(salinity_filename) %>% 
        resample(aet_annual, method='mode')
    
    #Combine all rasters
    all_annual = rast(list(
        aet_annual, cwd_annual, 
        ppt_annual, ppt_10_annual, ppt_11_annual, ppt_12_annual, 
        ppt_1_annual, ppt_2_annual, ppt_3_annual, ppt_4_annual,
        ppt_5_annual, ppt_6_annual, ppt_7_annual, ppt_8_annual, ppt_9_annual,
        tmx_annual, tmx_10_annual, tmx_11_annual, tmx_12_annual, 
        tmx_1_annual, tmx_2_annual, tmx_3_annual, tmx_4_annual,
        tmx_5_annual, tmx_6_annual, tmx_7_annual, tmx_8_annual, tmx_9_annual,
        tmn_annual, tdiff_annual, 
        ph, om, cec, drclass,
        bd_POLARIS, clay_POLARIS, om_POLARIS,
        ph_POLARIS, sand_POLARIS, silt_POLARIS, 
        salinity
    )) 
    #Rename the combined raster layers
    names(all_annual) = c(
        "aet", "cwd", 
        'ppt', 'ppt10', 'ppt11', 'ppt12',
        'ppt1', 'ppt2', 'ppt3', 'ppt4',
        'ppt5', 'ppt6', 'ppt7', 'ppt8', 'ppt9',
        'tmx', 'tmx10', 'tmx11', 'tmx12',
        'tmx1', 'tmx2', 'tmx3', 'tmx4',
        'tmx5', 'tmx6', 'tmx7', 'tmx8', 'tmx9',
        'tmn', 'tdiff',
        "ph", "om", "cec", "drclass",
        "bd_POLARIS", "clay_POLARIS", "om_POLARIS", 
        "ph_POLARIS", "sand_POLARIS", "silt_POLARIS", 
        "salinity"
    )
    
    #Vectorize the points
    points_y_vect = points %>% 
        filter(wy==this_wy) %>% 
        vect(geom=c("lon", "lat"), "EPSG:4326") %>% 
        terra::project(target_crs)
    
    #Extract the combined raster values at each point
    extracted = terra::extract(all_annual, points_y_vect, method = "simple", bind=T, xy=T) %>% 
        as.data.frame()
    return(extracted)
}

extractEnvWY = function(sp){
    
    occ_filename = paste0('data/1_occ/combined_spp_occ/', sp, '_lowFilter.csv')
    occ = read_csv(occ_filename) %>% 
        mutate(sp = sp, presence = 1) %>% 
        select(sp, lon, lat, month, year, presence)
    
    back_filename = paste0('data_v4_cvBG/2_background/back_', sp, '_5km_lowFilter_even.csv')
    back = read_csv(back_filename) %>% 
        mutate(sp = sp, presence = 0) %>% 
        select(sp, lon, lat, month, year, presence)
    
    out.dir = 'data_v4_cvBG/3_swd/wy'
    
    out.filename = paste0('data_v4_cvBG/3_swd/wy/', 'swd_', sp, '_soil200cm_lowFilter_wy.csv')
    # if(file.exists(out.filename)) return(NULL)
    
    #1) Read in the occurrence and background data 
    points = occ %>% 
        rbind(back) %>% 
        mutate(spp_id = row_number()) %>% 
        mutate(wy = water_year(my(paste0(month, '-', year)), origin = 10))
    
    target_crs = rast("data/0_env/bcm/bcmv8_historic/2000_2023_monthly/aet1999dec.tif") %>% 
        crs()
    
    #2) Create file grid of month-year-filename
    #Gather files
    bcm_filenames = list.files(
        "data/0_env/bcm/bcmv8_historic/2000_2023_monthly/",
        pattern = ".tif"
    )
    #Convert to data frame
    bcm_filename_df = tibble(
        filename = bcm_filenames
    ) %>%
        mutate(
            var = str_sub(filename, start = 1, end = 3),
            year = str_sub(filename, start = 4, end = 7) %>% as.integer(),
            month = str_sub(filename, start = 8, end = 10)
        ) %>%
        filter(var %in% c('aet', 'tmx', 'tmn', 'ppt', 'cwd')) %>%
        mutate(
            month_num = dplyr::recode(
                month,
                jan = 1,
                feb = 2,
                mar = 3,
                apr = 4,
                may = 5,
                jun = 6,
                jul = 7,
                aug = 8,
                sep = 9,
                oct = 10,
                nov = 11,
                dec = 12
            )
        ) 
    water_years = points %>% 
        pull(wy) %>% 
        unique() %>% 
        sort()
    
    spp_extract = map(
        .x = water_years,
        function(.x){
            message(.x)
            tryCatch({
                this_wy = as.numeric(as.character(.x))
                relevantWY(
                    points = points,
                    this_wy = this_wy, 
                    bcm_filename_df = bcm_filename_df,
                    target_crs = target_crs, sp = sp
                )
            }, error = function(e) {
                message(paste0(sp, .x))
                message(e)
                return(NULL)
            })
        },
        .progress = T
    ) 
    spp_extract_binded = spp_extract %>% 
        bind_rows() %>% 
        select(
            spp_id, x, y, month, year, 
            aet, tmx, tdiff, 
            ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
            ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
            tmx10, tmx11, tmx12, tmx1, tmx2, tmx3,
            tmx4, tmx5, tmx6, tmx7, tmx8, tmx9,
            cec, om, ph, drclass, salinity,
            presence
        ) %>% 
        distinct() %>% 
        mutate(
            drclass=as.factor(drclass),
            salinity=as.factor(salinity)
        ) %>% 
        rowwise() %>% 
        mutate(
            tmx_winter_mean = mean(c(tmx10, tmx11, tmx12, tmx1)),
            tmx_spring_mean = mean(c(tmx2, tmx3, tmx4, tmx5)),
            tmx_summer_mean = mean(c(tmx6, tmx7, tmx8, tmx9))
        ) %>% 
        group_by(spp_id) %>% 
        group_split() %>% 
        map(
            function(df){
                df %>% 
                    pivot_longer(tmx10:tmx9) %>%  
                    mutate(
                        tmx_max_month = name[which(value==max(value))] %>% 
                            str_split_i('tmx', 2) %>% 
                            as.numeric() %>% 
                            mean() %>% 
                            as.factor(),
                        tmx_min_month = name[which(value==min(value))] %>% 
                            str_split_i('tmx', 2) %>% 
                            as.numeric() %>% 
                            mean() %>% 
                            as.factor()
                    ) %>% 
                    pivot_wider(names_from = name, values_from = value) %>% 
                    pivot_longer(ppt10:ppt9) %>% 
                    mutate(
                        ppt_max_month = name[which(value==max(value))] %>% 
                            str_split_i('ppt', 2) %>% 
                            as.numeric() %>% 
                            mean() %>% 
                            as.factor(),
                        ppt_min_month = name[which(value==min(value))] %>% 
                            str_split_i('ppt', 2) %>% 
                            as.numeric() %>% 
                            mean() %>% 
                            as.factor()
                    ) %>% 
                    pivot_wider(names_from = name, values_from = value)
            }, .progress=T
        ) %>% 
        bind_rows() %>% 
        select(-c(tmx10:tmx9, spp_id)) %>% 
        ungroup() %>% 
        distinct()
    
    #Create full SWD set
    write_csv(spp_extract_binded, out.filename)
}

trainingTestingFunc = function(sp, in.filename, training.filename, testing.filename){
    
    combined = read_csv(in.filename)
    
    #Create training and testing sets
    set.seed(123)
    split = initial_split(combined, prop = 0.7, strata = "presence")
    training = training(split)
    testing = testing(split)
    write_csv(training, training.filename)
    write_csv(testing, testing.filename)
}
