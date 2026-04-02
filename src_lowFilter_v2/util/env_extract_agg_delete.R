#This version, we will be extracting the species-specific current (or previous) ecologically relevant period.

setwd("/home/tjsipin/will_tj")

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

relevantAgg = function(
        points, this_year, 
        bcm_filename_df,
        input_path='data/0_env/bcm/bcmv8_historic/2000_2023_monthly/', 
        relevant_months,
        target_crs, spp
){
    points_y = points %>% 
        filter(start_year == this_year)
    
    #DF of input filenames
    input_filenames = bcm_filename_df %>% 
        filter(
            start_year == this_year
        ) %>% 
        mutate(filename = paste0(input_path, "/", filename))
    #Subset to eco-relevant months
    input_filenames_relevant = input_filenames %>% 
        filter(
            start_year == this_year,
            month %in% relevant_months
        ) 
    if(nrow(input_filenames_relevant)==0) return(NULL)
    
    #For each variable, pull appropriate file names 
    aet_filenames = input_filenames_relevant %>% 
        filter(var == "aet") %>% 
        pull(filename)
    cwd_filenames = input_filenames_relevant %>% 
        filter(var == "cwd") %>% 
        pull(filename)
    ppt_filenames = input_filenames_relevant %>% 
        filter(var == "ppt") %>% 
        pull(filename)
    ppt_winter_filenames = input_filenames %>% 
        filter(var == "ppt") %>% 
        #Filter to months defined as winter
        filter(month_num %in% c(12, 1, 2)) %>% 
        pull(filename)
    tmx_filenames = input_filenames_relevant %>% 
        filter(var == "tmx") %>% 
        pull(filename)
    tmn_filenames = input_filenames_relevant %>% 
        filter(var == "tmn") %>% 
        pull(filename)
    tmx_summer_filenames = input_filenames %>% 
        filter(var == "tmx") %>% 
        #Filter to months defined as summer
        filter(month_num %in% 6:8) %>% 
        pull(filename)
    
    if(spp=='p_ciliata'){
        tmx_summer_filenames = bcm_filename_df %>% 
            filter(
                start_year == this_year - 1
            ) %>% 
            mutate(filename = paste0(input_path, "/", filename)) %>% 
            filter(var == "tmx") %>% 
            #Filter to months defined as summer
            filter(month_num %in% 6:8) %>% 
            pull(filename)
    }
    
    
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
    
    rast_filenames = list(
        aet_filenames,
        cwd_filenames, 
        ppt_filenames,
        ppt_winter_filenames,
        tmx_filenames,
        tmn_filenames,
        tmx_summer_filenames
    )
    
    for(f in rast_filenames){
        if(length(f)==0){
            message("Length of files for at least one BCM variable is 0.")
            return(NULL)
        } 
    }
    
    #Rasterize all sets of file names
    aet_agg = rast(aet_filenames) %>% 
        mean() %>% 
        crop(CA)
    cwd_agg = rast(cwd_filenames) %>% 
        sum() %>% 
        crop(CA)
    ppt_agg = rast(ppt_filenames) %>% 
        sum() %>% 
        crop(CA)
    ppt_winter_agg = rast(ppt_filenames) %>% 
        sum() %>% 
        crop(CA)
    tmx_agg = rast(tmx_filenames) %>% 
        mean() %>% 
        crop(CA)
    tmn_agg = rast(tmn_filenames) %>% 
        mean() %>% 
        crop(CA)
    tmx_summer_agg = rast(tmx_summer_filenames) %>% 
        mean() %>% 
        crop(CA)
    tdiff_agg = tmx_agg - tmn_agg
    
    ph = rast(ph_filename) %>% 
        crop(aet_agg)
    om = rast(om_filename) %>% 
        crop(aet_agg)
    cec = rast(cec_filename) %>% 
        crop(aet_agg)
    drclass = rast(drclass_filename) %>% 
        crop(aet_agg)
    
    bd_POLARIS = rast(bd_POLARIS_filename) %>% 
        crop(aet_agg)
    clay_POLARIS = rast(clay_POLARIS_filename) %>% 
        crop(aet_agg)
    om_POLARIS = rast(om_POLARIS_filename) %>% 
        crop(aet_agg)
    ph_POLARIS = rast(ph_POLARIS_filename) %>% 
        crop(aet_agg)
    sand_POLARIS = rast(sand_POLARIS_filename) %>% 
        crop(aet_agg)
    silt_POLARIS = rast(silt_POLARIS_filename) %>% 
        crop(aet_agg)
    
    #Combine all rasters
    all_agg = rast(list(
        aet_agg, cwd_agg, ppt_agg, ppt_winter_agg, 
        tmx_agg, tmn_agg, tmx_summer_agg, tdiff_agg, 
        ph, om, cec, drclass,
        bd_POLARIS, clay_POLARIS, om_POLARIS,
        ph_POLARIS, sand_POLARIS, silt_POLARIS
    )) 
    #Rename the combined raster layers
    names(all_agg) = c(
        "aet", "cwd", "ppt", "ppt_winter_sum", 
        "tmx", "tmn", "tmx_summer_mean", "tdiff", 
        "ph", "om", "cec", "drclass",
        "bd_POLARIS", "clay_POLARIS", "om_POLARIS", 
        "ph_POLARIS", "sand_POLARIS", "silt_POLARIS"
    )
    
    #Vectorize the points
    points_y_vect = points_y %>% 
        vect(geom=c("lon", "lat"), "EPSG:4326") %>% 
        project(target_crs)
    
    #Extract the combined raster values at each point
    extracted = terra::extract(all_agg, points_y_vect, method = "simple", bind=T, xy=T) %>% 
        as.data.frame()
    return(extracted)
}

extractEnvAgg = function(spp, occ_back=c('occ', 'back'), relevant_months){
    if(occ_back=='occ'){
        in.dir = paste0('data/1_occ/combined_spp_occ/')
        in.filename = paste0(in.dir, spp, '_lowFilter.csv')
        presence = 1
    } else{
        in.dir = paste0('data/2_background/')
        in.filename = paste0(in.dir, 'back_', spp, '_5km_lowFilter_even.csv')
        presence = 0
    }
    
    out.dir = 'data/3_swd/agg'
    
    points.filename = paste0(in.filename)
    out.filename = paste0('data/3_swd/agg/', 'swd_', spp, '_', occ_back, '_soil200cm_lowFilter_agg.csv')
    # if(file.exists(out.filename)) return(NULL)
    
    #Getting order of relevant months
    month_df = tibble(
        mon = month.abb %>% str_to_lower(),
        mm = 1:12
    ) %>% 
        arrange(mm) %>% 
        filter(
            (mon %in% relevant_months) | (mm %in% relevant_months)
        ) %>% 
        #jump = december to january for example
        mutate(
            jump = case_when(
                mm + 1 != mm[row_number() + 1] ~ T,
                .default = F
            )
        ) %>% 
        mutate(row = 1:n())
    
    if(sum(month_df$jump)==0){
        pre_jump_months = month_df$mm
    } else{
        pre_jump_months = month_df %>% 
            filter(!(row %in% 1:row[jump==T][1])) %>% 
            pull(mm)
        post_jump_months = month_df %>% 
            filter(row %in% 1:row[jump==T][1]) %>% 
            pull(mm)
    }
    
    #1) Read in the occurrence or background data 
    spp_occBack = read_csv(in.filename) %>% 
        select(lon, lat, month, year) %>% 
        mutate(
            presence = presence, 
            start_year = ifelse(month %in% pre_jump_months, year, year - 1)
        ) %>% 
        mutate(spp_id = row_number())
    
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
        ) %>% 
        #Species specific:
        mutate(
            # start_year = ifelse(month_num %in% pre_jump_months, year, year - 1)
            start_year = case_when(
                month_num %in% pre_jump_months ~ year,
                !(month_num %in% pre_jump_months) ~ year - 1
            )
        )
    
    spp_extract = map(
        .x = sort(unique(spp_occBack$start_year)),
        ~ relevantAgg(points = spp_occBack, this_year = .x,
                      bcm_filename_df = bcm_filename_df,
                      relevant_months = month_df$mon,
                      target_crs = target_crs, spp = spp),
        .progress = T
    ) 
    spp_extract_binded = spp_extract %>% 
        bind_rows() 
    
    #Create full SWD set
    write_csv(spp_extract_binded, out.filename)
}

trainingTestingFunc = function(sp, occ.filename, back.filename, training.filename, testing.filename){
    occ = read_csv(occ.filename)
    back = read_csv(back.filename)
    
    combined = occ %>% 
        rbind(back)
    
    #Create training and testing sets
    set.seed(123)
    split = initial_split(combined, prop = 0.7, strata = "presence")
    training = training(split)
    testing = testing(split)
    write_csv(training, training.filename)
    write_csv(testing, testing.filename)
}
