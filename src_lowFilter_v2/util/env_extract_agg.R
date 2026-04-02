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
        target_crs, sp
){
    if(this_year < min(bcm_filename_df$ecorelevant_year)) {
        message("Outside year bounds")
        return(NULL)
    }
    points_y = points %>% 
        filter(ecorelevant_year == this_year)
    #Gather files
    bcm_filenames = list.files(
        "data/0_env/bcm/bcmv8_historic/2000_2023_monthly/",
        pattern = ".tif"
    ) 
    
    
    # 3:18PM ------------------------------------------------------------------
    
    bcm_tibble = bcm_filename_df %>%
        mutate(filename = paste0(input_path, filename)) %>% 
        arrange(var, year, month_num) %>%
        group_split(var) %>%
        map(
            function(df){
                df_row = df %>%
                    mutate(
                        row = row_number(),
                        max_row = case_when(
                            (ecorelevant_year == this_year) & (month_num == relevant_months[length(relevant_months)]) ~ T,
                            .default = F
                        )
                    ) 
                if(sum(df_row$max_row)!=1) return(NULL)
                df_row %>%
                    filter(
                        row %in% seq(row[max_row] - 11, row[max_row])
                    ) %>%
                    ungroup() %>%
                    select(-row, -max_row)
            }
        ) %>%
        bind_rows() 
    if(nrow(bcm_tibble)==0) return(NULL)
    
    bcm_tibble_relevant = bcm_tibble %>% 
        filter(month_num %in% relevant_months)
    
    
    
    if(nrow(bcm_tibble_relevant)==0) return(NULL)
    
    #For each variable, pull appropriate file names 
    aet_filenames = bcm_tibble_relevant %>% 
        filter(var == "aet") %>% 
        pull(filename)
    cwd_filenames = bcm_tibble_relevant %>% 
        filter(var == "cwd") %>% 
        pull(filename)
    ppt_filenames = bcm_tibble_relevant %>% 
        filter(var == "ppt") %>% 
        pull(filename)
    ppt_winter_filenames = bcm_tibble %>% 
        filter(var == "ppt") %>% 
        #Filter to months defined as winter
        filter(month_num %in% c(12, 1, 2)) %>% 
        pull(filename)
    tmx_filenames = bcm_tibble_relevant %>% 
        filter(var == "tmx") %>% 
        pull(filename)
    tmn_filenames = bcm_tibble_relevant %>% 
        filter(var == "tmn") %>% 
        pull(filename)
    tmx_summer_filenames = bcm_tibble %>% 
        filter(var == "tmx") %>% 
        filter(
            month_num %in% 6:8
        ) %>% 
        #make sure summer temperature is for the first year in the ecorelevant period
        filter(
            year==min(year)
        ) %>% 
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
        resample(aet_agg)
    om = rast(om_filename) %>% 
        resample(aet_agg)
    cec = rast(cec_filename) %>% 
        resample(aet_agg)
    drclass = rast(drclass_filename) %>% 
        resample(aet_agg, method='mode')
    
    bd_POLARIS = rast(bd_POLARIS_filename) %>% 
        resample(aet_agg)
    clay_POLARIS = rast(clay_POLARIS_filename) %>% 
        resample(aet_agg)
    om_POLARIS = rast(om_POLARIS_filename) %>% 
        resample(aet_agg)
    ph_POLARIS = rast(ph_POLARIS_filename) %>% 
        resample(aet_agg)
    sand_POLARIS = rast(sand_POLARIS_filename) %>% 
        resample(aet_agg)
    silt_POLARIS = rast(silt_POLARIS_filename) %>% 
        resample(aet_agg)
    
    salinity = rast(salinity_filename) %>% 
        resample(aet_agg, method='mode')
    
    #Combine all rasters
    all_agg = rast(list(
        aet_agg, cwd_agg, ppt_agg, ppt_winter_agg, 
        tmx_agg, tmn_agg, tmx_summer_agg, tdiff_agg, 
        ph, om, cec, drclass,
        bd_POLARIS, clay_POLARIS, om_POLARIS,
        ph_POLARIS, sand_POLARIS, silt_POLARIS,
        salinity
    )) 
    #Rename the combined raster layers
    names(all_agg) = c(
        "aet", "cwd", "ppt", "ppt_winter_sum", 
        "tmx", "tmn", "tmx_summer_mean", "tdiff", 
        "ph", "om", "cec", "drclass",
        "bd_POLARIS", "clay_POLARIS", "om_POLARIS", 
        "ph_POLARIS", "sand_POLARIS", "silt_POLARIS",
        "salinity"
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

extractEnvAgg = function(sp, relevant_months){
    occ_filename = paste0('data/1_occ/combined_spp_occ/', sp, '_lowFilter.csv')
    occ = read_csv(occ_filename) %>% 
        mutate(sp = sp, presence = 1) %>% 
        select(sp, lon, lat, month, year, presence)
    
    back_filename = paste0('data/2_background/back_', sp, '_5km_lowFilter_even.csv')
    back = read_csv(back_filename) %>% 
        mutate(sp = sp, presence = 0) %>% 
        select(sp, lon, lat, month, year, presence)
    
    out.dir = 'data/3_swd/agg'
    
    out.filename = paste0('data/3_swd/agg/', 'swd_', sp, '_soil200cm_lowFilter_agg.csv')
    # if(file.exists(out.filename)) return(NULL)
    
    #Getting order of months
    jump = !identical(min(relevant_months):max(relevant_months), relevant_months)
    relevant_month_df = tibble(
        mm = relevant_months,
        mon = str_to_lower(month.abb)[relevant_months],
        ecorelevant=T
    )
    not_relevant_month_df = tibble(
        mm = (1:12)[!1:12 %in% relevant_months]
    ) %>% 
        mutate(
            mon = str_to_lower(month.abb)[mm],
            ecorelevant=F
        )
    month_df = relevant_month_df %>% 
        rbind(not_relevant_month_df) %>% 
        mutate(jump = jump) %>% 
        # arrange(mm) %>% 
        mutate(row = row_number()) %>%
        mutate(
            adjust_year = case_when(
                jump & (mm >= relevant_months[1]) ~ 0,
                jump & (mm < relevant_months[1]) ~ -1,
                !jump & (mm < relevant_months[1]) ~ -1,
                !jump & (mm > relevant_months[length(relevant_months)]) ~ 0,
                !jump & (mm %in% relevant_months) ~ 0
            )
        ) %>% 
        select(-row, -ecorelevant, -jump)
    
    #1) Read in the occurrence or background data 
    points = occ %>% 
        rbind(back) %>% 
        mutate(spp_id = row_number()) %>% 
        full_join(month_df, by=join_by(month==mm)) %>%
        mutate(
            ecorelevant_year = year + adjust_year
        )
    
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
        full_join(month_df, by=join_by(month_num==mm)) %>%
        mutate(
            ecorelevant_year = year + adjust_year
        )
    
    spp_extract = map(
        .x = sort(unique(points$ecorelevant_year)),
        function(.x){
            message(.x)
            relevantAgg(points = points, this_year = .x,
                        bcm_filename_df = bcm_filename_df,
                        relevant_months = relevant_months,
                        target_crs = target_crs, sp = sp)
        },
        .progress = T
    ) 
    spp_extract_binded = spp_extract %>% 
        bind_rows() 
    
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
