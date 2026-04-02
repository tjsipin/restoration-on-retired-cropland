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
library(paletteer)
library(scico)
library(RColorBrewer)
library(tidyterra)
library(stringr)

tidymodels_prefer()

central_valley = vect('data/central_valley/ds2632.gdb') %>% 
    project('epsg:3310')
CA = tigris::states() %>% 
    vect() %>% 
    project("epsg:3310") %>% 
    filter(STUSPS=="CA")

## sp to run
names = c(
    'a_menziesii',
    'c_pungens',
    'l_pentachaeta',
    'p_arborea',
    'p_ciliata',
    'a_polycarpa',
    'c_lasiophyllus'
)
## list of months
months = month.abb %>% str_to_lower()




# Create color palette
pl <- rev(scico(13, palette = 'bam'))

# Set color for each class
pal<-c(
    "0"=pl[1],
    "1"=pl[2],
    "2"=pl[3],
    "3"=pl[4],
    "4"=pl[5],
    "5"=pl[6],
    "6"=pl[7],
    "7"=pl[8],
    "8"=pl[9],
    "9"=pl[10],
    "10"=pl[11],
    "11"=pl[12],
    "12"=pl[13]
)

# Color for na values
na_col<-'grey50'
# Create custom theme
custom_theme<-theme_void()+
    theme(plot.background = element_rect(fill="#202020",color=NA))

#Code to plot thresholds per species per model-years, after summing by month
p10Visualize = function(sp, model){
    print(paste(sp, model, sep = " | "))
    sp_proper = str_replace(sp, "_", ". ") %>% 
        str_to_sentence()
    
    out.dir = paste0("data/5_figs/lowFilter/comprehensive/p10_sum/")
    dir.create(out.dir, recursive=T)
    out.png.filename.CA = paste0("data/5_figs/lowFilter/comprehensive/p10_sum/", sp, "_", model, "_lowFilter_sum_CA.png")
    out.png.filename.CV = paste0("data/5_figs/lowFilter/comprehensive/p10_sum/", sp, "_", model, "_lowFilter_sum_CV.png")
    in.rast.filename.CA = paste0("data/4_maxent_outputs/comprehensive/", sp, "/", "lowFilter/p10/p10_", sp, "_", model, "_monthly_CA.tif")
    in.rast.filename.CV = paste0("data/4_maxent_outputs/comprehensive/", sp, "/", "lowFilter/p10/p10_", sp, "_", model, "_monthly_CV.tif")
    in.rast.CA = rast(in.rast.filename.CA) %>% 
        mutate(sum = factor(sum, levels = as.integer(0:12)))
    in.rast.CV = rast(in.rast.filename.CV) %>% 
        mutate(sum = factor(sum, levels = as.integer(0:12)))
    
    ggplot() +
        geom_spatraster(data = in.rast.CA) + 
        scale_fill_manual(values=pal, na.value=na_col, drop=F)+
        theme_void() +
        labs(title = paste(sp_proper, model, sep = " — "), subtitle = paste0("Low filter") %>% str_to_sentence)
    ggsave(out.png.filename.CA, scale = 1)
    
    ggplot() +
        geom_spatraster(data = in.rast.CV) + 
        scale_fill_manual(values=pal, na.value=na_col, drop=F)+
        theme_void() +
        labs(title = paste(sp_proper, model, sep = " — "), subtitle = paste0("Low filter") %>% str_to_sentence)
    ggsave(out.png.filename.CV, scale = 1)
}

map(
    .x = names,
    .f = ~ p10Visualize(sp = .x, model = "2000_2023")
)
map(
    .x = names,
    .f = ~ p10Visualize(sp = .x, model = "MIROC45_2070_2099")
)
map(
    .x = names,
    .f = ~ p10Visualize(sp = .x, model = "MIROC85_2070_2099")
)




# Fig. 5: Sum of suitable pixels per month --------------------------------

speciesPerMonthHist = function(sp, crop_poly = c("CA",  "CV")){
    print(sp)
    input.filenames.hist = list.files(
        paste0("data/4_maxent_outputs/comprehensive/", sp, "/lowFilter/p10/p10_hist"),
        full.names=T
    )
    input.filenames.MIROC45 = list.files(
        paste0("data/4_maxent_outputs/comprehensive/", sp, "/lowFilter/p10/p10_MIROC45"),
        full.names=T
    )
    input.filenames.MIROC85 = list.files(
        paste0("data/4_maxent_outputs/comprehensive/", sp, "/lowFilter/p10/p10_MIROC85"),
        full.names=T
    )
    
    if(crop_poly == "CA"){
        input.rast.hist = rast(input.filenames.hist) %>%
            crop(CA, mask=T)
        input.rast.MIROC45 = rast(input.filenames.MIROC45) %>%
            crop(CA, mask=T)
        input.rast.MIROC85 = rast(input.filenames.MIROC85) %>%
            crop(CA, mask=T)
    } else if(crop_poly == "CV"){
        input.rast.hist = rast(input.filenames.hist) %>%
            crop(central_valley, mask=T)
        input.rast.MIROC45 = rast(input.filenames.MIROC45) %>%
            crop(central_valley, mask=T)
        input.rast.MIROC85 = rast(input.filenames.MIROC85) %>%
            crop(central_valley, mask=T)
    }
    
    
    ha_per_month.hist = map(
        1:12,
        function(i){
            print(i)
            rast = input.rast.hist[[i]]
            mon = names(rast) %>%
                str_sub(1, 3)
            #Get proportion of suitability
            num_ha = terra::global(rast, fun='sum', na.rm=T)$sum*0.027 # convert pixel to hectare
            
            out = tibble(
                mon=mon,
                mm = which(mon==months),
                # presence = rep(1, num_ha),
                num_ha = num_ha,
                model = "2000 to 2023"
            )
            return(out)
        }
    ) %>%
        bind_rows() %>%
        arrange(mm)
    ha_per_month.MIROC45 = map(
        1:12,
        function(i){
            print(i)
            rast = input.rast.MIROC45[[i]]
            mon = names(rast) %>%
                str_sub(1, 3)
            #Get proportion of suitability
            num_ha = terra::global(rast, fun='sum', na.rm=T)$sum*0.027 # convert pixel to hectare
            
            out = tibble(
                mon=mon,
                mm = which(mon==months),
                # presence = rep(1, num_ha),
                num_ha = num_ha,
                model = "2070 to 2099 (MIROC RCP 4.5)"
            )
            return(out)
        }
    ) %>%
        bind_rows() %>%
        arrange(mm)
    ha_per_month.MIROC85 = map(
        1:12,
        function(i){
            print(i)
            rast = input.rast.MIROC85[[i]]
            mon = names(rast) %>%
                str_sub(1, 3)
            #Get proportion of suitability
            num_ha = terra::global(rast, fun='sum', na.rm=T)$sum*0.027 # convert pixel to hectare
            
            out = tibble(
                mon=mon,
                mm = which(mon==months),
                # presence = rep(1, num_ha),
                num_ha = num_ha,
                model = "2070 to 2099 (MIROC RCP 8.5)"
            )
            return(out)
        }
    ) %>%
        bind_rows() %>%
        arrange(mm)
    
    ha_per_month.all = ha_per_month.hist %>%
        rbind(ha_per_month.MIROC45) %>%
        rbind(ha_per_month.MIROC85) %>%
        mutate(mon = ordered(str_to_sentence(mon), str_to_sentence(months)))
    
    sp_pretty = sp %>%
        str_replace('_', '. ') %>%
        str_to_sentence()
    crop_poly_pretty = ifelse(crop_poly=="CA", "CA", "the Great Valley")
    
    if(crop_poly=="CA"){
        ggplot(data=ha_per_month.all) +
            geom_histogram(
                aes(
                    x = mon,
                    y = num_ha,
                    fill = model,
                    group = model
                ),
                color = 'darkgrey',
                stat='identity',
                position = 'dodge',
                width = 0.7
            ) +
            theme_classic() +
            scale_fill_brewer(palette = 'Set1') +
            scale_y_continuous(labels = scales::comma, limits = c(0, 101000)) +
            # scale_y_continuous(labels = scales::comma) +
            labs(
                y = 'Hectares',
                fill = "Years (and Climate Scenario)",
                title = paste0('Suitability across ', crop_poly_pretty, ' (hectares)'),
                subtitle = sp_pretty
            ) +
            theme(
                text = element_text(
                    family="Times New Roman", size = 18
                ), 
                plot.subtitle = element_text(
                    face = 'italic', size = 14
                )
            )
    } else{
        ggplot(data=ha_per_month.all) +
            geom_histogram(
                aes(
                    x = mon,
                    y = num_ha,
                    fill = model,
                    group = model
                ),
                color = 'darkgrey',
                stat='identity',
                position = 'dodge',
                width = 0.7
            ) +
            theme_classic() +
            scale_fill_brewer(palette = 'Set1') +
            scale_y_continuous(labels = scales::comma, limits = c(0, 18000)) +
            # scale_y_continuous(labels = scales::comma) +
            labs(
                y = 'Hectares',
                fill = "Years (and Climate Scenario)",
                title = paste0('Suitability across ', crop_poly_pretty, ' (hectares)'),
                subtitle = sp_pretty
            ) +
            theme(
                text = element_text(
                    family="Times New Roman", size = 18
                ), 
                plot.subtitle = element_text(
                    face = 'italic', size = 14
                )
            )
    }
    output.filename = paste0('data/5_figs/lowFilter/comprehensive/ha_histograms/p10_', sp, '_monthly_barplots_', crop_poly, '.png')
    ggsave(output.filename, scale = 2, width = 7, height = 4, units = 'in', bg = 'white')
    
    ha_per_month.all_wider = ha_per_month.all %>% 
        pivot_wider(id_cols = 'model', names_from = c('mon'), values_from = 'num_ha') %>% 
        mutate(species = sp_pretty, .before = model)
    
    output.csv.filename = paste0('data/6_tables/lowFilter/comprehensive/ha_histograms/p10_', sp, '_suitable_ha.csv')
    write_csv(ha_per_month.all_wider, output.csv.filename)
}

dir.create('data/5_figs/lowFilter/comprehensive/ha_histograms/', recursive = T)
dir.create('data/6_tables/lowFilter/comprehensive/ha_histograms/', recursive = T)
map(
    .x = names,
    .f = ~speciesPerMonthHist(sp = .x, crop_poly = "CA"),
    .progress = T
)
map(
    .x = names,
    .f = ~speciesPerMonthHist(sp = .x, crop_poly = "CV")
)


predictionAggregate = function(sp, monthly_dist_path){
    in.path = paste0("data/4_maxent_outputs/comprehensive/",
                      sp, 
                      "/lowFilter/", monthly_dist_path)
    in.training.filename = paste0('data/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')
    in.files = list.files(
        in.path,
        full.names=T
    )
    in.rast = rast(in.files) %>% 
        mean() %>% 
        crop(central_valley, mask = T)
    in.training = read_csv(in.training.filename) %>% 
        select(lat, lon, presence) %>% 
        filter(presence == 1)
    
    training_vect = vect(in.training, geom=c('lon', 'lat'), 'epsg:4326') %>% 
        project(in.rast) %>% 
        crop(in.rast)
    
    sp_pretty = sp %>% 
        str_replace("_", ". ") %>% 
        str_to_sentence()
    if(str_detect(monthly_dist_path, "hist")){
        model_years_pretty = "Current conditions"
    } else if(str_detect(monthly_dist_path, "45")){
        model_years_pretty = "End-of-century conditions (MIROC RCP 4.5)"
    } else if(str_detect(monthly_dist_path, "85")){
        model_years_pretty = "End-of-century conditions (MIROC RCP 8.5)"
    }
    out.filename = paste0("data/4_maxent_outputs/comprehensive/figs/prediction_aggregate/", monthly_dist_path, "/", sp, '.png')
    dir.create(paste0("data/4_maxent_outputs/comprehensive/figs/prediction_aggregate/", monthly_dist_path))
    ggplot() + 
        geom_spatraster(
            data = in.rast
        ) +
        geom_spatvector(
            data = training_vect,
            color = 'green3',
            shape = '+',
            size = 4
        ) +
        labs(
            title = paste0("Mean annual suitability for ", sp_pretty),
            subtitle = model_years_pretty
        ) +
        scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        theme_void()
    
    
    ggsave(out.filename, bg='white')
}

predictionPerMonth = function(sp, monthly_dist_path){
    out.filename = paste0('data/5_figs/lowFilter/comprehensive/monthly_suitability/', monthly_dist_path, '/', sp, '_monthly_suitability.png')
    in.path = paste0("data/4_maxent_outputs/comprehensive//",
                     sp, 
                     "/lowFilter/", monthly_dist_path)
    in.swd.filename = paste0('data/3_swd/monthly/comprehensive/swd_', sp, '_occ_soil200cm_lowFilter_monthly_comprehensive.csv')
    in.files = list.files(
        in.path,
        full.names=T
    )
    in.rast = rast(in.files) %>% 
        crop(central_valley, mask = T)
    names(in.rast) = names(in.rast) %>% 
        str_split_i('_', 1) %>% 
        str_to_sentence()
    in.swd = read_csv(in.swd.filename) %>% 
        select(lat, lon) 
    swd_vect = vect(in.swd, geom=c('lon', 'lat'), 'epsg:4326') %>% 
        project(in.rast) 
    swd_vect_filter = terra::extract(in.rast, swd_vect, xy=T) %>% 
        select(x, y, everything()) %>% 
        vect(geom=c('x', 'y'), crs=crs(in.rast)) %>% 
        select(-ID) %>% 
        filter(!is.na(Apr))

        
    
    sp_pretty = sp %>% 
        str_replace("_", ". ") %>% 
        str_to_sentence()
    if(str_detect(monthly_dist_path, "hist")){
        model_years_pretty = "Current conditions"
    } else if(str_detect(monthly_dist_path, "45")){
        model_years_pretty = "End-of-century conditions (MIROC RCP 4.5)"
    } else if(str_detect(monthly_dist_path, "85")){
        model_years_pretty = "End-of-century conditions (MIROC RCP 8.5)"
    }
    ggplot() + 
        geom_spatraster(
            data = in.rast
        ) +
        geom_spatvector(
            data = swd_vect_filter,
            color = 'green3',
            shape = '+',
            alpha = 0.5,
            size = 4
        ) +
        labs(
            title = paste0("Mean suitability for ", sp_pretty),
            subtitle = model_years_pretty
        ) +
        scale_fill_gradientn(colours = c('navy', 'yellow', 'red'), limits=c(0, 1), na.value='transparent') +
        theme_void() + 
        facet_wrap(~lyr)
    
    ggsave(out.filename, background = )
}

dir.create("data/5_figs/lowFilter/comprehensive/monthly_suitability/monthly_dist_hist", recursive = T)
dir.create("data/5_figs/lowFilter/comprehensive/monthly_suitability/monthly_dist_MIROC45", recursive = T)
dir.create("data/5_figs/lowFilter/comprehensive/monthly_suitability/monthly_dist_MIROC85", recursive = T)

predictionPerMonth(sp = 'a_polycarpa', monthly_dist_path = 'monthly_dist_hist')
predictionPerMonth(sp = 'p_arborea', monthly_dist_path = 'monthly_dist_hist')
predictionPerMonth(sp = 'c_pungens', monthly_dist_path = 'monthly_dist_hist')
predictionPerMonth(sp = 'l_pentachaeta', monthly_dist_path = 'monthly_dist_hist')
predictionPerMonth(sp = 'p_ciliata', monthly_dist_path = 'monthly_dist_hist')
predictionPerMonth(sp = 'a_menziesii', monthly_dist_path = 'monthly_dist_hist')
predictionPerMonth(sp = 'c_lasiophyllus', monthly_dist_path = 'monthly_dist_hist')

predictionPerMonth(sp = 'a_polycarpa', monthly_dist_path = 'monthly_dist_MIROC45')
predictionPerMonth(sp = 'p_arborea', monthly_dist_path = 'monthly_dist_MIROC45')
predictionPerMonth(sp = 'c_pungens', monthly_dist_path = 'monthly_dist_MIROC45')
predictionPerMonth(sp = 'l_pentachaeta', monthly_dist_path = 'monthly_dist_MIROC45')
predictionPerMonth(sp = 'p_ciliata', monthly_dist_path = 'monthly_dist_MIROC45')
predictionPerMonth(sp = 'a_menziesii', monthly_dist_path = 'monthly_dist_MIROC45')
predictionPerMonth(sp = 'c_lasiophyllus', monthly_dist_path = 'monthly_dist_MIROC45')

predictionPerMonth(sp = 'a_polycarpa', monthly_dist_path = 'monthly_dist_MIROC85')
predictionPerMonth(sp = 'p_arborea', monthly_dist_path = 'monthly_dist_MIROC85')
predictionPerMonth(sp = 'c_pungens', monthly_dist_path = 'monthly_dist_MIROC85')
predictionPerMonth(sp = 'l_pentachaeta', monthly_dist_path = 'monthly_dist_MIROC85')
predictionPerMonth(sp = 'p_ciliata', monthly_dist_path = 'monthly_dist_MIROC85')
predictionPerMonth(sp = 'a_menziesii', monthly_dist_path = 'monthly_dist_MIROC85')
predictionPerMonth(sp = 'c_lasiophyllus', monthly_dist_path = 'monthly_dist_MIROC85')

envDistribution = function(sp){
    output.filepath = paste0('data/5_figs/lowFilter/comprehensive/env_distribution/monthly_comprehensive/')
    dir.create(output.filepath, recursive=T)
    output.filename = paste0(output.filepath, sp, '_envDistribution_.png')
    training = read_csv(paste0('data/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')) %>% 
        mutate(
            tdiff = tmx - tmn,
            presence = as.factor(presence)
        )
    training_long = training %>% 
        select(presence, aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff, cec, drclass, om, ph) %>% 
        rename(
            AET=aet, 
            `Max temperature`=tmx, 
            `Winter Precipitation`=ppt_winter_sum,
            `Max summer temperature`=tmx_summer_mean,
            `Temperature difference`=tdiff,
            `CEC`=cec,
            `Drainage class`=drclass,
            `Organic matter`=om,
            `pH`=ph
        ) %>% 
        pivot_longer(-presence) %>% 
        mutate(presence = ifelse(presence==1, 'Occurrence', 'Background'))
    
    pretty_sp = sp %>% 
        str_replace('_', '. ') %>% 
        str_to_sentence()
    
    ggplot(training_long) +
        geom_density(
            aes(
                x = value,
                group = presence,
                color = presence
            )
        ) + 
        facet_wrap(~name, scales='free') +
        theme_classic() +
        labs(
            title = "Comprehensive monthly model environmental density distributions",
            subtitle = pretty_sp,
            x = "",
            y = '',
            color = ''
        )
    ggsave(output.filename, scale=2)
}

map(names, envDistribution)


# AUCs table --------------------------------------------------------------

#Get AUCs
aucs = map(
    names,
    function(sp){
        training = read_csv(paste0('data/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv'))
        testing = read_csv(paste0('data/3_swd/monthly/comprehensive/testing_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv'))
        model = readRDS(paste0('data/4_maxent_outputs/comprehensive/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
        best_rm = readRDS(paste0('data/4_maxent_outputs/comprehensive/tuning/', sp, '_finalModelArgs_lowFilter.rds'))[4] %>% 
            str_split_i('=', 2)
        
        training_pred = training %>% 
            mutate(tdiff = tmx - tmn) %>% 
            select( 
                lon, lat, month, year, 
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph,
                presence
            ) %>% 
            filter(complete.cases(.)) %>% 
            mutate(pred = predict(model, .))
        
        testing_pred = testing %>% 
            mutate(tdiff = tmx - tmn) %>% 
            select(
                lon, lat, month, year, 
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph,
                presence
            ) %>% 
            filter(complete.cases(.)) %>% 
            mutate(pred = predict(model, .))
        
        training.auc = pROC::auc(training_pred$presence, training_pred$pred) %>% 
            as.numeric()
        testing.auc = pROC::auc(testing_pred$presence, testing_pred$pred) %>% 
            as.numeric()
        
        out = tibble(
            sp = sp,
            training.auc = training.auc,
            testing.auc = testing.auc,
            reg_mult = best_rm
        )
        return(out)
    }
) %>% 
    bind_rows()

write_csv(aucs, "data/6_tables/lowFilter/comprehensive/aucs.csv")

# Permutation importance --------------------------------------------------

getVImp = function(sp){
    mod.filename = paste0('data/4_maxent_outputs/comprehensive/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    sp_pretty = sp %>% 
        str_replace("_", ". ") %>% 
        str_to_sentence()
    mod = readRDS(mod.filename)
    permutation.imp = tibble(
        variable = rownames(mod@results),
        val = mod@results[,1]
    ) %>% 
        filter(str_detect(variable, 'permutation')) %>% 
        mutate(
            sp = sp_pretty,
            variable = str_split_i(variable, '.permutation', 1)
        ) %>% 
        pivot_wider(id_cols = sp, names_from = variable, values_from = val) %>% 
        select(sp, aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff, cec, drclass, om, ph)
    return(permutation.imp)
}

var_imps = map(
    names,
    getVImp
) %>% 
    bind_rows()
write_csv(var_imps, 'data/6_tables/lowFilter/comprehensive/perm_var_imps.csv')
