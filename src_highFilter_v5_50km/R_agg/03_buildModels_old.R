library(tidyverse)    ## always
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
library(SDMtune)
tidymodels_prefer()
dir.create('data_preEllen/4_maxent_outputs/agg/tuning/', recursive = T)
set.seed(123)

names = c(
    'a_polycarpa',
    'p_arborea',
    'c_pungens',
    'l_pentachaeta',
    'p_ciliata',
    'a_menziesii',
    'a_intermedia',
    'c_lasiophyllus',
    'l_californica',
    'l_gracilis'
)

#Get AUCs
old_aucs = map(
    1:length(names),
    function(i){
        if(!file.exists(paste0('data_preEllen/3_swd/agg/training_', names[i], '_soil200cm_lowFilter_agg.csv'))) return(NULL)
        training = read_csv(paste0('data_preEllen/3_swd/agg/training_', names[i], '_soil200cm_lowFilter_agg.csv'))
        testing = read_csv(paste0('data_preEllen/3_swd/agg/testing_', names[i], '_soil200cm_lowFilter_agg.csv'))
        model = readRDS(paste0('data_preEllen/4_maxent_outputs/agg/', names[i], '/lowFilter/model_training/', names[i], '_training_sdm.rds'))
        best_rm = readRDS(paste0('data_preEllen/4_maxent_outputs/agg/tuning/', names[i], '_finalModelArgs_lowFilter.rds'))[4] %>%
            str_split_i('=', 2)
        
        old_training_pred = training %>%
            select(
                x, y, month, year,
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, #salinity,
                presence
            ) %>%
            mutate(
                drclass=as.factor(drclass),
                # salinity=as.factor(salinity)
            ) %>%
            filter(complete.cases(.)) %>%
            as.data.frame() %>% 
            mutate(pred = predict(model, .))
        
        old_testing_pred = testing %>%
            select(
                x, y, month, year,
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, #salinity,
                presence
            ) %>%
            mutate(
                drclass=as.factor(drclass),
                # salinity=as.factor(salinity)
            ) %>%
            filter(complete.cases(.)) %>%
            as.data.frame() %>% 
            mutate(pred = predict(model, .))
        
        training.auc = pROC::auc(training_pred$presence, training_pred$pred) %>%
            as.numeric()
        testing.auc = pROC::auc(testing_pred$presence, testing_pred$pred) %>%
            as.numeric()
        
        out = tibble(
            sp = names[i],
            training.auc = training.auc,
            testing.auc = testing.auc,
            reg_mult = best_rm
        )
        return(out)
    }
) %>%
    bind_rows()


aucs_monthly = map(
    1:length(names),
    function(i){
        training = read_csv(paste0('data_preEllen/3_swd/agg/training_', names[i], '_soil200cm_lowFilter_agg.csv'))
        testing = read_csv(paste0('data_preEllen/3_swd/agg/testing_', names[i], '_soil200cm_lowFilter_agg.csv'))
        model = readRDS(paste0('data_preEllen/4_maxent_outputs/agg/', names[i], '/lowFilter/model_training/', names[i], '_training_sdm.rds'))
        best_rm = readRDS(paste0('data_preEllen/4_maxent_outputs/agg/tuning/', names[i], '_finalModelArgs_lowFilter.rds'))[4] %>%
            str_split_i('=', 2)
        
        training_pred = training %>%
            select(
                x, y, month, year,
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            ) %>%
            mutate(
                drclass=as.factor(drclass),
                salinity=as.factor(salinity)
            ) %>%
            filter(complete.cases(.)) %>%
            mutate(pred = predict(model, .))
        
        testing_pred = testing %>%
            select(
                x, y, month, year,
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            ) %>%
            mutate(
                drclass=as.factor(drclass),
                salinity=as.factor(salinity)
            ) %>%
            filter(complete.cases(.)) %>%
            mutate(pred = predict(model, .))
        
        training.auc = pROC::auc(training_pred$presence, training_pred$pred) %>%
            as.numeric()
        testing.auc = pROC::auc(testing_pred$presence, testing_pred$pred) %>%
            as.numeric()
        
        out = tibble(
            sp = names[i],
            training.auc = training.auc,
            testing.auc = testing.auc,
            reg_mult = best_rm
        )
        return(out)
    }
)
aucs_monthly %>%
    bind_rows()

variable_importance = map(
    names,
    function(sp){
        model = readRDS(paste0('data_preEllen/4_maxent_outputs/agg/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
        var_imp = model@results %>%
            as.data.frame() %>%
            mutate(variable = rownames(.)) %>%
            tibble() %>%
            filter(str_detect(variable, 'permutation.importance')) %>%
            mutate(
                variable = str_replace(variable, '.permutation.importance', ''),
                sp = sp
            ) %>%
            rename(importance = V1) %>%
            select(sp, variable, importance) %>%
            pivot_wider(names_from = variable, values_from = importance)
        return(var_imp)
    }
) %>%
    bind_rows()

pearson_correlation = map(
    names,
    function(sp){
        model = readRDS(paste0('data_preEllen/4_maxent_outputs/agg/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
        training = read_csv(paste0('data_preEllen/3_swd/agg/training_', sp, '_soil200cm_lowFilter_agg.csv'))
        p = training %>%
            filter(presence==1)
        a = training %>%
            filter(presence==0)
        cor = evaluate(model=model, p=p, a=a)@cor
        out = tibble(
            sp = sp,
            cor = cor
        )
    }
) %>%
    bind_rows()

presence_distribution = map(
    names,
    function(sp){
        training = read_csv(paste0('data_preEllen/3_swd/agg/training_', sp, '_soil200cm_lowFilter_agg.csv'))
        presence_absence = training %>%
            mutate(sp = sp) %>%
            group_by(sp, presence) %>%
            summarize(freq = n()/nrow(.)) %>%
            ungroup() %>%
            pivot_wider(names_from=presence, values_from=freq, names_prefix = 'presence')
    }
) %>%
    bind_rows()


#FINAL and IMPORTANT: Train final model using both training and testing data
map(
    names,
    function(sp){
        set.seed(123)
        print(sp)
        # Read in all points
        training <- read_csv(
            paste0("data_preEllen/3_swd/agg/", "/training_", sp, "_soil200cm_lowFilter_agg.csv")
        ) %>%
            mutate(
                drclass=as.factor(drclass),
                salinity=as.factor(salinity)
            ) %>%
            select(
                x, y, month, year,
                cwd, aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            )
        testing <- read_csv(
            paste0("data_preEllen/3_swd/agg/", "/testing_", sp, "_soil200cm_lowFilter_agg.csv")
        ) %>%
            mutate(
                drclass=as.factor(drclass),
                salinity=as.factor(salinity)
            ) %>%
            select(
                x, y, month, year,
                cwd, aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            )
        
        full = training %>%
            rbind(testing)
        
        ## Select env predictors for model (presence + background)
        x <- full %>%
            # dplyr::select(aet:ph)
            select(aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff, cec, drclass, om, ph, salinity) %>%
            mutate(
                drclass=as.factor(drclass),
                salinity=as.factor(salinity)
            )
        
        
        ## Specify occurrence data
        p <- full$presence
        
        
        #reading in best args
        args = readRDS(paste0("data_preEllen/4_maxent_outputs/agg/tuning/", sp, "_finalModelArgs_lowFilter.rds"))
        
        ## Path to save results
        out.path <- paste0('data_preEllen/4_maxent_outputs/agg//', sp, '/lowFilter/model/')
        dir.create(out.path, recursive=T)
        
        ## Final Model Creation
        model <- maxent(x=x, p=p, path=out.path, args=args)
        
        ## Save model for pred rasters
        saveRDS(model, file = paste0(out.path, "/", sp, "_final_sdm.rds"))
    }
)
