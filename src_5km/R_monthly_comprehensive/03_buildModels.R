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
dir.create('data_v2/4_maxent_outputs/comprehensive/tuning/', recursive = T)
set.seed(123)

## sp to run
names = c(
    'a_menziesii',
    'a_intermedia',
    'c_pungens',
    'l_pentachaeta',
    'p_arborea',
    'p_ciliata',
    'a_polycarpa',
    'c_lasiophyllus',
    'l_californica',
    'l_gracilis'
)

map(
    names,
    function(sp){
        enm_eval.filename = paste0('data_v2/4_maxent_outputs/comprehensive/tuning/', sp, "_tuned_args_lowFilter.rds")
        enm_eval_res.filename = paste0('data_v2/4_maxent_outputs/comprehensive/tuning/', sp, "_tuned_args_res_lowFilter.rds")
        
        # if(file.exists(enm_eval.filename)){
        #     print(paste0("enm eval exists: ", sp))
        #     enm_eval = readRDS(enm_eval.filename)
        #     enm_eval_res = eval.results(enm_eval)
        #     saveRDS(enm_eval_res, enm_eval_res.filename)
        #     return(NULL)
        # }
        cat("Processing:", sp, "\n")
        swd <- read_csv(paste0('data_v2/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')) %>% 
            mutate(tdiff = tmx - tmn) %>% 
            select(
                lon, lat, month, year, 
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            ) %>% 
            mutate(drclass = as.factor(drclass))
        
        swdOcc <- swd %>% filter(presence == 1)
        swdBack <- swd %>% filter(presence == 0)
        
        ### occ and background can only be x/y
        occ <- swdOcc %>% 
            dplyr::select(lon, lat)
        
        bg <- swdBack %>% 
            dplyr::select(lon, lat)
        
        block <- get.block(occ, bg, orientation = "lon_lat")
        ## check for even number in each group
        print(table(block$occs.grp))
        
        # Evaluating --------------------------------------
        occs.z <- swdOcc %>% 
            dplyr::select(!c(month, year, presence))
        
        bg.z <- swdBack %>% 
            dplyr::select(!c(month, year, presence))
        
        
        #ENMeval fitting
        enm_eval <- ENMevaluate(occs.z, bg=bg.z, algorithm = "maxnet",
                                tune.args = list(fc = c("LQP"),
                                                 rm = seq(0.5, 4, by=0.5)), 
                                partitions = "block")
        saveRDS(enm_eval, enm_eval.filename)
        
        enm_eval_res = eval.results(enm_eval)
        saveRDS(enm_eval_res, enm_eval_res.filename)
        return(NULL)
    }
)

map(
    names,
    function(sp){
        res.filename = paste0('data_v2/4_maxent_outputs/comprehensive/tuning/', sp, "_tuned_args_res_lowFilter.rds")
        args.filename = paste0("data_v2/4_maxent_outputs/comprehensive/tuning/", sp, "_finalModelArgs_lowFilter.rds")
        res = readRDS(res.filename)
        best_res = res %>%
            select(rm, auc.train, auc.val.avg, auc.diff.avg, or.10p.avg) %>% 
            rowwise() %>% 
            mutate(
                # Following the Radosavljevic paper:
                score = (1-auc.val.avg)*auc.diff.avg*(or.10p.avg)
            ) %>% 
            tibble() %>% 
            filter(score == min(score))
        args <- c(
            'jackknife=TRUE', 
            'autofeature=TRUE', 
            'responsecurves=TRUE', 
            paste0('betamultiplier=', best_res$rm),
            'linear=TRUE',
            'quadratic=TRUE',
            'product=TRUE',
            'threshold=FALSE',
            'hinge=FALSE', 
            'maximumiterations=100000', 
            'writeplotdata=TRUE')
        saveRDS(args, args.filename)
        
        res %>% 
            select(rm, or.10p.avg, auc.val.avg, auc.diff.avg, AICc) %>% 
            pivot_longer(-rm) %>% 
            ggplot() + 
            geom_line(
                aes(x = rm, y = value) 
            ) +
            facet_wrap(~name, scales='free')
    }
) 

map(
    names,
    function(sp){
        set.seed(123)
        print(sp)    
        # Read in all points
        training <- read_csv(paste0('data_v2/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')) %>% 
            mutate(tdiff = tmx - tmn) %>% 
            select(
                lon, lat, month, year, 
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            ) %>% 
            mutate(drclass = as.factor(drclass))
        
        ## Select env predictors for model (presence + background)
        x <- training %>% 
            # dplyr::select(aet:ph) 
            select(aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff, cec, drclass, om, ph, salinity)
        
        ## Specify occurrence data
        p <- training$presence
        
        
        #reading in best args
        args = readRDS(paste0("data_v2/4_maxent_outputs/comprehensive/tuning/", sp, "_finalModelArgs_lowFilter.rds"))
        
        ## Path to save results
        out.path <- paste0('data_v2/4_maxent_outputs/comprehensive/', sp, '/lowFilter/model_training/')
        dir.create(out.path, recursive=T)
        
        ## Final Model Creation
        model <- maxent(x=x, p=p, path=out.path, args=args)
        
        ## Save model for pred rasters
        saveRDS(model, file = paste0(out.path, "/", sp, "_training_sdm.rds"))
    }
)

#Get AUCs
aucs = map(
    names,
    function(sp){
        training = read_csv(paste0('data_v2/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv'))
        testing = read_csv(paste0('data_v2/3_swd/monthly/comprehensive/testing_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv'))
        model = readRDS(paste0('data_v2/4_maxent_outputs/comprehensive/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
        best_rm = readRDS(paste0('data_v2/4_maxent_outputs/comprehensive/tuning/', sp, '_finalModelArgs_lowFilter.rds'))[4] %>% 
            str_split_i('=', 2)
        
        training_pred = training %>% 
            mutate(tdiff = tmx - tmn) %>% 
            select( 
                lon, lat, month, year, 
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            ) %>% 
            filter(complete.cases(.)) %>% 
            mutate(pred = predict(model, .))
        
        testing_pred = testing %>% 
            mutate(tdiff = tmx - tmn) %>% 
            select(
                lon, lat, month, year, 
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
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

aucs_monthly = map(
    names,
    function(sp){
        training = read_csv(paste0('data_v2/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv'))
        testing = read_csv(paste0('data_v2/3_swd/monthly/comprehensive/testing_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv'))
        model = readRDS(paste0('data_v2/4_maxent_outputs/comprehensive/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
        best_rm = readRDS(paste0('data_v2/4_maxent_outputs/comprehensive/tuning/', sp, '_finalModelArgs_lowFilter.rds'))[4] %>% 
            str_split_i('=', 2)
        
        training_pred = training %>% 
            mutate(tdiff = tmx - tmn) %>% 
            select( 
                lon, lat, month, year, 
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            ) %>% 
            filter(complete.cases(.)) %>% 
            mutate(pred = predict(model, .)) %>% 
            group_by(month) %>% 
            summarize(AUC = ifelse(length(unique(presence))==1, NA, pROC::auc(presence, pred) %>% as.numeric())) %>% 
            mutate(partition = "Training")
        
        testing_pred = testing %>% 
            mutate(tdiff = tmx - tmn) %>% 
            select( 
                lon, lat, month, year, 
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            ) %>% 
            filter(complete.cases(.)) %>% 
            mutate(pred = predict(model, .)) %>% 
            group_by(month) %>% 
            summarize(AUC = ifelse(length(unique(presence))==1, NA, pROC::auc(presence, pred) %>% as.numeric())) %>% 
            mutate(partition = "Testing")
        
        out = training_pred %>% 
            full_join(testing_pred) %>% 
            mutate(
                sp = sp,
                reg_mult = best_rm
            ) %>% 
            pivot_wider(id_cols = c(sp, reg_mult, partition), names_from = month, names_prefix = "Month_", values_from = AUC)
        
        return(out)
    }
) %>% 
    bind_rows()

variable_importance = map(
    names,
    function(sp){
        model = readRDS(paste0('data_v2/4_maxent_outputs/comprehensive/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
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
        model = readRDS(paste0('data_v2/4_maxent_outputs/comprehensive/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
        training = read_csv(paste0('data_v2/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')) %>% 
            mutate(tdiff=tmx-tmn)
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
        training = read_csv(paste0('data_v2/3_swd/monthly/comprehensive/training_', sp, '_soil200cm_lowFilter_monthly_comprehensive.csv')) 
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
            paste0("data_v2/3_swd/monthly/comprehensive/training_", sp, "_soil200cm_lowFilter_monthly_comprehensive.csv")
        ) %>%
            mutate(tdiff = tmx - tmn) %>% 
            mutate(drclass = as.factor(drclass)) %>%
            select(
                lon, lat, month, year,
                cwd, aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            )
        testing <- read_csv(
            paste0("data_v2/3_swd/monthly/comprehensive/testing_", sp, "_soil200cm_lowFilter_monthly_comprehensive.csv")
        ) %>%
            mutate(tdiff = tmx - tmn) %>% 
            mutate(drclass = as.factor(drclass)) %>%
            select(
                lon, lat, month, year,
                cwd, aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
            )
        
        full = training %>%
            rbind(testing) 
        
        ## Select env predictors for model (presence + background)
        x <- full %>%
            select(aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff, cec, drclass, om, ph, salinity)
        
        
        ## Specify occurrence data
        p <- full$presence
        
        
        #reading in best args
        args = readRDS(paste0("data_v2/4_maxent_outputs/comprehensive/tuning/", sp, "_finalModelArgs_lowFilter.rds"))
        
        ## Path to save results
        out.path <- paste0('data_v2/4_maxent_outputs/comprehensive/', sp, '/lowFilter/model/')
        dir.create(out.path, recursive=T)
        
        ## Final Model Creation
        model <- maxent(x=x, p=p, path=out.path, args=args)
        
        ## Save model for pred rasters
        saveRDS(model, file = paste0(out.path, "/", sp, "_final_sdm.rds"))
    }
)