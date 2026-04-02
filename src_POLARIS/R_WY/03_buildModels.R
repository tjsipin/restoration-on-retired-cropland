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
dir.create('data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/', recursive = T)
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

map(
    names,
    function(sp){
        set.seed(123)
        enm_eval.filename = paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/', sp, "_tuned_args_lowFilter.rds")
        enm_eval_res.filename = paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/', sp, "_tuned_args_res_lowFilter.rds")
        
        # if(file.exists(enm_eval.filename)){
        #     print(paste0("enm eval exists: ", sp))
        #     enm_eval = readRDS(enm_eval.filename)
        #     enm_eval_res = eval.results(enm_eval)
        #     saveRDS(enm_eval_res, enm_eval_res.filename)
        #     return(NULL)
        # }
        cat("Processing:", sp, "\n")
        swd <- read_csv(paste0('data_v5_50km/3_swd/wy/', '/training_', sp, '_soil200cm_lowFilter_wy.csv')) %>% 
            select(
                x, y, month, year, 
                aet, tmx, tmx_summer_mean, tdiff, 
                ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
                ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
                tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
                bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
                presence
            ) %>% 
            mutate(
                salinity=as.factor(salinity),
                
            ) 
        swdOcc <- swd %>% filter(presence == 1)
        swdBack <- swd %>% filter(presence == 0)
        
        ### occ and background can only be x/y
        occ <- swdOcc %>% 
            dplyr::select(x, y)
        
        bg <- swdBack %>% 
            dplyr::select(x, y)
        
        block <- get.block(occ, bg, orientation = "lon_lat")
        ## check for even number in each group
        print(table(block$occs.grp))
        
        #   1   2   3   4 
        # 503 502 502 502 
        
        # Evaluating --------------------------------------
        occs.z <- swdOcc %>% 
            dplyr::select(!c(month, year, presence))
        
        bg.z <- swdBack %>% 
            dplyr::select(!c(month, year, presence))
        
        
        #ENMeval fitting
        enm_eval <- ENMevaluate(occs=occs.z, bg=bg.z, algorithm = "maxnet",
                                tune.args = list(fc = c("LQP"),
                                                 rm = seq(0.5, 4, by=0.5)), 
                                partitions = "block", partition.settings = list(orientation='lat_lon'))
        saveRDS(enm_eval, enm_eval.filename)
        
        enm_eval_res = eval.results(enm_eval)
        saveRDS(enm_eval_res, enm_eval_res.filename)
        return(NULL)
    }
)

#Select best model arguments
best_model_args = map(
    names,
    function(sp){
        res.filename = paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/', sp, "_tuned_args_res_lowFilter.rds")
        args.filename = paste0("data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/", sp, "_finalModelArgs_lowFilter.rds")
        res = readRDS(res.filename) %>%
            select(rm, auc.train, auc.val.avg, auc.diff.avg, or.10p.avg, AICc) %>%
            rowwise() %>%
            mutate(
                # Following the Radosavljevic paper:
                score = (1-auc.val.avg)*auc.diff.avg*(or.10p.avg)
            ) %>%
            tibble()
        best_res = res %>%
            filter(score == min(score, na.rm=T))
            # filter(or.10p.avg == min(or.10p.avg, na.rm=T))
        args <- c(
            'jackknife=TRUE',
            'autofeature=TRUE',
            'responsecurves=TRUE',
            paste0('betamultiplier=', best_res$rm),
            # 'betamultiplier=2.5',
            'linear=TRUE',
            'quadratic=TRUE',
            'product=TRUE',
            'threshold=FALSE',
            'hinge=FALSE',
            'maximumiterations=100000',
            'writeplotdata=TRUE')
        saveRDS(args, args.filename)

        res %>%
            select(rm, auc.train, auc.val.avg, or.10p.avg, auc.diff.avg, AICc, score) %>%
            pivot_longer(-rm) %>%
            ggplot() +
            geom_line(
                aes(x = rm, y = value)
            ) +
            facet_wrap(~name, scales='free')
    }
)

best_model_args_tibble = map(
    names,
    function(sp){
        res.filename = paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/', sp, "_tuned_args_res_lowFilter.rds")
        res = readRDS(res.filename) %>%
            select(rm, auc.train, auc.val.avg, auc.diff.avg, or.10p.avg, AICc) %>%
            rowwise() %>%
            mutate(
                species = sp,
                # Following the Radosavljevic paper:
                score = (1-auc.val.avg)*auc.diff.avg*(or.10p.avg)
            ) %>%
            tibble()
        return(res)
    }
) %>%
    bind_rows()

best_model_args_tibble %>%
    pivot_longer(-c(species, rm)) %>%
    ggplot() +
    geom_line(aes(x = rm, y = value, color = species)) +
    facet_wrap(~name, scales='free')

#Create training model for validation
map(
    names,
    function(sp){
        set.seed(123)
        print(sp)
        # Read in all points
        training <- read_csv(
            paste0("data_v5_50km/3_swd/wy/", "/training_", sp, "_soil200cm_lowFilter_wy.csv")
        ) %>%
            mutate(
                salinity=as.factor(salinity)
            ) %>%
            select(
                x, y, month, year,
                aet, tmx, tmx_summer_mean, tdiff, 
                ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
                ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
                tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
                bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
                presence
            )

        ## Select env predictors for model (presence + background)
        x <- training %>%
            mutate(
                salinity=as.factor(salinity)
            ) %>%
            # dplyr::select(aet:ph)
            select(
                aet, tmx, tmx_summer_mean, tdiff, 
                ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
                ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
                tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
                bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity
            )


        ## Specify occurrence data
        p <- training$presence


        #reading in best args
        args = readRDS(paste0("data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/", sp, "_finalModelArgs_lowFilter.rds"))

        ## Path to save results
        out.path <- paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/', sp, '/lowFilter/model_training/')
        dir.create(out.path, recursive=T)

        ## Final Model Creation
        model <- maxent(x=x, p=p, path=out.path, args=args)

        ## Save model for pred rasters
        saveRDS(model, file = paste0(out.path, "/", sp, "_training_sdm.rds"))
    }
)

#Get AUCs
aucs = map(
    1:length(names),
    function(i){
        training = read_csv(paste0('data_v5_50km/3_swd/wy/training_', names[i], '_soil200cm_lowFilter_wy.csv'))
        testing = read_csv(paste0('data_v5_50km/3_swd/wy/testing_', names[i], '_soil200cm_lowFilter_wy.csv'))
        model = readRDS(paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/', names[i], '/lowFilter/model_training/', names[i], '_training_sdm.rds'))
        best_rm = readRDS(paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/', names[i], '_finalModelArgs_lowFilter.rds'))[4] %>%
            str_split_i('=', 2)

        training_pred = training %>%
            select(
                x, y, month, year,
                aet, tmx, tmx_summer_mean, tdiff, 
                ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
                ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
                tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
                bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
                presence
            ) %>%
            mutate(
                salinity=as.factor(salinity)
            ) %>%
            filter(complete.cases(.)) %>%
            as.data.frame() %>% 
            mutate(pred = predict(model, .))

        testing_pred = testing %>%
            select(
                x, y, month, year,
                aet, tmx, tmx_summer_mean, tdiff, 
                ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
                ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
                tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
                bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
                presence
            ) %>%
            mutate(
                salinity=as.factor(salinity)
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

aucs

# aucs_monthly = map(
#     1:length(names),
#     function(i){
#         training = read_csv(paste0('data_v5_50km/3_swd/wy/training_', names[i], '_soil200cm_lowFilter_wy.csv'))
#         testing = read_csv(paste0('data_v5_50km/3_swd/wy/testing_', names[i], '_soil200cm_lowFilter_wy.csv'))
#         model = readRDS(paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/', names[i], '/lowFilter/model_training/', names[i], '_training_sdm.rds'))
#         best_rm = readRDS(paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/', names[i], '_finalModelArgs_lowFilter.rds'))[4] %>%
#             str_split_i('=', 2)
# 
#         training_pred = training %>%
#             select(
#                 x, y, month, year,
#                 aet, tmx, tmx_summer_mean, tdiff, 
#                 ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
#                 ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
#                 tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
#                 bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
#                 presence
#             ) %>%
#             mutate(
#                 salinity=as.factor(salinity)
#             ) %>%
#             filter(complete.cases(.)) %>%
#             mutate(pred = predict(model, .))
# 
#         testing_pred = testing %>%
#             select(
#                 x, y, month, year,
#                 aet, tmx, tmx_summer_mean, tdiff, 
#                 ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
#                 ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
#                 tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
#                 bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
#                 presence
#             ) %>%
#             mutate(
#                 salinity=as.factor(salinity)
#             ) %>%
#             filter(complete.cases(.)) %>%
#             mutate(pred = predict(model, .))
# 
#         training.auc = pROC::auc(training_pred$presence, training_pred$pred) %>%
#             as.numeric()
#         testing.auc = pROC::auc(testing_pred$presence, testing_pred$pred) %>%
#             as.numeric()
# 
#         out = tibble(
#             sp = names[i],
#             training.auc = training.auc,
#             testing.auc = testing.auc,
#             reg_mult = best_rm
#         )
#         return(out)
#     }
# ) %>% bind_rows()

variable_importance = map(
    names,
    function(sp){
        model = readRDS(paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
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
        model = readRDS(paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy/', sp, '/lowFilter/model_training/', sp, '_training_sdm.rds'))
        training = read_csv(paste0('data_v5_50km/3_swd/wy/training_', sp, '_soil200cm_lowFilter_wy.csv'))
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
        training = read_csv(paste0('data_v5_50km/3_swd/wy/training_', sp, '_soil200cm_lowFilter_wy.csv'))
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
            paste0("data_v5_50km/3_swd/wy/", "/training_", sp, "_soil200cm_lowFilter_wy.csv")
        ) %>%
            mutate(
                salinity=as.factor(salinity)
            ) %>%
            select(
                x, y, month, year,
                aet, tmx, tmx_summer_mean, tdiff, 
                ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
                ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
                tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
                bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
                presence
            )
        testing <- read_csv(
            paste0("data_v5_50km/3_swd/wy/", "/testing_", sp, "_soil200cm_lowFilter_wy.csv")
        ) %>%
            mutate(
                salinity=as.factor(salinity)
            ) %>%
            select(
                x, y, month, year,
                aet, tmx, tmx_summer_mean, tdiff, 
                ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
                ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
                tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
                bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
                presence
            )

        full = training %>%
            rbind(testing)

        ## Select env predictors for model (presence + background)
        x <- full %>%
            # dplyr::select(aet:ph)
            select(
                aet, tmx, tmx_summer_mean, tdiff, 
                ppt10, ppt11, ppt12, ppt1, ppt2, ppt3,
                ppt4, ppt5, ppt6, ppt7, ppt8, ppt9,
                tmx_max_month, tmx_min_month, ppt_max_month, ppt_min_month,
                bd_POLARIS, clay_POLARIS, om_POLARIS, ph_POLARIS, silt_POLARIS, salinity,
            ) %>%
            mutate(
                salinity=as.factor(salinity)
            )


        ## Specify occurrence data
        p <- full$presence


        #reading in best args
        args = readRDS(paste0("data_v5_50km_POLARIS/4_maxent_outputs/wy/tuning/", sp, "_finalModelArgs_lowFilter.rds"))

        ## Path to save results
        out.path <- paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy//', sp, '/lowFilter/model/')
        dir.create(out.path, recursive=T)

        ## Final Model Creation
        model <- maxent(x=x, p=p, path=out.path, args=args)

        ## Save model for pred rasters
        saveRDS(model, file = paste0(out.path, "/", sp, "_final_sdm.rds"))
    }
)


map(
    names,
    function(sp){
        out.path <- paste0('data_v5_50km_POLARIS/4_maxent_outputs/wy//', sp, '/lowFilter/model/')
        file.exists(paste0(out.path, "/", sp, "_final_sdm.rds"))
    }
)