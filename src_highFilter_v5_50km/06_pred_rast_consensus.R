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
library(ggnewscale) #citation("ggnewscale")
library(cowplot)
library(ggspatial)
library(ggpubr)
library(kableExtra) #tables
library(paletteer) #tables pretty
library(SDMtune) #VIP
library(wacolors)
tidymodels_prefer()

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

consensusAUC = function(sp){
    er_mod_filename = paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    wy_mod_filename = paste0('data_v2/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    
    er_testing_filename = paste0('data_v2/3_swd/agg/testing_', sp, '_soil200cm_lowFilter_agg.csv')
    wy_testing_filename = paste0('data_v2/3_swd/wy/testing_', sp, '_soil200cm_lowFilter_wy.csv')
    
    er_mod = readRDS(er_mod_filename)
    wy_mod = readRDS(wy_mod_filename)
    
    er_testing = read_csv(er_testing_filename) %>% 
        mutate(xy=paste0(x,', ', y))
    wy_testing = read_csv(wy_testing_filename) %>% 
        mutate(xy=paste0(x,', ', y))
    
    common_xy = tibble(
        xy = base::intersect(er_testing$xy, wy_testing$xy)
    ) %>% 
        mutate(pred_id = row_number())
    
    er_testing_xy = er_testing %>% 
        right_join(common_xy)
    wy_testing_xy = wy_testing %>% 
        right_join(common_xy)
    
    er_testing_pred = er_testing_xy %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(er_mod, .)) %>% 
        select(pred_id, pred, presence)
    wy_testing_pred = wy_testing_xy %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(wy_mod, .)) %>% 
        select(pred_id, pred, presence)
    
    common_pred = er_testing_pred %>% 
        rbind(wy_testing_pred) %>% 
        group_by(pred_id, presence) %>% 
        filter(n() == 2) %>% 
        summarize(pred = mean(pred)) %>% 
        ungroup() 
    
    er_auc = pROC::auc(er_testing_pred$presence, er_testing_pred$pred) %>%
        as.numeric()
    wy_auc = pROC::auc(wy_testing_pred$presence, wy_testing_pred$pred) %>%
        as.numeric()
    common_auc = pROC::auc(common_pred$presence, common_pred$pred) %>%
        as.numeric()
    sp_pretty = sp %>% 
        str_replace('_', '. ') %>% 
        str_to_sentence()
    out = tibble(
        Species = sp_pretty,
        er_auc = er_auc,
        wy_auc = wy_auc,
        consensus_auc = common_auc
    )
}
consensusFOR = function(sp){
    er_mod_filename = paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    wy_mod_filename = paste0('data_v2/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    
    er_training_filename = paste0('data_v2/3_swd/agg/training_', sp, '_soil200cm_lowFilter_agg.csv')
    wy_training_filename = paste0('data_v2/3_swd/wy/training_', sp, '_soil200cm_lowFilter_wy.csv')
    er_testing_filename = paste0('data_v2/3_swd/agg/testing_', sp, '_soil200cm_lowFilter_agg.csv')
    wy_testing_filename = paste0('data_v2/3_swd/wy/testing_', sp, '_soil200cm_lowFilter_wy.csv')
    
    er_mod = readRDS(er_mod_filename)
    wy_mod = readRDS(wy_mod_filename)
    
    er_training = read_csv(er_training_filename) %>% 
        mutate(xy=paste0(x,', ', y))
    wy_training = read_csv(wy_training_filename) %>% 
        mutate(xy=paste0(x,', ', y))
    er_testing = read_csv(er_testing_filename) %>% 
        mutate(xy=paste0(x,', ', y))
    wy_testing = read_csv(wy_testing_filename) %>% 
        mutate(xy=paste0(x,', ', y))
    
    common_training_xy = tibble(
        xy = base::intersect(er_training$xy, wy_training$xy)
    ) %>% 
        mutate(pred_id = row_number())
    common_testing_xy = tibble(
        xy = base::intersect(er_testing$xy, wy_testing$xy)
    ) %>% 
        mutate(pred_id = row_number())
    
    er_training_xy = er_training %>% 
        right_join(common_training_xy)
    wy_training_xy = wy_training %>% 
        right_join(common_training_xy)
    
    er_testing_xy = er_testing %>% 
        right_join(common_testing_xy)
    wy_testing_xy = wy_testing %>% 
        right_join(common_testing_xy)
    
    er_training_pred = er_training_xy %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(er_mod, .)) %>% 
        select(pred_id, pred, presence)
    wy_training_pred = wy_training_xy %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(wy_mod, .)) %>% 
        select(pred_id, pred, presence)
    common_training_pred = er_training_pred %>% 
        rbind(wy_training_pred) %>% 
        group_by(pred_id, presence) %>% 
        filter(n() == 2) %>% 
        summarize(pred = mean(pred)) %>% 
        ungroup() 
    
    er_testing_pred = er_testing_xy %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(er_mod, .)) %>% 
        select(pred_id, pred, presence)
    wy_testing_pred = wy_testing_xy %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(wy_mod, .)) %>% 
        select(pred_id, pred, presence)
    common_testing_pred = er_testing_pred %>% 
        rbind(wy_testing_pred) %>% 
        group_by(pred_id, presence) %>% 
        filter(n() == 2) %>% 
        summarize(pred = mean(pred)) %>% 
        ungroup() 
    
    er_p10_thresh = er_training_pred %>% 
        filter(presence==1) %>% 
        pull(pred) %>% 
        quantile(.1)
    wy_p10_thresh = wy_training_pred %>% 
        filter(presence==1) %>% 
        pull(pred) %>% 
        quantile(.1)
    common_p10_thresh = common_training_pred %>% 
        filter(presence==1) %>% 
        pull(pred) %>% 
        quantile(.1)
    
    er_testing_pred_p10 = er_testing_pred %>% 
        mutate(
            pred_thresh = ifelse(pred > er_p10_thresh, 1, 0) %>% 
                as.factor(),
            presence = presence %>% 
                as.factor()
        )
    wy_testing_pred_p10 = wy_testing_pred %>% 
        mutate(
            pred_thresh = ifelse(pred > wy_p10_thresh, 1, 0) %>% 
                as.factor(),
            presence = presence %>% 
                as.factor()
        )
    common_testing_pred_p10 = common_testing_pred %>% 
        mutate(
            pred_thresh = ifelse(pred > common_p10_thresh, 1, 0) %>% 
                as.factor(),
            presence = presence %>% 
                as.factor()
        )
    
    er_conf = conf_mat(er_testing_pred_p10, truth=presence, estimate=pred_thresh)
    er_FOR = er_conf$table[1,2]/(er_conf$table[1,1] + er_conf$table[1,2])
    wy_conf = conf_mat(wy_testing_pred_p10, truth=presence, estimate=pred_thresh)
    wy_FOR = wy_conf$table[1,2]/(wy_conf$table[1,1] + wy_conf$table[1,2])
    common_conf = conf_mat(common_testing_pred_p10, truth=presence, estimate=pred_thresh)
    common_FOR = common_conf$table[1,2]/(common_conf$table[1,1] + common_conf$table[1,2])
    
    sp_pretty = sp %>% 
        str_replace('_', '. ') %>% 
        str_to_sentence()
    out = tibble(
        Species = sp_pretty,
        er_FOR = er_FOR,
        wy_FOR = wy_FOR,
        consensus_FOR = common_FOR
    )
}

aucs_v2 = map(names, consensusAUC) %>% 
    bind_rows() %>% 
    rename(`Ecorelevant AUC`=er_auc, `Water year AUC`=wy_auc, `Consensus AUC`=consensus_auc) %>% 
    pivot_longer(-Species)
ggplot(aucs_v2) +
    geom_bar(aes(x=Species, y=value, fill=name, group=name), stat='identity', position='dodge') +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    theme_light() +
    labs(y='AUC', fill='')

FORs_v2 = map(names, consensusFOR) %>% 
    bind_rows() %>% 
    rename(`Ecorelevant FOR`=er_FOR, `Water year FOR`=wy_FOR, `Consensus FOR`=consensus_FOR) %>% 
    pivot_longer(-Species)
ggplot(FORs_v2) +
    geom_bar(aes(x=Species, y=value, fill=name, group=name), stat='identity', position='dodge') +
    scale_y_continuous(breaks = seq(0, 0.1, 0.005)) +
    theme_light() +
    labs(y='FOR', fill='')

aucs = aucs_v2 %>% 
    mutate(model = 'v2') %>% 
    rbind(
        aucs_cvBG %>% 
            mutate(model='CV BG')
    )
ggplot(aucs) +
    geom_bar(aes(x=model, y=value, fill=name, group=name), stat='identity', position='dodge') +
    scale_y_continuous(breaks = seq(0, 1, 0.05)) +
    theme_light() +
    labs(y='auc', fill='') +
    facet_wrap(~Species, nrow = 1)

FORs = FORs_v2 %>% 
    mutate(model = 'v2') %>% 
    rbind(
        FORs_cvBG %>% 
            mutate(model='CV BG')
    )
ggplot(FORs) +
    geom_bar(aes(x=model, y=value, fill=name, group=name), stat='identity', position='dodge') +
    scale_y_continuous(breaks = seq(0, 0.1, 0.005)) +
    theme_light() +
    labs(y='FOR', fill='') +
    facet_wrap(~Species, nrow = 1)

distributionSpPlotWithoutLegendER = function(sp){
    print(sp)
    output.filename = paste0("data_v2/5_figs/lowFilter/agg/fig2_", sp, "_2000_2023.png")
    input.filename = paste0("data_v2/4_maxent_outputs/agg/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_agg.tif")
    sppOccs.filename = paste0("data/1_occ/combined_spp_occ/", sp, "_lowFilter.csv")
    sppOccs = read_csv(sppOccs.filename)
    sppOccs_vect = sppOccs %>%
        vect(geom=c('lon', 'lat'), crs='epsg:4326')
    
    input.rast = rast(input.filename)
    rast_CA = input.rast %>%
        crop(CA, mask=T)
    rast_CV = input.rast %>%
        crop(central_valley, mask=T)
    
    vect_CV = rast_CV %>%
        as.polygons()
    sp_pretty = sp %>%
        str_replace("_", ". ") %>%
        str_to_sentence()
    
    # plot_CV = ggplot() +
    #     geom_spatraster(data=rast_CV) +
    #     geom_spatvector(data=sppOccs_vect, shape='plus') +
    #     # scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
    #     scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
    #     labs(
    #         fill = "Predicted suitability",
    #         title = sp_pretty
    #     ) +
    #     theme_void() +
    #     theme(
    #         # legend.position = 'bottom',
    #         legend.position = 'none',
    #         legend.direction = 'horizontal',
    #         legend.key.width = unit(1, 'cm'),
    #         legend.title.position = 'top',
    #         text = element_text("Times", face = 'italic')
    #     )
    #
    # CV = ifel(!is.na(rast_CV), 0, NA)
    
    plot_CA = ggplot() +
        geom_spatraster(data=rast_CA) +
        geom_spatvector(data=sppOccs_vect, shape=6, color='green', fill = 'transparent', size=1, alpha=0.2) +
        scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        new_scale_fill() +
        #Greyed out CV
        # geom_spatvector(data=vect_CV, na.rm=T, fill='transparent') +
        # scale_fill_gradientn(colours = c('grey20', 'grey21'), na.value='transparent') +
        theme_void() +
        labs(title = sp_pretty) +
        theme(legend.position = "none")
    # plot_inset = ggdraw(plot = plot_CV) +
    #     draw_plot(plot_CA, x = 0.48, y = 0.54, width = 0.4, height = 0.4)
    # plot_inset
    
    
    # ggsave(output.filename,
    #        scale=1.25, bg = 'white', units = 'px', width = 1080, height = 1920)
    
    return(plot_CA)
}

Fig_2_ER = map(names, distributionSpPlotWithoutLegendER)
cowplot::plot_grid(plotlist = Fig_2_ER, nrow = 2)
ggsave('data_v2/5_figs/lowFilter/agg/fig2_sppOccs.png', width = 4800, height = 3200, units = "px", bg = 'white')

distributionSpPlotWithoutLegendWY = function(sp){
    print(sp)
    output.filename = paste0("data_v2/5_figs/lowFilter/wy/fig2_", sp, "_2000_2023.png")
    input.filename = paste0("data_v2/4_maxent_outputs/wy/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_wy.tif")
    sppOccs.filename = paste0("data/1_occ/combined_spp_occ/", sp, "_lowFilter.csv")
    sppOccs = read_csv(sppOccs.filename)
    sppOccs_vect = sppOccs %>%
        vect(geom=c('lon', 'lat'), crs='epsg:4326')
    
    input.rast = rast(input.filename)
    rast_CA = input.rast %>%
        crop(CA, mask=T)
    rast_CV = input.rast %>%
        crop(central_valley, mask=T)
    
    vect_CV = rast_CV %>%
        as.polygons()
    sp_pretty = sp %>%
        str_replace("_", ". ") %>%
        str_to_sentence()
    
    # plot_CV = ggplot() +
    #     geom_spatraster(data=rast_CV) +
    #     geom_spatvector(data=sppOccs_vect, shape='plus') +
    #     # scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
    #     scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
    #     labs(
    #         fill = "Predicted suitability",
    #         title = sp_pretty
    #     ) +
    #     theme_void() +
    #     theme(
    #         # legend.position = 'bottom',
    #         legend.position = 'none',
    #         legend.direction = 'horizontal',
    #         legend.key.width = unit(1, 'cm'),
    #         legend.title.position = 'top',
    #         text = element_text("Times", face = 'italic')
    #     )
    #
    # CV = ifel(!is.na(rast_CV), 0, NA)
    
    plot_CA = ggplot() +
        geom_spatraster(data=rast_CA) +
        geom_spatvector(data=sppOccs_vect, shape=6, color='green', fill = 'transparent', size=1, alpha=0.2) +
        scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        new_scale_fill() +
        #Greyed out CV
        # geom_spatvector(data=vect_CV, na.rm=T, fill='transparent') +
        # scale_fill_gradientn(colours = c('grey20', 'grey21'), na.value='transparent') +
        theme_void() +
        labs(title = sp_pretty) +
        theme(legend.position = "none")
    # plot_inset = ggdraw(plot = plot_CV) +
    #     draw_plot(plot_CA, x = 0.48, y = 0.54, width = 0.4, height = 0.4)
    # plot_inset
    
    
    # ggsave(output.filename,
    #        scale=1.25, bg = 'white', units = 'px', width = 1080, height = 1920)
    
    return(plot_CA)
}

Fig_2_WY = map(names, distributionSpPlotWithoutLegendWY)
cowplot::plot_grid(plotlist = Fig_2_WY, nrow = 2)
ggsave('data_v2/5_figs/lowFilter/wy/fig2_sppOccs.png', width = 4800, height = 3200, units = "px", bg = 'white') 



# Test on CV BG data ------------------------------------------------------

externalTestingAUC = function(sp){
    model_agg.filename = paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    model_agg = readRDS(model_agg.filename)
    model_wy.filename = paste0('data_v2/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    model_wy = readRDS(model_wy.filename)
    
    #Internal testing
    internal_data_agg.filename = paste0('data_v2/3_swd/agg/testing_', sp, '_soil200cm_lowFilter_agg.csv')
    internal_data_agg = read_csv(internal_data_agg.filename)
    internal_data_wy.filename = paste0('data_v2/3_swd/wy/testing_', sp, '_soil200cm_lowFilter_wy.csv')
    internal_data_wy = read_csv(internal_data_wy.filename)
    
    internal_pred_agg = internal_data_agg %>% 
        as.data.frame() %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(model_agg, .))
    internal_pred_wy = internal_data_wy %>% 
        as.data.frame() %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(model_wy, .))
    
    internal_auc_agg = pROC::auc(internal_pred_agg$presence, internal_pred_agg$pred)
    internal_auc_wy = pROC::auc(internal_pred_wy$presence, internal_pred_wy$pred)
    
    #External testing
    external_data_agg.filename = paste0('data_cvBG/3_swd/agg/swd_', sp, '_soil200cm_lowFilter_agg.csv')
    external_data_agg = read_csv(external_data_agg.filename)
    external_data_wy.filename = paste0('data_cvBG/3_swd/wy/swd_', sp, '_soil200cm_lowFilter_wy.csv')
    external_data_wy = read_csv(external_data_wy.filename)
    
    external_pred_agg = external_data_agg %>% 
        as.data.frame() %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(model_agg, .))
    external_pred_wy = external_data_wy %>% 
        as.data.frame() %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(model_wy, .))
    
    external_auc_agg = pROC::auc(external_pred_agg$presence, external_pred_agg$pred)
    external_auc_wy = pROC::auc(external_pred_wy$presence, external_pred_wy$pred)
    
    out = tibble(
        sp = sp, 
        internal_auc_agg = internal_auc_agg[[1]],
        internal_auc_wy = internal_auc_wy[[1]],
        external_auc_agg = external_auc_agg[[1]],
        external_auc_wy = external_auc_wy[[1]]
    )
}

external_testing_auc_v2 = map(names, externalTestingAUC) %>% 
    bind_rows()

external_testing_auc_v2 %>% 
    pivot_longer(-sp) %>% 
    mutate(
        spec = ifelse(str_detect(name, '_agg'), 'Ecorelevant', 'WY'),
        testing_partition = ifelse(str_detect(name, 'internal'), 'Internal', 'External')
    ) %>% 
    ggplot() +
    geom_bar(aes(x=spec, y=value, fill=testing_partition, group=testing_partition), stat='identity', position='dodge') +
    facet_wrap(~sp, nrow=2)