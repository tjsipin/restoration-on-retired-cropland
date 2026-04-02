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
    er_mod_filename = paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    wy_mod_filename = paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    
    er_testing_filename = paste0('data_cvBG/3_swd/agg/testing_', sp, '_soil200cm_lowFilter_agg.csv')
    wy_testing_filename = paste0('data_cvBG/3_swd/wy/testing_', sp, '_soil200cm_lowFilter_wy.csv')
    
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
    er_mod_filename = paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    wy_mod_filename = paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
    
    er_training_filename = paste0('data_cvBG/3_swd/agg/training_', sp, '_soil200cm_lowFilter_agg.csv')
    wy_training_filename = paste0('data_cvBG/3_swd/wy/training_', sp, '_soil200cm_lowFilter_wy.csv')
    er_testing_filename = paste0('data_cvBG/3_swd/agg/testing_', sp, '_soil200cm_lowFilter_agg.csv')
    wy_testing_filename = paste0('data_cvBG/3_swd/wy/testing_', sp, '_soil200cm_lowFilter_wy.csv')
    
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

aucs_cvBG = map(names, consensusAUC) %>% 
    bind_rows() %>% 
    rename(`Ecorelevant AUC`=er_auc, `Water year AUC`=wy_auc, `Consensus AUC`=consensus_auc) %>% 
    pivot_longer(-Species)
ggplot(aucs_cvBG) +
    geom_bar(aes(x=Species, y=value, fill=name, group=name), stat='identity', position='dodge') +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    theme_light() +
    labs(y='AUC', fill='')

FORs_cvBG = map(names, consensusFOR) %>% 
    bind_rows() %>% 
    rename(`Ecorelevant FOR`=er_FOR, `Water year FOR`=wy_FOR, `Consensus FOR`=consensus_FOR) %>% 
    pivot_longer(-Species)
ggplot(FORs_cvBG) +
    geom_bar(aes(x=Species, y=value, fill=name, group=name), stat='identity', position='dodge') +
    scale_y_continuous(breaks = seq(0, 0.1, 0.005)) +
    theme_light() +
    labs(y='FOR', fill='')

consensusRast = function(sp, model_years){
    
    if(model_years=='2000_2023'){
        er_mod_path = paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_hist/')
        wy_mod_path = paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/monthly_dist_hist/')
        out_path = paste0('data_cvBG/4_maxent_outputs/consensus/', sp, '/lowFilter/monthly_dist_hist/')
    } else if(model_years=='MIROC45_2070_2099'){
        er_mod_path = paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_MIROC45/')
        wy_mod_path = paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/monthly_dist_MIROC45/')
        out_path = paste0('data_cvBG/4_maxent_outputs/consensus/', sp, '/lowFilter/monthly_dist_MIROC45/')
    } else if(model_years=='MIROC85_2070_2099'){
        er_mod_path = paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_MIROC85/')
        wy_mod_path = paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/monthly_dist_MIROC85/')
        out_path = paste0('data_cvBG/4_maxent_outputs/consensus/', sp, '/lowFilter/monthly_dist_MIROC85/')
    } 
    dir.create(out_path, recursive=T)
    out_filename = paste0(out_path, sp, '_', model_years, '_consensus.tif')
    
    er_rast_filename = paste0(
        er_mod_path, sp, '_', model_years, '_agg.tif'
    )
    wy_rast_filename = paste0(
        wy_mod_path, sp, '_', model_years, '_wy.tif'
    )
    
    er_rast = rast(er_rast_filename)
    wy_rast = rast(wy_rast_filename)
    
    consensus_rast = rast(list(er_rast, wy_rast)) %>% 
        mean()
    
    writeRaster(consensus_rast, out_filename)
}

map(names, ~consensusRast(.x, '2000_2023'))
map(names, ~consensusRast(.x, 'MIROC45_2070_2099'))
map(names, ~consensusRast(.x, 'MIROC85_2070_2099'))

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


varImpPlot = function(sp){
    erModel = readRDS(paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds'))
    wyModel = readRDS(paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds'))
    
    erVIP = erModel@results %>% 
        as.data.frame() %>%
        mutate(variable = rownames(.)) %>%
        tibble() %>%
        filter(str_detect(variable, 'permutation.importance')) %>%
        mutate(
            variable = str_replace(variable, '.permutation.importance', ''),
            sp = sp,
            model = 'Ecorelevant'
        ) %>%
        rename(importance = V1) %>%
        select(sp, model, variable, importance) %>% 
        arrange(desc(importance))
    wyVIP = wyModel@results %>% 
        as.data.frame() %>%
        mutate(variable = rownames(.)) %>%
        tibble() %>%
        filter(str_detect(variable, 'permutation.importance')) %>%
        mutate(
            variable = str_replace(variable, '.permutation.importance', ''),
            sp = sp,
            model = 'WY'
        ) %>%
        rename(importance = V1) %>%
        select(sp, model, variable, importance) %>% 
        arrange(desc(importance))
    combineVIP = erVIP %>% 
        rbind(wyVIP) %>% 
        mutate(sp_pretty = sp %>% str_replace('_', '. ') %>% str_to_sentence())
    return(combineVIP)
}

vips = map(names, varImpPlot) %>% 
    bind_rows()

ggplot(data=vips) + 
    geom_bar(aes(x=importance, y=variable, fill=model), stat = 'identity', position = 'dodge') +
    facet_wrap(~sp_pretty, nrow=2) +
    theme_light()

partialDependencePlot = function(sp){
    erModel = readRDS(paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds'))
    wyModel = readRDS(paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds'))
    
    erPDP_files = list.files(
        paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/model/plots'),
        pattern='.dat',
        full.names = T
    ) 
    erPDP_files = erPDP_files[!str_detect(erPDP_files, '_only')]
    
    wyPDP_files = list.files(
        paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/model/plots'),
        pattern='.dat',
        full.names = T
    ) 
    wyPDP_files = wyPDP_files[!str_detect(wyPDP_files, '_only')]
    
    erPDP = map(
        erPDP_files,
        function(fi){
            read_csv(fi) %>% 
                mutate(model = 'Ecorelevant')
        }
    ) %>% 
        bind_rows()
    wyPDP = map(
        wyPDP_files,
        function(fi){
            read_csv(fi) %>% 
                mutate(model = 'WY')
        }
    ) %>% 
        bind_rows()
    
    allPDP = erPDP %>% 
        bind_rows(wyPDP) %>% 
        group_by(variable) %>% 
        mutate(n_mod = paste0('n models = ', length(unique(model)))) %>% 
        ungroup()
    
    sp_pretty = sp %>% 
        str_replace('_', '. ') %>% 
        str_to_sentence()
    p = ggplot(allPDP) +
        geom_line(
            aes(x=x, y=y,color=model)
        ) +
        facet_wrap(~variable + n_mod, scales='free_x') +
        labs(
            title = paste0(sp_pretty, ' PDP')
        )
    return(p)
}
partialDependencePlot('l_californica')
partialDependencePlot('l_pentachaeta')
partialDependencePlot('a_menziesii')

# Plot --------------------------------------------------------------------

dir.create('data_cvBG/5_figs/lowFilter/consensus', recursive=T)

set.seed(123)
central_valley = vect('data/central_valley/ds2632.gdb') %>% 
    project('epsg:3310')
CA = tigris::states() %>% 
    vect() %>% 
    project('epsg:3310') %>% 
    filter(STUSPS=="CA")

# Figure 2: current distributions across species ---------------------------

distributionSpPlotWithoutLegend = function(sp){
    print(sp)
    output.filename = paste0("data_cvBG/5_figs/lowFilter/consensus/fig2_", sp, "_2000_2023.png")
    input.filename = paste0("data_cvBG/4_maxent_outputs/consensus/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_consensus.tif")
    input.rast = rast(input.filename)
    rast_CA = input.rast %>%
        crop(CA, mask=T)
    rast_CV = input.rast %>%
        crop(central_valley, mask=T)
    sp_pretty = sp %>%
        str_replace("_", ". ") %>%
        str_to_sentence()
    
    plot_CV = ggplot() +
        geom_spatraster(data=rast_CV) +
        # scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        labs(
            fill = "Predicted suitability",
            title = sp_pretty
        ) +
        theme_void() +
        theme(
            # legend.position = 'bottom',
            legend.position = 'none',
            legend.direction = 'horizontal',
            legend.key.width = unit(1, 'cm'),
            legend.title.position = 'top',
            text = element_text("Times", face = 'italic')
        )
    
    CV = ifel(!is.na(rast_CV), 0, NA)
    
    plot_CA = ggplot() +
        geom_spatraster(data=rast_CA) +
        scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        new_scale_fill() +
        #Greyed out CV
        geom_spatraster(data=CV, na.rm=T) +
        scale_fill_gradientn(colours = c('grey20', 'grey21'), na.value='transparent') +
        theme_void() +
        theme(legend.position = "none")
    plot_inset = ggdraw(plot = plot_CV) +
        draw_plot(plot_CA, x = 0.48, y = 0.54, width = 0.4, height = 0.4)
    plot_inset
    
    ggsave(output.filename,
           scale=1.25, bg = 'white', units = 'px', width = 1080, height = 1920)
    
    return(plot_inset)
}

Fig_2 = map(names, distributionSpPlotWithoutLegend)

#Extract legend
sp='l_gracilis'
input.legend.filename = paste0("data_cvBG/4_maxent_outputs/consensus/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_consensus.tif")
input.legend.rast = rast(input.legend.filename)
legend.rast = input.legend.rast %>%
    crop(central_valley, mask=T)
legend.rast_CA = input.legend.rast %>%
    crop(CA, mask=T)
sp_pretty = sp %>%
    str_to_sentence() %>%
    str_replace("_", ". ")
legend.plot = ggplot() +
    geom_spatraster(data=legend.rast) +
    scale_fill_gradientn("Predicted suitability", colours = c('navy', 'lightblue', 'lightyellow', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
    theme(text = element_text(family = 'Times'))  +
    labs(
        fill = "Predicted suitability",
        title = sp_pretty
    ) +
    theme_void() +
    theme(
        legend.position = 'right',
        # legend.position = 'none',
        legend.direction = 'vertical',
        legend.key.width = unit(1, 'cm'),
        legend.title.position = 'top',
        text = element_text("Times")
    )

legend = legend.plot %>%
    ggpubr::get_legend() %>%
    as_ggplot() +
    # ggspatial::annotation_scale(data = legend.rast, plot_unit = 'km') +
    ggspatial::annotation_north_arrow(
        data = legend.rast,
        location = "br", which_north = "true",
        pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
        style = ggspatial::north_arrow_nautical(
            fill = c("grey40", "white"),
            line_col = "grey20",
            text_family = "ArcherPro Book"
        )
    )

# legend_arrow = cowplot::plot_grid(legend, arrow, nrow = 1)

Fig_2_legend = Fig_2
Fig_2_legend[[11]] = legend
cowplot::plot_grid(plotlist = Fig_2_legend, nrow = 2)
ggsave('data_cvBG/5_figs/lowFilter/consensus/fig2.png', width = 2800, height = 2180, units = "px", bg = 'white')




# Difference rasters ------------------------------------------------------

central_valley = vect('data/central_valley/ds2632.gdb') %>% 
    project('epsg:3310')
ref_rast = rast('data/0_env/bcm/bcmv8_historic/2000_2023_monthly/aet1999dec.tif')
rast_CV = ref_rast %>%
    crop(central_valley, mask=T) 
rast_CV = ifel(!is.na(rast_CV), 1, NA)
names(rast_CV) = 'value'
vect_CV = rast_CV %>% 
    as.polygons()

diff.rasters = map(
    names,
    function(sp){
        input.raster.agg = rast(paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_hist/', sp, '_2000_2023_agg.tif'))
        names(input.raster.agg) = paste0(sp, '_agg')
        input.raster.wy = rast(paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/monthly_dist_hist/', sp, '_2000_2023_wy.tif'))
        names(input.raster.wy) = paste0(sp, '_wy')
        diff.raster = input.raster.agg - input.raster.wy
        sp_pretty = sp %>% 
            str_replace('_', '. ') %>% 
            str_to_sentence()
        names(diff.raster) = paste0(sp_pretty)
        return(diff.raster)
    }
) %>% 
    rast()

ggplot() +
    geom_spatraster(data=diff.rasters) +
    geom_spatvector(data=vect_CV, fill='transparent') +
    scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(-1, 1), na.value='transparent') +
    theme_void() + 
    labs(
        fill = "Difference in predicted suitability",
        title = "Ecorelevant (CV BG) - WY (CV BG)"
    ) +
    facet_wrap(~lyr) +
    theme(strip.text = element_text(face='italic'))



# Present points overlaid fig 2 -------------------------------------------

distributionSpPlotWithoutLegendER = function(sp){
    print(sp)
    output.filename = paste0("data_cvBG/5_figs/lowFilter/agg/fig2_", sp, "_2000_2023.png")
    input.filename = paste0("data_cvBG/4_maxent_outputs/agg/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_agg.tif")
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
ggsave('data_cvBG/5_figs/lowFilter/agg/fig2_sppOccs.png', width = 4800, height = 3200, units = "px", bg = 'white')

distributionSpPlotWithoutLegendWY = function(sp){
    print(sp)
    output.filename = paste0("data_cvBG/5_figs/lowFilter/wy/fig2_", sp, "_2000_2023.png")
    input.filename = paste0("data_cvBG/4_maxent_outputs/wy/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_wy.tif")
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
ggsave('data_cvBG/5_figs/lowFilter/wy/fig2_sppOccs.png', width = 4800, height = 3200, units = "px", bg = 'white')



distributionSpPlotWithoutLegendER_CV = function(sp){
    print(sp)
    output.filename = paste0("data_cvBG/5_figs/lowFilter/agg/fig2_", sp, "_2000_2023.png")
    input.filename = paste0("data_cvBG/4_maxent_outputs/agg/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_agg.tif")
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
    
    plot_CV = ggplot() +
        geom_spatraster(data=rast_CV) +
        geom_spatvector(data=sppOccs_vect, shape='plus') +
        # scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        labs(
            fill = "Predicted suitability",
            title = sp_pretty
        ) +
        theme_void() +
        theme(
            # legend.position = 'bottom',
            legend.position = 'none',
            legend.direction = 'horizontal',
            legend.key.width = unit(1, 'cm'),
            legend.title.position = 'top',
            text = element_text("Times", face = 'italic')
        )

    # CV = ifel(!is.na(rast_CV), 0, NA)
    # 
    # plot_CA = ggplot() +
    #     geom_spatraster(data=rast_CA) +
    #     geom_spatvector(data=sppOccs_vect, shape=6, color='green', fill = 'transparent', size=1, alpha=0.2) +
    #     scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
    #     new_scale_fill() +
    #     #Greyed out CV
    #     geom_spatvector(data=vect_CV, na.rm=T, fill='transparent') +
    #     scale_fill_gradientn(colours = c('grey20', 'grey21'), na.value='transparent') +
    #     theme_void() +
    #     labs(title = sp_pretty) +
    #     theme(legend.position = "none")
    # plot_inset = ggdraw(plot = plot_CV) +
    #     draw_plot(plot_CA, x = 0.48, y = 0.54, width = 0.4, height = 0.4)
    # plot_inset
    
    
    # ggsave(output.filename,
    #        scale=1.25, bg = 'white', units = 'px', width = 1080, height = 1920)
    
    return(plot_CV)
}

Fig_2_ER_CV = map(names, distributionSpPlotWithoutLegendER_CV)
cowplot::plot_grid(plotlist = Fig_2_ER_CV, nrow = 2)
ggsave('data_cvBG/5_figs/lowFilter/agg/fig2_sppOccs_CV.png', width = 4800, height = 3600, units = "px", bg = 'white')

distributionSpPlotWithoutLegendWY_CV = function(sp){
    print(sp)
    output.filename = paste0("data_cvBG/5_figs/lowFilter/wy/fig2_", sp, "_2000_2023.png")
    input.filename = paste0("data_cvBG/4_maxent_outputs/wy/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_wy.tif")
    sppOccs.filename = paste0("data/1_occ/combined_spp_occ/", sp, "_lowFilter.csv")
    sppOccs = read_csv(sppOccs.filename)
    sppOccs_vect = sppOccs %>%
        vect(geom=c('lon', 'lat'), crs='epsg:4326')
    
    input.rast = rast(input.filename)
    # rast_CA = input.rast %>%
    #     crop(CA, mask=T)
    rast_CV = input.rast %>%
        crop(central_valley, mask=T)
    
    vect_CV = rast_CV %>%
        as.polygons()
    sppOccs_vect = sppOccs_vect %>% 
        crop(vect_CV)
    sp_pretty = sp %>%
        str_replace("_", ". ") %>%
        str_to_sentence()
    
    plot_CV = ggplot() +
        geom_spatraster(data=rast_CV) +
        geom_spatvector(data=sppOccs_vect, shape='plus') +
        # scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
        labs(
            fill = "Predicted suitability",
            title = sp_pretty
        ) +
        theme_void() +
        theme(
            # legend.position = 'bottom',
            legend.position = 'none',
            legend.direction = 'horizontal',
            legend.key.width = unit(1, 'cm'),
            legend.title.position = 'top',
            text = element_text("Times", face = 'italic')
        )

    # CV = ifel(!is.na(rast_CV), 0, NA)
    # 
    # plot_CA = ggplot() +
    #     geom_spatraster(data=rast_CA) +
    #     geom_spatvector(data=sppOccs_vect, shape=6, color='green', fill = 'transparent', size=1, alpha=0.2) +
    #     scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(0, 1), na.value='transparent') +
    #     new_scale_fill() +
    #     #Greyed out CV
    #     geom_spatvector(data=vect_CV, na.rm=T, fill='transparent') +
    #     scale_fill_gradientn(colours = c('grey20', 'grey21'), na.value='transparent') +
    #     theme_void() +
    #     labs(title = sp_pretty) +
    #     theme(legend.position = "none")
    # plot_inset = ggdraw(plot = plot_CV) +
    #     draw_plot(plot_CA, x = 0.48, y = 0.54, width = 0.4, height = 0.4)
    # plot_inset
    
    
    # ggsave(output.filename,
    #        scale=1.25, bg = 'white', units = 'px', width = 1080, height = 1920)
    
    return(plot_CV)
}

Fig_2_WY_CV = map(names, distributionSpPlotWithoutLegendWY_CV)
cowplot::plot_grid(plotlist = Fig_2_WY_CV, nrow = 2)
ggsave('data_cvBG/5_figs/lowFilter/wy/fig2_sppOccs_CV.png', width = 4800, height = 3200, units = "px", bg = 'white')



# Difference rasters ------------------------------------------------------

diffRast = function(sp, model_years){
    sp_pretty = sp %>% str_replace('_', '. ') %>% str_to_sentence()
    if(model_years=='2000_2023'){
        er_mod_path = paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_hist/')
        wy_mod_path = paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/monthly_dist_hist/')
        out_path = paste0('data_cvBG/4_maxent_outputs/consensus/', sp, '/lowFilter/monthly_dist_hist/')
    } else if(model_years=='MIROC45_2070_2099'){
        er_mod_path = paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_MIROC45/')
        wy_mod_path = paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/monthly_dist_MIROC45/')
        out_path = paste0('data_cvBG/4_maxent_outputs/consensus/', sp, '/lowFilter/monthly_dist_MIROC45/')
    } else if(model_years=='MIROC85_2070_2099'){
        er_mod_path = paste0('data_cvBG/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_MIROC85/')
        wy_mod_path = paste0('data_cvBG/4_maxent_outputs/wy/', sp, '/lowFilter/monthly_dist_MIROC85/')
        out_path = paste0('data_cvBG/4_maxent_outputs/consensus/', sp, '/lowFilter/monthly_dist_MIROC85/')
    } 
    dir.create(out_path, recursive=T)
    out_filename = paste0(out_path, sp, '_', model_years, '_consensus.tif')
    
    er_rast_filename = paste0(
        er_mod_path, sp, '_', model_years, '_agg.tif'
    )
    wy_rast_filename = paste0(
        wy_mod_path, sp, '_', model_years, '_wy.tif'
    )
    
    er_rast = rast(er_rast_filename)
    wy_rast = rast(wy_rast_filename)
    
    diff_rast = er_rast - wy_rast
    names(diff_rast) = sp_pretty
    return(diff_rast)
}

diff_rasts = rast(map(names, ~diffRast(.x, '2000_2023')))
ref_rast = rast('data/0_env/bcm/bcmv8_historic/2000_2023_monthly/aet1999dec.tif')
central_valley = vect('data/central_valley/ds2632.gdb') %>% 
    project('epsg:3310')

rast_CV = ref_rast %>%
    crop(central_valley, mask=T) 
rast_CV = ifel(!is.na(rast_CV), 1, NA)
names(rast_CV) = 'value'
vect_CV = rast_CV %>% 
    as.polygons()

ggplot() + 
    geom_spatraster(data=diff_rasts) +
    geom_spatvector(data=vect_CV, fill='transparent') +
    scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(-1, 1), na.value='transparent') +
    theme_void() + 
    labs(
        fill = "Difference",
        title = "Ecorelevant (CV BG) - WY (CV BG)"
    ) +
    facet_wrap(~lyr, nrow=2)
