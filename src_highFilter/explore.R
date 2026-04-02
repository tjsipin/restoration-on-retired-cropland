library(terra)
library(tidyterra)
library(dplyr)
library(stringr)
library(readr)
library(vip)

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
    erModel = readRDS(paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds'))
    wyModel = readRDS(paste0('data_v2/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds'))
    
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
        rbind(wyVIP)
    return(combineVIP)
}

vips = map(names, varImpPlot) %>% 
    bind_rows()

ggplot(data=vips) + 
    geom_bar(aes(x=importance, y=variable, fill=model), stat = 'identity', position = 'dodge') +
    facet_wrap(~sp)

partialDependencePlot = function(sp){
    erModel = readRDS(paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds'))
    wyModel = readRDS(paste0('data_v2/4_maxent_outputs/wy/', sp, '/lowFilter/model/', sp, '_final_sdm.rds'))
    
    erPDP_files = list.files(
        paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/model/plots'),
        pattern='.dat',
        full.names = T
    ) 
    erPDP_files = erPDP_files[!str_detect(erPDP_files, '_only')]
    
    wyPDP_files = list.files(
        paste0('data_v2/4_maxent_outputs/wy/', sp, '/lowFilter/model/plots'),
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

a_menziesii_pdp = partialDependencePlot('a_menziesii')
a_intermedia_pdp = partialDependencePlot('a_intermedia')
c_pungens_pdp = partialDependencePlot('c_pungens')
l_pentachaeta_pdp = partialDependencePlot('l_pentachaeta')
p_arborea_pdp = partialDependencePlot('p_arborea')
p_ciliata_pdp = partialDependencePlot('p_ciliata')
a_polycarpa_pdp = partialDependencePlot('a_polycarpa')
c_lasiophyllus_pdp = partialDependencePlot('c_lasiophyllus')
l_californica_pdp = partialDependencePlot('l_californica')
l_gracilis_pdp = partialDependencePlot('l_gracilis')