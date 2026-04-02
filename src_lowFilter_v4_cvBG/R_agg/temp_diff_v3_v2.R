library(terra)
library(tidyterra)
library(dplyr)
library(stringr)
library(readr)

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

ref_rast = rast('data_bufferBG/0_env/bcm/bcmv8_historic/2000_2023_monthly/aet1999dec.tif')
central_valley = vect('data_bufferBG/central_valley/ds2632.gdb') %>% 
    project('epsg:3310')

rast_CV = ref_rast %>%
    crop(central_valley, mask=T) 
rast_CV = ifel(!is.na(rast_CV), 1, NA)
names(rast_CV) = 'value'
vect_CV = rast_CV %>% 
    as.polygons()

diff.rasters = map(
    names,
    function(sp){
        input.raster.v2 = rast(paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_hist/', sp, '_2000_2023_agg.tif'))
        names(input.raster.v2) = paste0(sp, '_v2')
        input.raster.v3 = rast(paste0('data_bufferBG/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_hist/', sp, '_2000_2023_agg.tif'))
        names(input.raster.v3) = paste0(sp, '_v3')
        diff.raster = input.raster.v3 - input.raster.v2
        names(diff.raster) = paste0(sp, '_diff')
        return(diff.raster)
    }
) %>% 
    rast()

occsPerSp = map(
    names,
    function(sp){
        input.csv = paste0('data_bufferBG/1_occ/combined_spp_occ/', sp, '_lowFilter.csv') %>% 
            read_csv() %>% 
            mutate(sp = sp)
        return(input.csv)
    }
) %>% 
    bind_rows() %>% 
    vect(geom=c('lon', 'lat'), crs='epsg:4326') %>% 
    project('epsg:3310')

ggplot() +
    geom_spatraster(data=diff.rasters) +
    geom_spatvector(data=vect_CV, fill='transparent') +
    geom_spatvector(data=occsPerSp, geom='cross') +
    scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(-1, 1), na.value='transparent') +
    theme_void() + 
    labs(
        fill = "Difference in predicted suitability",
        title = "Buffer BG - Convex Hull",
        subtitle = 'Ecorelevant agg'
    ) +
    facet_wrap(~lyr)

diff.rasters.v2 = map(
    names,
    function(sp){
        input.raster.agg = rast(paste0('data_v2/4_maxent_outputs/agg/', sp, '/lowFilter/monthly_dist_hist/', sp, '_2000_2023_agg.tif'))
        names(input.raster.agg) = paste0(sp, '_agg')
        input.raster.wy = rast(paste0('data_v2/4_maxent_outputs/wy/', sp, '/lowFilter/monthly_dist_hist/', sp, '_2000_2023_wy.tif'))
        names(input.raster.wy) = paste0(sp, '_wy')
        diff.raster = input.raster.agg - input.raster.wy
        names(diff.raster) = paste0(sp, '_diff')
        return(diff.raster)
    }
) %>% 
    rast()
ggplot() +
    geom_spatraster(data=diff.rasters.v2) +
    geom_spatvector(data=vect_CV, fill='transparent') +
    scale_fill_gradientn(colours = c('navy', 'lightblue', 'lightyellow2', 'goldenrod1', 'red'), limits=c(-1, 1), na.value='transparent') +
    theme_void() + 
    labs(
        fill = "Difference in predicted suitability",
        title = "Ecorelevant (v2) - WY (v2)"
    ) +
    facet_wrap(~lyr)
