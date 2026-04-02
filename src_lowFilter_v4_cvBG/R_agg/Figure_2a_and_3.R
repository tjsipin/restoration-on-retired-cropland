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
library(shades)
library(ggnewscale) #citation("ggnewscale")
library(cowplot)
tidymodels_prefer()

central_valley = vect('data/ds2632/ds2632.gdb') %>% 
    project('epsg:3310')
CA = tigris::states() %>% 
    vect() %>% 
    project('epsg:3310') %>% 
    filter(STUSPS=="CA")

## sp to run
names <- c(
    "a_polycarpa",
    "p_arborea",
    "c_pungens",
    "l_pentachaeta",
    "p_ciliata",
    "a_menziesii",
    "c_lasiophyllus"
)
## list of months
months = c("jan","feb", "mar","apr", "may","jun",
           "jul","aug","sep","oct","nov", "dec")

sp_highlow_grid = expand_grid(
    sp = names
)

## P10 threshold raster per species and model-years combination
getThresholdRasts = function(sp, model_years=c("2000_2023", "MIROC45_2070_2099", "MIROC85_2070_2099")){
    #Read in training data and model
    training.filename = paste0('data/agg/swd/training_', sp, '_soil200cm_lowFilter072925.csv')
    model.filename = paste0("data/agg/maxent_outputs_072925/", sp, "/lowFilter/model/", sp, "_sdm072925.rds")
    
    training = read_csv(training.filename)
    model = readRDS(model.filename)
    
    #Subset training data to occurrences only
    training_1pred = training %>% 
        filter(presence==1) %>% 
        mutate(tdiff = tmx - tmn) %>% 
        filter(complete.cases(.)) %>% 
        mutate(pred = predict(model, .)) 
    
    #Get the 10th percentile of occurrence point predictions
    p10 = quantile(training_1pred$pred, 0.1)
    
    if(str_detect(model_years, "MIROC45")){
        monthly_dist_str = "monthly_dist_MIROC45/"
    } else if(str_detect(model_years, "MIROC85")){
        monthly_dist_str = "monthly_dist_MIROC85/"
    } else if(str_detect(model_years, "2023")){
        monthly_dist_str = "monthly_dist_hist/"
    }
    sp_theme = sp %>% 
        str_replace("_", ". ") %>% 
        str_to_sentence()
    
    input.filenames = list.files(
        paste0('data/agg/maxent_outputs_072925/', sp, '/lowFilter/', monthly_dist_str),
        full.names = T
    )
    
    output.rast = rast(input.filenames)
    output.rast = ifel(output.rast > p10, 1, 0)
    output.rast = output.rast %>% sum()
    # output.rast = output.rast %>% 
    #     crop(central_valley, mask=T)
    
    dir.create(paste0('data/agg/maxent_outputs_072925/', sp, '/lowFilter/p10/'), recursive=T, showWarnings=F)
    output.filename = paste0('data/agg/maxent_outputs_072925/', sp, '/lowFilter/p10/', '/p10_', sp, '_', model_years, '_ecorelevant.tif')
    writeRaster(output.rast, output.filename, overwrite = T)
}

map(
    .x = names,
    .f = ~ getThresholdRasts(sp = .x, '2000_2023')
) 
map(
    .x = names,
    .f = ~ getThresholdRasts(sp = .x, 'MIROC45_2070_2099')
)
map(
    .x = names,
    .f = ~ getThresholdRasts(sp = .x, 'MIROC85_2070_2099')
)


# Create color palette
pl <- rev(scico(8, palette = 'managua'))

# Set color for each class
pal<-c(
    "0"=pl[1],
    "1"=pl[2],
    "2"=pl[3],
    "3"=pl[4],
    "4"=pl[5],
    "5"=pl[6],
    "6"=pl[7],
    "7"=pl[8]
)

# Color for na values
na_col<-'grey50'
# Create custom theme
custom_theme<-theme_void()+
    theme(plot.background = element_rect(fill="#202020",color=NA))

# Figure 3 ----------------------------------------------------------------

#Code to plot sum of suitable species per model-years
p10VisualizeAgg = function(model){
    print(model)
    out.dir = paste0("data/agg/maxent_outputs_072925/p10/")
    dir.create(out.dir, recursive=T)
    out.rast.filename.CA = paste0("data/agg/maxent_outputs_072925/p10/", model, "_lowFilter_sum_CA_072925.tif")
    out.rast.filename.CV = paste0("data/agg/maxent_outputs_072925/p10/", model, "_lowFilter_sum_CV_072925.tif")
    out.png.filename = paste0("data/agg/maxent_outputs_072925/figs/", model, "_lowFilter_sum_072925.png")
    input_rast.filenames = list.files(
        paste0("data/agg/maxent_outputs_072925/"),
        pattern=paste0("p10_.*", model, '_ecorelevant'),
        recursive = T,
        full.names=T
    )
    input_rast_CA = rast(input_rast.filenames) %>% 
        sum() %>% 
        mutate(sum = factor(sum, levels = as.integer(0:7))) %>% 
        crop(CA, mask=T)
    
    input_rast_CV = input_rast_CA %>% 
        crop(central_valley, mask=T) 
    
    writeRaster(input_rast_CA, out.rast.filename.CA, overwrite = T)
    writeRaster(input_rast_CV, out.rast.filename.CV, overwrite = T)
    
    if(str_detect(model, "2000_2023")){
        model_pretty = "2000 to 2023"
        subtitle = "Current conditions\nEcorelevant period aggregate models"
    } else if(str_detect(model, "MIROC45")){
        model_pretty = "end of century"
        subtitle = "MIROC 4.5\nEcorelevant period aggregate models"
    } else if(str_detect(model, "MIROC85")){
        model_pretty = "end of century"
        subtitle = "MIROC 8.5\nEcorelevant period aggregate models"
    }
    
    plot_CV = ggplot() +
        geom_spatraster(data = input_rast_CV, show.legend=T) +
        scale_fill_manual("Number of species suitable", values=pal, drop=F, na.translate=F)+ #drop=F, na.translate=F to show all 8 colors in legend
        labs(title = paste0("Sum of suitable species for ", model_pretty), 
             subtitle = subtitle) +
        theme_void() +
        theme(
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            legend.key.width = unit(1, 'cm'),
            text = element_text("Times")
        ) +
        guides(fill=guide_legend(
            nrow=1,byrow=TRUE, 
            theme = theme(legend.text.position = 'bottom', legend.title.position = 'top')
        ))
    CV = ifel(!is.na(rast_CV), 0, NA)
    
    plot_CA = ggplot() +
        geom_spatraster(data=input_rast_CA, show.legend = F) +
        scale_fill_manual(values=pal, drop=F, na.translate=F)+ #drop=F, na.translate=F to show all 8 colors in legend
        guides(fill="none") +
        new_scale_fill() +
        #Greyed out CV
        geom_spatraster(data=CV, na.rm=T) +
        scale_fill_gradientn(colours = c('grey20', 'grey21'), na.value='transparent') +
        theme_void() +
        theme(legend.position = "none")
    plot_inset = ggdraw(plot = plot_CV) +
        draw_plot(plot_CA, x = 0.48, y = 0.48, width = 0.4, height = 0.4)
    plot_inset
    
    ggsave(
        out.png.filename, 
        scale=1.25, bg = 'white', units = 'px', width = 1440, height = 1920
    )
}

p10VisualizeAgg(model = "2000_2023")
p10VisualizeAgg(model = "MIROC45_2070_2099")
p10VisualizeAgg(model = "MIROC85_2070_2099")


# Figure 2A: current distributions across species ---------------------------

distributionSpPlot = function(sp){
    print(sp)
    output.filename = paste0("data/agg/maxent_outputs_072925/figs/", sp, "_2000_2023_ecorelevant_fig2A.png")
    input.filename = paste0("data/agg/maxent_outputs_072925/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_ecorelevant.tif")
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
        scale_fill_gradientn(colours = c('orange4', 'white', 'navy'), limits=c(0, 1), na.value='transparent') +
        labs(
            fill = "Predicted suitability",
            title = "Current species distribution",
            subtitle = sp_pretty
        ) +
        theme_void() +
        theme(
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            legend.key.width = unit(1, 'cm'),
            legend.title.position = 'top',
            text = element_text("Times")
        ) 
    
    CV = ifel(!is.na(rast_CV), 0, NA)
        
    plot_CA = ggplot() +
        geom_spatraster(data=rast_CA) +
        scale_fill_gradientn(colours = c('orange4', 'white', 'navy'), limits=c(0, 1), na.value='transparent') +
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
    
}

map(names, distributionSpPlot)
