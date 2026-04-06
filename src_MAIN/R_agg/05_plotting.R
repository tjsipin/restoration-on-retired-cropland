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

dir.create('data_MAIN/5_figs/lowFilter/agg', recursive=T)

set.seed(123)
central_valley = vect('data/central_valley/ds2632.gdb') %>% 
    project('epsg:3310')
CA = tigris::states() %>% 
    vect() %>% 
    project('epsg:3310') %>% 
    filter(STUSPS=="CA")

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
# Figure 2: current distributions across species ---------------------------

distributionSpPlotWithoutLegend = function(sp){
    print(sp)
    output.filename = paste0("data_MAIN/5_figs/lowFilter/agg/fig2_", sp, "_2000_2023.png")
    input.filename = paste0("data_MAIN/4_maxent_outputs/agg/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_agg.tif")
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
input.legend.filename = paste0("data_MAIN/4_maxent_outputs/agg/", sp, "/lowFilter/monthly_dist_hist/", sp, "_2000_2023_agg.tif")
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
ggsave('data_MAIN/5_figs/lowFilter/agg/fig2.png', width = 2800, height = 2180, units = "px", bg = 'white')




# Figure 3 ---------------------------------------------------------------

#'@author Madi Calbert
# Species and file root (matches your screenshot/naming)
species <- names
agg_root <- "data_MAIN/4_maxent_outputs/agg"

# Read & dissolve the Great Valley ecoregion (once)
gv_path <- "data/central_valley/ds2632.gdb/"
sf_use_s2(TRUE)
gv_raw <- st_read(gv_path, quiet = TRUE)
if (is.na(st_crs(gv_raw))) st_crs(gv_raw) <- 4326
gv_raw <- gv_raw |>
    st_zm(drop = TRUE, what = "ZM") |>
    st_cast("MULTIPOLYGON", warn = FALSE) |>
    st_make_valid()
sf_use_s2(FALSE)
gv_union <- gv_raw |>
    mutate(one = 1L) |>
    group_by(one) |>
    summarise(.groups = "drop") |>
    st_make_valid()

# Change coder: 0=no detect, 1=Loss, 2=Gain, 3=Stable
calc_change_pair <- function(h, f) {
    terra::ifel(is.na(h) | is.na(f), NA,
                terra::ifel(h == 1 & f == 0, 1,
                            terra::ifel(h == 0 & f == 1, 2,
                                        terra::ifel(h == 1 & f == 1, 3, 0))))
}

# Build one eco-relevant change layer for a species+scenario (keeps zeros)
change_layer_for <- function(spp, scenario = c("MIROC45","MIROC85")) {
    scenario <- match.arg(scenario)
    p10_dir  <- file.path(agg_root, spp, "lowFilter", "p10")
    f_hist   <- file.path(p10_dir, sprintf("p10_%s_2000_2023_agg.tif", spp))
    f_fut    <- file.path(p10_dir, sprintf("p10_%s_%s_2070_2099_agg.tif", spp, scenario))
    stopifnot(file.exists(f_hist), file.exists(f_fut))

    r_hist <- terra::rast(f_hist)
    r_fut  <- terra::rast(f_fut)
    if (!terra::compareGeom(r_hist, r_fut, stopOnError = FALSE)) {
        r_fut <- terra::resample(r_fut, r_hist, method = "near")
    }

    # Transform GV to the raster CRS
    crs_wkt   <- terra::crs(r_hist)
    target_crs <- tryCatch(sf::st_crs(crs_wkt), error = function(e) sf::st_crs(3310))
    if (is.na(target_crs)) target_crs <- sf::st_crs(3310)
    gv_t  <- sf::st_transform(gv_union, target_crs)
    gv_tv <- terra::vect(gv_t)

    chg <- calc_change_pair(r_hist, r_fut)                   # 0..3 (keeps 0)
    chg <- terra::mask(terra::crop(chg, gv_tv), gv_tv)       # clip to GV
    chg <- as.factor(chg)                                    # categorical
    names(chg) <- spp
    chg
}

# Build stacks across species for each scenario
chg_45 <- terra::rast(lapply(species, change_layer_for, scenario = "MIROC45"))
chg_85 <- terra::rast(lapply(species, change_layer_for, scenario = "MIROC85"))

# Transform GV for plotting in each CRS
gv45 <- st_transform(gv_union, st_crs(terra::crs(chg_45)))
gv85 <- st_transform(gv_union, st_crs(terra::crs(chg_85)))

# Palette + labels
pal_fig3 <- c("0"="#9C9C9C","1"="#A2C7FF","2"="#2E8B57","3"="#F6D26D")
labs <- c("0"="Unsuitable Habitat",
          "1"="Habitat Loss",
          "2"="Habitat Gain",
          "3"="No Change")
facet_labs <- setNames(gsub("_", " ", species), species)

# Choose the font family
tnr_family <- "Times"

facet_labs <- c(
    "a_menziesii"    = "A. menziesii",
    'a_intermedia'   = "A. intermedia",
    "c_pungens"      = "C. pungens",
    "l_pentachaeta"  = "L. pentachaeta",
    "p_arborea"      = "P. arborea",
    "p_ciliata"      = "P. ciliata",
    "a_polycarpa"    = "A. polycarpa",
    "c_lasiophyllus" = "C. lasiophyllus",
    "l_californica"  = "L. californica",
    "l_gracilis"     = "L. gracilis"
)

ggRastLayer = function(i, chg, gv, show_compass = F, model=c('RCP4.5', 'RCP8.5')){
    r = chg[[i]]
    if(model=='RCP 4.5'){
        pretty_title = "Current to MIROC RCP 4.5 Distribution"
    } else if(model=="RCP 8.5"){
        pretty_title = "Current to MIROC RCP 8.5 Distribution"
    }
    p = ggplot() +
        geom_sf(data = gv, fill = "grey90", color = "grey80", linewidth = 0.2) +
        geom_spatraster(data = r, na.rm = TRUE) +
        facet_wrap(~lyr, ncol = 4, labeller = labeller(lyr = facet_labs)) +
        scale_fill_manual(name = NULL, values = pal_fig3, breaks = names(pal_fig3),
                          labels = NULL, na.translate = FALSE) +
        coord_sf(expand = FALSE) +
        theme_void(base_size = 14, base_family = tnr_family) +
        theme(
            text = element_text(family = tnr_family, size = 16),
            strip.text = element_text(family = tnr_family, face = "italic",
                                      hjust = 0, vjust = 1.15, size = 14),
            strip.clip = 'off',
            legend.position = "bottom",
            legend.text  = element_text(family = tnr_family, size = 16),
            legend.title = element_text(family = tnr_family, size = 16),
            legend.key.height = unit(6, "pt"),
            legend.key.width  = unit(10, "pt"),
            legend.spacing.y  = unit(2, "pt")
        )
    if(show_compass){
        p = p +
            ggspatial::annotation_scale(location = "bl", width_hint = 0.25,
                                        text_cex = 0.8, line_width = 0.3,
                                        pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm")) +
            ggspatial::annotation_north_arrow(
                location = "bl", which_north = "true",
                style = ggspatial::north_arrow_nautical(
                    fill = c("grey35","grey85"), line_col = "grey10"
                ),
                height = unit(1.1, "cm"), width = unit(1.1, "cm"),
                pad_x = unit(0.25, "cm"), pad_y = unit(1.10, "cm")
            )
    }
    return(p)
}
draw_label_theme <- function(label, theme = NULL, element = "text", title_size, ...) {
    theme <- theme(
        text = element_text(family = tnr_family),
        plot.title = element_text(family = tnr_family, margin = margin(b = 5), size = title_size)
    )
    if (!element %in% names(theme)) {
        stop("Element must be a valid ggplot theme element name")
    }

    elements <- ggplot2::calc_element(element, theme)

    cowplot::draw_label(label,
                        fontfamily = elements$family,
                        fontface = elements$face,
                        colour = elements$color,
                        size = elements$size,
                        ...
    )
}


# --- Make a vertical legend grob to place on the RIGHT -----------------------
vertical_legend_fig3 = (
    ggplot() +
        geom_sf(data = gv45, fill = "grey90", color = "grey80", linewidth = 0.2) +
        geom_spatraster(data = chg_45, na.rm = TRUE) +
        facet_wrap(~lyr, ncol = 4, labeller = labeller(lyr = facet_labs)) +
        scale_fill_manual(name = NULL, values = pal_fig3, breaks = names(pal_fig3),
                          labels = labs, na.translate = FALSE) +
        coord_sf(expand = FALSE) +
        theme_void(base_size = 17, base_family = tnr_family) +
        theme(
            legend.position      = "right",
            legend.box           = "vertical",
            legend.justification = "center",
            legend.key.width     = unit(1.0, "cm"),
            legend.key.height    = unit(0.6, "cm"),
            legend.title         = element_text(hjust = 0.5, size = 20),
            legend.text          = element_text(size = 14),
            legend.box.margin    = margin(t = 0, r = 16, b = 0, l = 8),
            plot.margin          = margin(t = 0, r = 16, b = 0, l = 0),
            text                 = element_text(family = "Times")
        ) +
        guides(fill = guide_legend(
            title.position = "top", title.hjust = 0.5,
            ncol = 1, byrow = FALSE, label.position = "right"
        ))
) %>%
    cowplot::get_legend()


p_45_1 = ggRastLayer(1, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_2 = ggRastLayer(2, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_3 = ggRastLayer(3, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_4 = ggRastLayer(4, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_5 = ggRastLayer(5, chg=chg_45, gv=gv45, model='RCP4.5', show_compass = T)
p_45_6 = ggRastLayer(6, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_7 = ggRastLayer(7, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_8 = ggRastLayer(8, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_9 = ggRastLayer(9, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_10 = ggRastLayer(10, chg=chg_45, gv=gv45, model='RCP4.5')
p_45_cowplot = cowplot::plot_grid(
    p_45_1, p_45_2, p_45_3, p_45_4,
    p_45_5, p_45_6, p_45_7, p_45_8, p_45_9, p_45_10, vertical_legend_fig3,
    ncol = 4, rel_widths = 1, align = "h"
)
p_45_title = ggdraw() +
    draw_label_theme("Current to MIROC RCP 4.5 Distribution (Ecorelevant model)", element = 'plot.title', title_size = 20)
p_45_final = cowplot::plot_grid(
    p_45_title,
    p_45_cowplot,
    ncol = 1, rel_heights = c(0.05, 1), align = 'h'
)
ggsave(
    'data_MAIN/5_figs/lowFilter/agg/fig3_RCP45.png',
    plot = p_45_final,
    width = 3400, height = 3000, units = "px", scale = 1,
    bg = "white", limitsize = FALSE
)

p_85_1 = ggRastLayer(1, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_2 = ggRastLayer(2, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_3 = ggRastLayer(3, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_4 = ggRastLayer(4, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_5 = ggRastLayer(5, chg=chg_85, gv=gv85, model='RCP8.5', show_compass = T)
p_85_6 = ggRastLayer(6, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_7 = ggRastLayer(7, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_8 = ggRastLayer(8, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_9 = ggRastLayer(9, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_10 = ggRastLayer(10, chg=chg_85, gv=gv85, model='RCP8.5')
p_85_cowplot = cowplot::plot_grid(
    p_85_1, p_85_2, p_85_3, p_85_4,
    p_85_5, p_85_6, p_85_7, p_85_8, p_85_9, p_85_10, vertical_legend_fig3,
    ncol = 4, rel_widths = 1, align = "h"
)
p_85_title = ggdraw() +
    draw_label_theme("Current to MIROC RCP 8.5 Distribution (Ecorelevant model)", element = 'plot.title', title_size = 20)
p_85_final = cowplot::plot_grid(
    p_85_title,
    p_85_cowplot,
    ncol = 1, rel_heights = c(0.05, 1), align = 'h'
)
ggsave(
    'data_MAIN/5_figs/lowFilter/agg/fig3_RCP85.png',
    plot = p_85_final,
    width = 3400, height = 3000, units = "px",
    bg = "white", limitsize = FALSE
)

# Figure 4 ----------------------------------------------------------------
#' @author Madi, TJ

# --- Palette (8 discrete colors for classes 0..7) ----------------------------
# Palette: blue → yellow → orange → red (8 steps for 0..7)
suit_grad <- colorRampPalette(c(
    "#1e3e95",  # deep blue
    "#8fb0d9",  # light blue
    "#F2D544",  # yellow (midpoint)
    "#F5CA96",  # light orange
    "#F1853F",  # orange
    "#CC331A"   # red
))
pal_fig4 <- setNames(suit_grad(11), as.character(0:10))


# --- Helper: nice labels -----------------------------------------------------
###
# Panel letters
.panel_letter <- function(model) {
    if (grepl("2000_2023", model)) "A"
    else if (grepl("MIROC45", model)) "B"
    else "C"
}

# --- Build ONE panel (legend-free), compass optional -------------------------
p10Panel <- function(model, show_compass = FALSE) {
    message("Rendering: ", model)
    base_dir = 'data_MAIN/4_maxent_outputs/agg/p10/'
    input_rast_CA_name <- paste0(base_dir, model, "_lowFilter_sum_CA.tif")
    input_rast_CV_name <- paste0(base_dir, model, "_lowFilter_sum_CV.tif")

    input_rast_CA = rast(input_rast_CA_name)
    input_rast_CV = rast(input_rast_CV_name)

    names(input_rast_CA) <- "sum"
    names(input_rast_CV) <- "sum"

    # categorical classes 0..7 so legend shows all bins
    input_rast_CA <- input_rast_CA %>%
        mutate(sum = factor(sum, levels = 0:10))
    input_rast_CV <- input_rast_CV %>%
        mutate(sum = factor(sum, levels = 0:10))

    # grey mask over CA (so CV pops)
    CV_mask <- terra::ifel(!is.na(input_rast_CV[[1]]), 0, NA)

    p_cv <- ggplot() +
        geom_spatraster(data = input_rast_CV, aes(fill = sum), show.legend = TRUE) +
        scale_fill_manual(
            name = "# of species",
            values = pal_fig4, drop = FALSE, na.translate = FALSE
        ) +
        theme_void()

    if (isTRUE(show_compass)) {
        p_cv <- p_cv +
            ggspatial::annotation_scale(location = "bl", width_hint = 0.25,
                                        text_cex = 0.8, line_width = 0.3,
                                        pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm")) +
            ggspatial::annotation_north_arrow(
                location = "bl", which_north = "true",
                style = ggspatial::north_arrow_nautical(
                    fill = c("grey35","grey85"), line_col = "grey10"
                ),
                height = unit(1.1, "cm"), width = unit(1.1, "cm"),
                pad_x = unit(0.25, "cm"), pad_y = unit(1.10, "cm")
            )
    }

    p_ca <- ggplot() +
        geom_spatraster(data = input_rast_CA, aes(fill = sum), show.legend = FALSE) +
        scale_fill_manual(values = pal_fig4, drop = FALSE, na.translate = FALSE) +
        guides(fill = "none") +
        ggnewscale::new_scale_fill() +
        geom_spatraster(data = CV_mask, na.rm = TRUE, show.legend = FALSE) +
        scale_fill_gradientn(colours = c("grey20", "grey21"),
                             na.value = "transparent", guide = "none") +
        theme_void()

    # compose (legend removed so all panels are same size)
    cowplot::ggdraw() +
        cowplot::draw_plot(p_cv + theme(legend.position = "none"),
                           x = 0, y = 0, width = 1, height = 1) +
        cowplot::draw_label(.panel_letter(model),
                            x = 0.985, y = 0.985, hjust = 1, vjust = 1,
                            fontface = "bold", size = 18, fontfamily = "Times") +
        cowplot::draw_plot(p_ca, x = 0.48, y = 0.48, width = 0.4, height = 0.4)
}

# --- Make a vertical legend grob to place on the RIGHT -----------------------
make_vertical_legend <- function(model_for_scale = "MIROC85_2070_2099") {
    base_dir = 'data_MAIN/4_maxent_outputs/agg/p10/'
    input_rast_CV_name <- paste0(base_dir, model_for_scale, "_lowFilter_sum_CV.tif")
    input_rast_CV = rast(input_rast_CV_name)
    names(input_rast_CV) <- "sum"
    input_rast_CV <- input_rast_CV %>%
        mutate(sum = factor(sum, levels = 0:10))


    base <- ggplot() +
        geom_spatraster(data = input_rast_CV, aes(fill = sum), show.legend = TRUE) +
        scale_fill_manual(
            name = "# of species",
            values = pal_fig4, drop = FALSE, na.translate = FALSE
        ) +
        theme_void() +
        theme(
            legend.position      = "right",
            legend.box           = "vertical",
            legend.justification = "center",
            legend.key.width     = unit(1.0, "cm"),
            legend.key.height    = unit(0.6, "cm"),
            legend.title         = element_text(hjust = 0.5, size = 12),
            legend.text          = element_text(size = 10),
            legend.box.margin    = margin(t = 0, r = 16, b = 0, l = 8),
            plot.margin          = margin(t = 0, r = 16, b = 0, l = 0),
            text                 = element_text(family = "Times")
        ) +
        guides(fill = guide_legend(
            title.position = "top", title.hjust = 0.5,
            ncol = 1, byrow = FALSE, label.position = "right"
        ))

    cowplot::get_legend(base)
}

# --- Build panels (A with compass; B plain; C plain) --------------------------
p1 <- p10Panel("2000_2023",         show_compass = TRUE)
p2 <- p10Panel("MIROC45_2070_2099", show_compass = FALSE)
p3 <- p10Panel("MIROC85_2070_2099", show_compass = FALSE)

# combine the three panels (equal widths)
row3 <- cowplot::plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1,1,1), align = "h")

# build one vertical legend and attach it to the RIGHT of the combined row
legend_right <- make_vertical_legend("MIROC85_2070_2099")
final_fig <- cowplot::plot_grid(row3, legend_right, ncol = 2,
                                rel_widths = c(1, 0.20), align = "h")

# save (wider canvas so nothing clips)
final_path <- paste0('data_MAIN/5_figs/lowFilter/agg/fig4.png')
ggsave(final_path, plot = final_fig,
       width = 2800, height = 2000, units = "px",
       bg = "white", limitsize = FALSE)

final_fig

# In results, aggregated version of Figure 5 ------------------------------

suitable_hectares_CA = map(
    names,
    function(sp){
        print(sp)
        input.filename.hist = paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/', '/p10_', sp, '_2000_2023_agg.tif')
        input.filename.RCP45 = paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/', '/p10_', sp, '_MIROC45_2070_2099_agg.tif')
        input.filename.RCP85 = paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/', '/p10_', sp, '_MIROC85_2070_2099_agg.tif')
        
        
        
        input.rast.hist = rast(input.filename.hist) %>%
            crop(CA, mask = T)
        input.rast.RCP45 = rast(input.filename.RCP45) %>%
            crop(CA, mask = T)
        input.rast.RCP85 = rast(input.filename.RCP85) %>%
            crop(CA, mask = T)
        
        num_ha.hist = expanse(input.rast.hist, unit='ha', byValue=T)
        num_ha.RCP45 = expanse(input.rast.RCP45, unit='ha', byValue=T)
        num_ha.RCP85 = expanse(input.rast.RCP85, unit='ha', byValue=T)
        
        sp_pretty = sp %>%
            str_to_sentence() %>%
            str_replace("_", ". ")
        
        ha.hist = num_ha.hist$area[num_ha.hist$value==1]
        ha.RCP45 = num_ha.RCP45$area[num_ha.RCP45$value==1]
        ha.RCP85 = num_ha.RCP85$area[num_ha.RCP85$value==1]
        
        percent.hist = round((ha.hist/total_ca_ha)*100, 2)
        percent.RCP45 = round((ha.RCP45/total_ca_ha)*100, 2)
        percent.RCP85 = round((ha.RCP85/total_ca_ha)*100, 2)
        
        out = tibble(
            Species = sp_pretty,
            `Current conditions (hectares)` = round(ha.hist),
            `Future conditions MIROC RCP 4.5 (hectares)` = round(ha.RCP45),
            `Future conditions MIROC RCP 8.5 (hectares)` = round(ha.RCP85),
            `Current conditions (percent)` = percent.hist,
            `Future conditions MIROC RCP 4.5 (percent)` = percent.RCP45,
            `Future conditions MIROC RCP 8.5 (percent)` = percent.RCP85,
            
            `Change under RCP 4.5 (hectares)` =  round(ha.RCP45 - ha.hist),
            `Change under RCP 8.5 (hectares)` =  round(ha.RCP85 - ha.hist),
            `Percent change under RCP 4.5` = round(100*(ha.RCP45 - ha.hist)/(abs(ha.hist)), 2),
            `Percent change under RCP 8.5` = round(100*(ha.RCP85 - ha.hist)/(abs(ha.hist)), 2),
        )
        
        return(out)
    }
) %>%
    bind_rows()

dir.create('data_MAIN/6_tables/lowFilter/agg', recursive=T)
write_csv(suitable_hectares_CA, "data_MAIN/6_tables/lowFilter/agg/suitable_hectares_agg_CA.csv")
suitable_hectares_CA = read_csv("data_MAIN/6_tables/lowFilter/agg/suitable_hectares_agg_CA.csv")

ref_rast = rast('data/0_env/bcm/bcmv8_historic/2000_2023_monthly/aet1999dec.tif')
cv_rast = ref_rast %>%
    crop(central_valley, mask = T)
cv_rast = ifel(is.na(cv_rast), NA, 1) %>%
    as.polygons()
total_cv_ha = cv_rast %>% expanse('ha') #4675858

#In results, aggregated version of Figure 4
suitable_hectares_CV = map(
    names,
    function(sp){
        print(sp)
        input.filename.hist = paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/', '/p10_', sp, '_2000_2023_agg.tif')
        input.filename.RCP45 = paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/', '/p10_', sp, '_MIROC45_2070_2099_agg.tif')
        input.filename.RCP85 = paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/', '/p10_', sp, '_MIROC85_2070_2099_agg.tif')
        
        
        
        input.rast.hist = rast(input.filename.hist) %>%
            crop(cv_rast, mask = T)
        input.rast.RCP45 = rast(input.filename.RCP45) %>%
            crop(cv_rast, mask = T)
        input.rast.RCP85 = rast(input.filename.RCP85) %>%
            crop(cv_rast, mask = T)
        
        num_ha.hist = expanse(input.rast.hist, unit='ha', byValue=T)
        num_ha.RCP45 = expanse(input.rast.RCP45, unit='ha', byValue=T)
        num_ha.RCP85 = expanse(input.rast.RCP85, unit='ha', byValue=T)
        
        sp_pretty = sp %>%
            str_to_sentence() %>%
            str_replace("_", ". ")
        
        ha.hist = num_ha.hist$area[num_ha.hist$value==1]
        ha.RCP45 = num_ha.RCP45$area[num_ha.RCP45$value==1]
        ha.RCP85 = num_ha.RCP85$area[num_ha.RCP85$value==1]
        
        percent.hist = round((ha.hist/total_cv_ha)*100, 2)
        percent.RCP45 = round((ha.RCP45/total_cv_ha)*100, 2)
        percent.RCP85 = round((ha.RCP85/total_cv_ha)*100, 2)
        
        out = tibble(
            Species = sp_pretty,
            `Current conditions (hectares)` = round(ha.hist),
            `Future conditions MIROC RCP 4.5 (hectares)` = round(ha.RCP45),
            `Future conditions MIROC RCP 8.5 (hectares)` = round(ha.RCP85),
            `Current conditions (percent)` = percent.hist,
            `Future conditions MIROC RCP 4.5 (percent)` = percent.RCP45,
            `Future conditions MIROC RCP 8.5 (percent)` = percent.RCP85,
            
            `Change under RCP 4.5 (hectares)` =  round(ha.RCP45 - ha.hist),
            `Change under RCP 8.5 (hectares)` =  round(ha.RCP85 - ha.hist),
            `Percent change under RCP 4.5` = round(100*(ha.RCP45 - ha.hist)/(abs(ha.hist)), 2),
            `Percent change under RCP 8.5` = round(100*(ha.RCP85 - ha.hist)/(abs(ha.hist)), 2),
        )
        
        return(out)
    }
) %>%
    bind_rows()

dir.create('data_MAIN/6_tables/lowFilter/agg')
write_csv(suitable_hectares_CV, "data_MAIN/6_tables/lowFilter/agg/suitable_hectares_agg_CV.csv")

# AUCs table --------------------------------------------------------------
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
aucs = map(
    1:10,
    function(i){
        training = read_csv(paste0('data_MAIN/3_swd/agg/training_', names[i], '_soil200cm_lowFilter_agg.csv'))
        testing = read_csv(paste0('data_MAIN/3_swd/agg/testing_', names[i], '_soil200cm_lowFilter_agg.csv'))
        model = readRDS(paste0('data_MAIN/4_maxent_outputs/agg/', names[i], '/lowFilter/model_training/', names[i], '_training_sdm.rds'))
        best_rm = readRDS(paste0('data_MAIN/4_maxent_outputs/agg/tuning/', names[i], '_finalModelArgs_lowFilter.rds'))[4] %>%
            str_split_i('=', 2)

        training_pred = training %>%
            select(
                x, y, month, year,
                aet, tmx, ppt_winter_sum, tmx_summer_mean, tdiff,
                cec, drclass, om, ph, salinity,
                presence
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
            filter(complete.cases(.)) %>%
            mutate(pred = predict(model, .))

        training.auc = pROC::auc(training_pred$presence, training_pred$pred) %>%
            as.numeric()
        testing.auc = pROC::auc(testing_pred$presence, testing_pred$pred) %>%
            as.numeric()

        out = tibble(
            sp = names[i],
            # training.auc = training.auc,
            testing.auc = round(testing.auc, 3),
            reg_mult = best_rm
        )
        return(out)
    }
) %>%
    bind_rows()

write_csv(aucs, 'data_MAIN/6_tables/lowFilter/agg/aucs_lowFilter.csv')


# Permutation importance --------------------------------------------------

getVImp = function(sp){
    mod.filename = paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
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
        select(sp, everything())
    return(permutation.imp)
}

var_imps = map(
    names,
    getVImp
) %>%
    bind_rows()
write_csv(var_imps, 'data_MAIN/6_tables/lowFilter/agg/perm_var_imps_lowFilter.csv')

sp_pretty_c = names %>%
    str_replace('_', '. ') %>%
    str_to_sentence()
vipLong = function(sp){
    mod.filename = paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/model/', sp, '_final_sdm.rds')
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
            sp = sp_pretty %>%
                factor(levels=sp_pretty_c, ordered=T),
            variable = str_split_i(variable, '.permutation', 1) %>%
                factor(levels=c(
                    'aet', 'tmx', 'ppt_winter_sum', 'tmx_summer_mean', 'tdiff',
                    'cec', 'drclass', 'om', 'ph', 'salinity'
                ), ordered=T)
        )
    return(permutation.imp)
}

var_imps_long = map(
    names, vipLong
) %>%
    bind_rows()

ggplot(var_imps_long) +
    geom_bar(
        aes(x = val, y=variable),
        stat='identity'
    ) +
    facet_wrap(~sp, nrow=2) +
    theme_light() +
    labs(y='', x='', title = 'Variable permutation importance (Ecorelevant model)') +
    theme(
        text=element_text(family='times', size=12),
        strip.text=element_text(face='italic')
    )
ggsave('data_MAIN/5_figs/lowFilter/agg/figS6_varImp.png', width = 3200, height = 2400, units = 'px', bg='white')


# Histogram of species per pixel ------------------------------------------

current_p10_rasts = map(
    names,
    function(sp){
        list.files(
            paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/'),
            pattern = '_2000_2023_agg.tif',
            full.names=T
        )
    }
) %>% 
    unlist() %>%
    rast() %>% 
    sum() %>% 
    crop(central_valley, mask = T)
names(current_p10_rasts) = "Suitable_Current"
suitable_sp_count.current = current_p10_rasts %>% 
    values() %>% 
    table() %>% 
    data.frame()
names(suitable_sp_count.current) = c('n_sp', 'Pixels')
suitable_sp_count.current = suitable_sp_count.current %>% 
    mutate(
        ha = Pixels*0.027,
        Model = "2000 to 2023"
    )

rcp45_p10_rasts = map(
    names,
    function(sp){
        list.files(
            paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/'),
            pattern = '45',
            full.names=T
        )
    }
) %>% 
    unlist() %>% 
    rast() %>% 
    sum() %>% 
    crop(central_valley, mask = T)
names(rcp45_p10_rasts) = "Suitable_RCP4.5"
suitable_sp_count.rcp45 = rcp45_p10_rasts %>% 
    values() %>% 
    table() %>% 
    data.frame()
names(suitable_sp_count.rcp45) = c('n_sp', 'Pixels')
suitable_sp_count.rcp45 = suitable_sp_count.rcp45 %>% 
    mutate(
        ha = Pixels*0.027, 
        Model = "2070 to 2099 (MIROC RCP 4.5)"
    )

rcp85_p10_rasts = map(
    names,
    function(sp){
        list.files(
            paste0('data_MAIN/4_maxent_outputs/agg/', sp, '/lowFilter/p10/'),
            pattern = '85',
            full.names=T
        )
    }
) %>% 
    unlist() %>% 
    rast() %>% 
    sum()
names(rcp85_p10_rasts) = "Suitable_RCP8.5"
suitable_sp_count.rcp85 = rcp85_p10_rasts %>% 
    values() %>% 
    table() %>% 
    data.frame()
names(suitable_sp_count.rcp85) = c('n_sp', 'Pixels')
suitable_sp_count.rcp85 = suitable_sp_count.rcp85 %>% 
    mutate(
        ha = Pixels*0.027, 
        Model = "2070 to 2099 (MIROC RCP 8.5)"
    )

suitable_sp_count.all = suitable_sp_count.current %>% 
    rbind(suitable_sp_count.rcp45) %>% 
    rbind(suitable_sp_count.rcp85)

ggplot(data=suitable_sp_count.all, aes(fill=Model, group=Model)) +
    geom_histogram(
        aes(x=n_sp, y = ha), 
        color = 'darkgrey',
        stat='identity',
        position = 'dodge'
    ) +
    theme_classic() +
    scale_fill_wa_d('flag') +
    scale_y_continuous(labels = scales::comma, limits = c(0, 30000)) +
    labs(
        x = 'Number of Species',
        y = 'Hectares',
        fill = "Years (and Climate Scenario)",
        title = paste0('Area per number of suitable species')
    ) +
    theme(
        text = element_text(
            family="Times New Roman", size = 18
        ), 
        plot.title = element_text(
            family="Times New Roman", size = 24
        )
    )
ggsave("data_MAIN/5_figs/lowFilter/agg/figS4.png", width = 2040, height = 1440, units='px', scale = 2)
suitable_sp_count.all_wider = suitable_sp_count.all %>% 
    select(-Pixels) %>% 
    pivot_wider(names_from='n_sp', values_from=c('ha'), names_prefix = 'N_sp_', values_fill = 0)
write_csv(suitable_sp_count.all_wider, "data_MAIN/6_tables/lowFilter/agg/tableS9_nsp_ha.csv")
