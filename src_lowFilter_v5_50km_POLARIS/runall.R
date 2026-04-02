setwd('/home/tjsipin/will_tj')

source('src_lowFilter/02_spp_occ_background.R')

source('src_lowFilter/R_agg/03_buildModels.R')
source('src_lowFilter/R_agg/04_predictions.R')
source('src_lowFilter/R_agg/05_plotting.R')

source('src_lowFilter/R_monthly_comprehensive/03_buildModels.R')
source('src_lowFilter/R_monthly_comprehensive/04_predictions.R')
source('src_lowFilter/R_monthly_comprehensive/05_plotting.R')

source('src_lowFilter/R_monthly_selective/03_buildModels.R')
source('src_lowFilter/R_monthly_selective/04_predictions.R')
source('src_lowFilter/R_monthly_selective/05_plotting.R')