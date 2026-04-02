setwd('/home/tjsipin/will_tj')

source('src_lowFilter_v5_50km/02_spp_occ_background.R')

source('src_lowFilter_v5_50km/R_WY/03_buildModels.R')
source('src_lowFilter_v5_50km/R_WY/04_predictions.R')
source('src_lowFilter_v5_50km/R_WY/05_plotting.R')

source('src_lowFilter_v5_50km/R_agg/03_buildModels.R')
source('src_lowFilter_v5_50km/R_agg/04_predictions.R')
source('src_lowFilter_v5_50km/R_agg/05_plotting.R')

source('src_lowFilter_v5_50km/R_monthly_comprehensive/03_buildModels.R')
source('src_lowFilter_v5_50km/R_monthly_comprehensive/04_predictions.R')
source('src_lowFilter_v5_50km/R_monthly_comprehensive/05_plotting.R')

source('src_lowFilter_v5_50km/R_monthly_selective/03_buildModels.R')
source('src_lowFilter_v5_50km/R_monthly_selective/04_predictions.R')
source('src_lowFilter_v5_50km/R_monthly_selective/05_plotting.R')