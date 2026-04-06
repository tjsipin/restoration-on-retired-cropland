list.files(
    'data_MAIN/3_swd/monthly/',
    pattern = '_comprehensive',
    full.names = T
) %>%  
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )

list.files(
    'data_v2/3_swd/monthly/comprehensive/',
    pattern = '_comprehensive',
    full.names = T
) %>% 
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("monthly/comprehensive//", "monthly/") %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )

list.files(
    'data_v4_cvBG/3_swd/monthly/comprehensive/',
    pattern = '_comprehensive',
    full.names = T
) %>% 
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("monthly/comprehensive//", "monthly/") %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )

list.files(
    'data_v5_50km_highFilter/3_swd/monthly/comprehensive/',
    pattern = '_comprehensive',
    full.names = T
) %>% 
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("monthly/comprehensive//", "monthly/") %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )

list.files(
    'data_v5_50km_POLARIS/3_swd/monthly/comprehensive/',
    pattern = '_comprehensive',
    full.names = T
) %>% 
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("monthly/comprehensive//", "monthly/") %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )


#####

list.files(
    'data_MAIN/',
    pattern = '_comprehensive',
    full.names = T
) %>%  
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )

list.files(
    'data_v2/',
    pattern = 'comprehensive',
    all.files = T, 
    full.names = T,
    recursive = T,
    include.dirs = T
) %>% 
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("comprehensive", "monthly") %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )

list.files(
    'data_v4_cvBG/',
    pattern = 'comprehensive',
    all.files = T, 
    full.names = T,
    recursive = T,
    include.dirs = T
) %>% 
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("comprehensive", "monthly") %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )

list.files(
    'data_v5_50km_highFilter/4_maxent_outputs/monthly/comprehensive/',
    pattern = '_comprehensive',
    full.names = T
) %>% 
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("monthly/comprehensive//", "monthly/") %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )

list.files(
    'data_v5_50km_POLARIS/4_maxent_outputs/monthly/comprehensive/',
    pattern = '_comprehensive',
    full.names = T
) %>% 
    map(
        function(fi){
            new_filename = fi %>% 
                str_replace("monthly/comprehensive//", "monthly/") %>% 
                str_replace("_comprehensive", "")
            file.rename(fi, new_filename)
        }
    )
