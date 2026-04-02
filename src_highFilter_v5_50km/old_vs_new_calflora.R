calflora_path_new = "data/1_occ/calflora_downloads_030626/"
calflora_path_old = "data/1_occ/calflora_downloads_122625/"


map(
    1:length(names),
    function(i) {
        print(names[i])
        ## Read in data
        input.filename.new = list.files(
            paste0(calflora_path_new),
            pattern = paste0(names[i], '_calflora_download'),
            full.names = T
        )
        input.filename.old = list.files(
            paste0(calflora_path_old),
            pattern = paste0(names[i], '_calflora_download'),
            full.names = T
        )
        df.new <- read_csv(input.filename.new) %>%
            janitor::clean_names()
        df.old <- read_csv(input.filename.old) %>%
            janitor::clean_names()
        
        
        ## filter and keep coords
        df_filter.new <- df.new %>%
            ## "LOW" filter criteria
            filter(location_quality %in% c("high", "medium"),
                   accuracy_square_meters <= 72900,
                   date >= "1999-10-01") 
        df_filter.old <- df.old %>%
            ## "LOW" filter criteria
            filter(location_quality %in% c("high", "medium"),
                   accuracy_square_meters <= 72900,
                   date >= "1999-10-01") 
        
        df_repro_filter.new = df_filter.new %>%
            filter(phenology_code %in% c("F", "R") | is.na(phenology_code)) %>%
            mutate(phenology_code = as.character(phenology_code))
        df_repro_filter.old = df_filter.old %>%
            filter(phenology_code %in% c("F", "R") | is.na(phenology_code)) %>%
            mutate(phenology_code = as.character(phenology_code))
        
        out = tibble(
            species = names[i],
            new = nrow(df_repro_filter.new),
            old = nrow(df_repro_filter.old)
        )
        return(out)
    }
) %>% 
    bind_rows()
