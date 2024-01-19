library(sf)
library(soilDB)
library(mapview)
library(terra)
library(tidyverse)

test_ws <- st_read("C:/Users/laura/OneDrive/Documents/000UGA/Research/5-N-EWN/10-Chapter 3/dam-removal/1-input-data/1-geospatial/watersheds/veaziewatershed.shp")
mapview(test_ws)

bb <- sf::st_bbox(test_ws)

mid_x <- (bb["xmin"] + bb["xmax"]) / 2
mid_y <- (bb["ymin"] + bb["ymax"]) / 2

bbox1 <- st_bbox(c(xmin = bb["xmin"] %>% as.numeric(), 
                   xmax = mid_x %>% as.numeric(), 
                   ymin = bb["ymin"]%>% as.numeric(), 
                   ymax = mid_y%>% as.numeric()), 
                 crs = st_crs(test_ws))

bbox2 <- st_bbox(c(xmin = mid_x %>% as.numeric(), 
                            xmax = bb["xmax"] %>% as.numeric(), 
                            ymin = bb["ymin"]%>% as.numeric(), 
                            ymax = mid_y%>% as.numeric()), 
                          crs = st_crs(test_ws))

bbox3 <- st_bbox(c(xmin = mid_x %>% as.numeric(), 
                   xmax = bb["xmax"] %>% as.numeric(), 
                   ymin =mid_y%>% as.numeric(), 
                   ymax = bb["ymax"]%>% as.numeric()), 
                 crs = st_crs(test_ws))

bbox4 <- st_bbox(c(xmin = bb["xmin"] %>% as.numeric(), 
                   xmax = mid_x %>% as.numeric(), 
                   ymin = mid_y%>% as.numeric(), 
                   ymax = bb["ymax"]%>% as.numeric()), 
                 crs = st_crs(test_ws))
mapview(test_ws)+
mapview(bbox1)+
  mapview(bbox2)+
  mapview(bbox3)+
  mapview(bbox4)

bbox_ls <- list(bbox1, bbox2, bbox3, bbox4)
top_extract_all <- NULL
top_comp_all <- NULL

for(i in 1:length(bbox_ls)){
  soil <- try((soilDB::mukey.wcs(aoi = bbox_ls[[i]],
                                 db = 'gssurgo',
                                 quiet = TRUE)))
  
  # should build a chunking method for this
  if(class(soil) == 'try-error'){
    
    this_var_tib <- tibble(site_code = site,
                           year = NA,
                           var =  names(nrcs_var_name),
                           val = NA)
    
    return(generate_ms_exception(glue('{s} is too large',
                                      s = site)))
  }
  
  mukey_values <- terra::as.data.frame(soil) %>% distinct() %>% pull(mukey)
  
  # dbovendry_r g/cm3, om_r %
  nrcs_var_name <- c("om_r", "dbovendry_r", "kffact", "kwfact")
  #### Grab soil vars
  
  # Query Soil Data Acess (SDA) to get the component key (cokey) in each mukey.
  #### each map unit made up of componets but components do not have
  #### spatial informaiton associted with them, but they do include information
  #### on the percentage of each mukey that is made up of each component.
  #### This informaiton on compositon is held in the component table in the
  #### comppct_r column, givin in a percent
  
  mukey_sql <- soilDB::format_SQL_in_statement(mukey_values)
  component_sql <- sprintf("SELECT cokey, mukey, compname, comppct_r, majcompflag FROM component WHERE mukey IN %s", mukey_sql)
  component <- soilDB::SDA_query(component_sql)
  
  # Check is soil data is available
  if(length(unique(component$compname)) == 1 &&
     unique(component$compname) == 'NOTCOM'){
    
    this_var_tib <- tibble(site_code = site,
                           year = NA,
                           var =  names(nrcs_var_name),
                           val = NA)
    
    return(generate_ms_exception(glue('No data was retrived for {s}',
                                      s = site)))
  }
  
  cokey <- unique(component[,1])
  cokey <- data.frame(ID=cokey)
  
  # Query SDA for the componets to get infromation on their horizons from the
  #### chorizon table. Each componet is made up of soil horizons (layer of soil
  #### vertically) identified by a chkey.
  #### Informaiton about the depth of each horizon is needed to calculate weighted
  #### averges of any parameter for the whole soil column
  #### hzdept_r = depth to top of horizon
  #### hzdepb_r = depth to bottom of horizon
  #### om_r = percent organic matter
  cokey_sql <- soilDB::format_SQL_in_statement(cokey$ID)
  chorizon_sql <- sprintf(paste0('SELECT cokey, chkey, hzname, desgnmaster, hzdept_r, hzdepb_r, ',
                                 paste(nrcs_var_name, collapse = ', '),
                                 ' FROM chorizon WHERE cokey IN %s'), cokey_sql)
  
  full_soil_data <- soilDB::SDA_query(chorizon_sql)
  
  # as I understand it, a chkey is a horizon (a soil layer)
  # a cokey is a component, a group of soil layers
  # a mukey is a map unit, a geographic area
  
  # this calculates the carbon in g/cm3 in the top 5 cm of soil by component
  top_extract <- full_soil_data %>% 
    mutate(kffact = as.numeric(kffact), 
           kwfact = as.numeric(kwfact)) %>% 
    filter(hzdept_r < 5) %>% 
    mutate(depth_calc = if_else(hzdept_r == 0, if_else(hzdepb_r < 5, hzdepb_r, 5), 5-hzdept_r), 
           carb_prop = om_r * 10^-2 * dbovendry_r * (depth_calc/5), 
           kw_prop = kwfact* (depth_calc/5), 
           kf_prop = kffact * (depth_calc/5)) %>% 
    group_by(cokey) %>% 
    summarize(carb = sum(carb_prop), 
              kw = sum(kw_prop),
              kf = sum(kf_prop)) %>% 
    mutate(bbox = i)
  
  # this calculates a weighted average of components by map unit
  top_comp <- top_extract %>% full_join(component, by = "cokey") %>% 
    group_by(mukey) %>%
    mutate(comppct_r_sum = sum(comppct_r)) %>%
    summarise(carb_value_mukey = sum(carb*(comppct_r/comppct_r_sum)), 
              kw_value_mukey = sum(kw*(comppct_r/comppct_r_sum)), 
              kf_value_mukey = sum(kf*(comppct_r/comppct_r_sum))) %>% 
    mutate(bbox = i)
  
  
  # Calculate weighted average for the entire soil column. This involves 3 steps.
  #### First the component's weighted average of all horizones. Second, the
  #### weighted avergae of each component in a map unit. And third, the weighted
  #### average of all mukeys in the watershed (weighted by their area)
  
  
  # Watershed weighted average
  site_boundary_p <- sf::st_transform(test_ws, crs = sf::st_crs(soil))
  
  soil_masked <- terra::mask(soil, site_boundary_p)
  
  assign(paste0("soil_masked", i), soil_masked)
  top_extract_all <- top_extract_all %>% bind_rows(top_extract)
  top_comp_all <- top_comp_all %>% bind_rows(top_comp)
  
}

# something funky is happening with the merge process so I am going to proceed to see if I actually need to merge
rlist <- list(soil_masked1, soil_masked2, soil_masked3, soil_masked4)
rsrc <- sprc(rlist)

watershed_mukey_values <- NULL

for(i in 1:length(rlist)){
  temp <- terra::as.data.frame(rlist[[i]]) %>% 
    filter(!is.na(mukey)) %>%
    group_by(mukey) %>%
    summarise(n = n())
  
  watershed_mukey_values <- watershed_mukey_values %>% bind_rows(temp)
}
  
# I think you have to sum after left joining and removing NAs
watershed_mukey_sum <- watershed_mukey_values %>% group_by(mukey) %>% summarize(count = sum(n)) 

carb_mukey <- watershed_mukey_sum %>% mutate(mukey_ch = as.character(mukey)) %>% 
  left_join(top_comp_all %>% mutate(mukey_ch = as.character(mukey)) %>% dplyr::select(mukey_ch, carb_value_mukey) %>% distinct(), by = "mukey_ch") %>% 
  filter(is.na(carb_value_mukey)==FALSE) %>% 
  mutate(total = sum(count), prop = count/total) %>% 
  mutate(carb_w = prop * carb_value_mukey)
  
kw_mukey <- watershed_mukey_sum %>% mutate(mukey_ch = as.character(mukey)) %>% 
  left_join(top_comp_all %>% mutate(mukey_ch = as.character(mukey)) %>% dplyr::select(mukey_ch, kw_value_mukey) %>% distinct(), by = "mukey_ch") %>% 
  filter(is.na(kw_value_mukey)==FALSE) %>% 
  mutate(total = sum(count), prop = count/total) %>% 
  mutate(kw_w = prop * kw_value_mukey)

kf_mukey <- watershed_mukey_sum %>% mutate(mukey_ch = as.character(mukey)) %>% 
  left_join(top_comp_all %>% mutate(mukey_ch = as.character(mukey)) %>% dplyr::select(mukey_ch, kf_value_mukey) %>% distinct(), by = "mukey_ch") %>% 
  filter(is.na(kf_value_mukey)==FALSE) %>% 
  mutate(total = sum(count), prop = count/total) %>% 
  mutate(kf_w = prop * kf_value_mukey)

watershed_mukey <- watershed_mukey_sum %>% mutate(mukey_ch = as.character(mukey)) %>% 
  left_join(top_comp_all %>% mutate(mukey_ch = as.character(mukey)) %>% dplyr::select(mukey_ch, carb_value_mukey, kw_value_mukey, kf_value_mukey) %>% distinct(), by = "mukey_ch")

watershed_carb <- sum(carb_mukey$carb_w, na.rm = T)
watershed_kw <- sum(kw_mukey$kw_w, na.rm = T)
watershed_kf <- sum(kf_mukey$kf_w, na.rm = T)

print(paste("soil OM (g/m2):", watershed_carb * 5 * 10^4))
print(paste("kw:", watershed_kw))
print(paste("kf:", watershed_kf))


# 8% is NA for carbon
carb_na_mukey <- watershed_mukey %>% filter(is.na(carb_value_mukey))
print(paste("soil OM NA percent:", (sum(carb_na_mukey$count)/29302684)*100))

kw_na_mukey <- watershed_mukey %>% filter(is.na(kw_value_mukey))
print(paste("kw NA percent:", (sum(kw_na_mukey$count)/29302684)*100))

kf_na_mukey <- watershed_mukey %>% filter(is.na(kf_value_mukey))
print(paste("kf NA percent:", (sum(kf_na_mukey$count)/29302684)*100))

out_df <- data.frame(soc = watershed_carb * 5 * 10^4, kw = watershed_kw, kf = watershed_kf)
write.csv(out_df, "1-input-data/6-clow-data/veazie_soil.csv", row.names = FALSE)
