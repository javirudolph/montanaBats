#
# MONTANA COVARIATES FOR GRID CELLS
# SCRIPT BY F. JAVIERA RUDOLPH
# QUESTIONS: f.javierarudolph at gmail.com
# LAST UPDATED: January 13, 2023


# Building a data frame with information on covariates for each of the grid cells in Montana

# libraries
library(sf)
library(dplyr)
library(ggplot2)

# Join jurisdictions info with the montana covariates and the NABat grids
mt_jurisdictions <- readRDS("datafiles/jurisdictions/mt_jurisdictions.RDS")

# NABAT covariates
# source: https://www.sciencebase.gov/catalog/item/620e6f3bd34e6c7e83baa603

NABat_covariates <- read_sf("datafiles/nabat_covariates/NABat_grid_covariates/NABat_grid_covariates.shp") %>%
  filter(., admin1 == "Montana") %>% 
  select(GRTS_ID, karst, p_forest, p_wetland, mean_temp,
         precip, DEM_max, physio_div, dist_mines, starts_with("eco"))

# change crs to match the jurisdiction data
mt_covs_NAD83 <- st_transform(NABat_covariates, st_crs(mt_jurisdictions))

juris_list <- unique(mt_jurisdictions$jurisdiction)
jusis_list_sub <- juris_list[c(5, 7, 10, 13, 17, 18, 20)]
focus_jurisdictions <- mt_jurisdictions %>% filter(jurisdiction %in% jusis_list_sub)

mt_covariates <- st_join(mt_covs_NAD83, focus_jurisdictions, largest = TRUE)

# CHECKS
# ecoregions

ggplot(data = mt_covariates) + geom_sf(aes(fill = eco1_name))
ggplot(data = mt_covariates) + geom_sf(aes(fill = eco2_name))
ggplot(data = mt_covariates) + geom_sf(aes(fill = eco3_name))

# jurisdictions
ggplot(data = mt_covariates) + geom_sf(aes(fill = jurisdiction))

# environmental variables
ggplot(data = mt_covariates) + geom_sf(aes(fill = as.factor(karst)))
ggplot(data = mt_covariates) + geom_sf(aes(fill = p_forest))
ggplot(data = mt_covariates) + geom_sf(aes(fill = p_wetland))
ggplot(data = mt_covariates) + geom_sf(aes(fill = mean_temp))
ggplot(data = mt_covariates) + geom_sf(aes(fill = precip))
ggplot(data = mt_covariates) + geom_sf(aes(fill = DEM_max))
ggplot(data = mt_covariates) + geom_sf(aes(fill = physio_div))
ggplot(data = mt_covariates) + geom_sf(aes(fill = dist_mines))

saveRDS(mt_covariates, file = "datafiles/mt_covariates.RDS")

