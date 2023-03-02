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

# Descriptions for these covariates and how they were aggregated at the NABat 10x10 cell level are in the report
# https://ecos.fws.gov/ServCat/Reference/Profile/144755
# downloaded report in the data folder for this RProj

# Physiographic diversity on page 8 of the report - it's a measure of landscape complexity: ruggedness
# source: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143619

# Ecoregions could be very useful as they describe habitat features 
# descriptions for ecoregions are here: https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states 
# and the paper: https://link.springer.com/article/10.1007/s00267-014-0364-1
# in the report they mention these can be used to account for spatial autocorrelation
# seems like in their model they use these as the intercept

# karst - grid cells that intersect karst polygons get a 1, https://link.springer.com/article/10.1007/s10040-016-1519-3
# karst indicator was included to capture the effect of caves

# riverlake variable is a 1 if the grid intersects rivers or shorelines

# For detection-level covariates they extracted the data from https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1852
# These depend on the specific night


# change crs to match the jurisdiction data
mt_covs_NAD83 <- st_transform(NABat_covariates, st_crs(mt_jurisdictions))

juris_list <- unique(mt_jurisdictions$jurisdiction)
jusis_list_sub <- juris_list[c(5, 7, 10, 13, 17, 18, 20)]
focus_jurisdictions <- mt_jurisdictions %>% filter(jurisdiction %in% jusis_list_sub)

mt_covariates_w_jurisdictions <- st_join(mt_covs_NAD83, focus_jurisdictions, largest = TRUE)

# Now we are incorporating the three regions we used for the expert elicitation

three_regions_raw <- st_read("datafiles/MT_three_regions_googleEarth.kml")

# make sure coordinate reference systems match
three_regions <- st_transform(three_regions_raw, crs = st_crs(mt_covariates))
st_crs(three_regions)

# join with the covariates dataframe to assign each grid cell to one of the three regions

mt_covariates <- st_join(mt_covariates_w_jurisdictions, three_regions, largest = TRUE) %>% 
  rename(region = Name) %>% 
  select(-Description)

# CHECKS
# Three regions
ggplot(data = mt_covariates) + geom_sf(aes(fill = region))

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
