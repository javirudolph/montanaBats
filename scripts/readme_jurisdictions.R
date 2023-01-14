#
# MONTANA JURISDICTIONS SPATIAL DATA
# SCRIPT BY F. JAVIERA RUDOLPH
# QUESTIONS: f.javierarudolph at gmail.com
# LAST UPDATED: January 13, 2023

# Working on this script to show where I downloaded data to set jurisdictions

# libraries
library(sf)
library(dplyr)
library(ggplot2)

# the folder 'MTPublicLands_SHP' has the shapefiles for the different jurisdictions
# I downloaded it from the Montana state library
# link: https://mslservices.mt.gov/Geographic_Information/Data/DataList/datalist_Details.aspx?did=%7B60b5a8b0-b272-11e2-9e96-0800200c9a66%7D

# load the shapefiles

raw_mtpublands <- read_sf("datafiles/jurisdictions/MTPublicLands_SHP/Montana_PublicLands/Montana_PublicLands.shp")

# getting a warning
self_inter <- raw_mtpublands %>% 
  mutate(validity = st_is_valid(raw_mtpublands, reason = TRUE)) %>% 
  filter(validity != "Valid Geometry")

# compare:
self_inter %>% 
  ggplot() + 
  geom_sf()

st_make_valid(self_inter) %>% 
  ggplot() +
  geom_sf()

# make those geometries valid again:
mtpublands <- st_make_valid(raw_mtpublands)

# To visualize independent jurisdictions
owners <- mtpublands$OWNER %>%unique()
# change the value of j to visualize a different jurisdiction
j<-18
mtpublands %>% 
  filter(OWNER == owners[j]) %>% 
  ggplot() +
  geom_sf(fill = "darkslateblue", color = "darkslateblue") +
  labs(title = owners[j])

# Missing the tribal lands
# Found a different shape file:
# https://mslservices.mt.gov/Geographic_Information/Data/DataList/datalist_Details.aspx?did={341205DA-7668-4119-9D21-0D1C8AFCF5F1}

# Map the tribal jurisdictions

mtreservations <- read_sf("datafiles/jurisdictions/MontanaReservations_shp/MontanaReservations.shp")

st_transform(mt_covariates, st_crs(mtreservations)) %>% 
  st_join(., mtreservations) %>%
  ggplot() +
  geom_sf(aes(fill = NAME, color = NAME))



# Put these together in the same thing:
# check they have the same crs
st_crs(mtpublands)
st_crs(mtreservations)

mtpublands %>% 
  select(OWNER, geometry) %>% 
  rename(jurisdiction = OWNER) -> A

# not sure if we need to focus on each independent reservation
# for now, I am not
mtreservations %>% 
  mutate(jurisdiction = "Reservations") %>% 
  select(jurisdiction, geometry) -> B

mt_jurisdictions <- rbind(A, B)
saveRDS(mt_jurisdictions, file = "datafiles/jurisdictions/mt_jurisdictions.RDS")
















