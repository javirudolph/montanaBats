# testing simulation for Montana


# Libraries
library(tidyverse)
library(sf)

# Things that I need:
# geometry for the state with polygons for the grid cells

# NABAT covariates
# source: https://www.sciencebase.gov/catalog/item/620e6f3bd34e6c7e83baa603

mt_covariates <- read_sf("datafiles/nabat_covariates/NABat_grid_covariates/NABat_grid_covariates.shp") %>%
  filter(., admin1 == "Montana") %>% 
  select(GRTS_ID, karst, p_forest, p_wetland, mean_temp,
         precip, DEM_max, physio_div, dist_mines, starts_with("eco"))


# Visualize the grids
mt_covariates %>% 
  ggplot() +
  geom_sf()

# Get state outline
library(maps)
mt_outline <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>%
  filter(ID == "montana")

mt_covariates %>% 
  ggplot() +
  geom_sf(data = mt_outline) +
  geom_sf()

# Map the tribal jurisdictions

mtreservations <- read_sf("datafiles/jurisdictions/MontanaReservations_shp/MontanaReservations.shp")

st_transform(mt_covariates, st_crs(mtreservations)) %>% 
  st_join(., mtreservations) %>%
  ggplot() +
  geom_sf(aes(fill = NAME, color = NAME))

# Map the public lands

mtpublands <- read_sf("datafiles/jurisdictions/MTPublicLands_SHP/Montana_PublicLands/Montana_PublicLands.shp")

# getting a warning
self_inter <- mtpublands %>% 
  mutate(validity = st_is_valid(mtpublands, reason = TRUE)) %>% 
  filter(validity != "Valid Geometry")

self_inter %>% 
  ggplot() + 
  geom_sf()

st_make_valid(self_inter) %>% 
  ggplot() +
  geom_sf()

# make those geometries valid again:
valid_mtpublands <- st_make_valid(mtpublands)

st_transform(mt_covariates, st_crs(valid_mtpublands)) %>% 
  st_join(., valid_mtpublands) %>% 
  ggplot() +
  geom_sf(aes(fill = OWNER, color = OWNER))


# Put these together in the same thing:

valid_mtpublands %>% 
  select(OWNER, geometry) %>% 
  rename(jurisdiction = OWNER) -> A

mtreservations %>% 
  mutate(jurisdiction = "Reservations") %>% 
  select(jurisdiction, geometry) -> B

juris_geom <- rbind(A,B)
st_crs(juris_geom)

mt_covs_NAD83 <- st_transform(mt_covariates, st_crs(juris_geom))

juris_plus_covs <- st_join(mt_covs_NAD83, juris_geom)


juris_plus_covs %>% 
  filter(is.na(jurisdiction)) -> test

test %>% 
  ggplot() + 
  geom_sf()

# Jurisdictions we are focused on and working with:
unique(juris_plus_covs$jurisdiction)
focus_juris <- c("Montana Fish, Wildlife, and Parks", "US Bureau of Land Management", "US Forest Service",
                 "Reservations", "US Fish and Wildlife Service", "National Park Service")

focus_jurs_covs <- st_join(mt_covs_NAD83, juris_geom %>% 
          filter(jurisdiction %in% focus_juris))

focus_jurs_covs %>% 
  ggplot() + 
  geom_sf(aes(color = jurisdiction, fill = jurisdiction))

