# Working on this script to show where I downloaded data to set jurisdictions

# the folder 'MTPublicLands_SHP' has the shapefiles for the different jurisdictions
# I downloaded it from the Montana state library
# link: https://mslservices.mt.gov/Geographic_Information/Data/DataList/datalist_Details.aspx?did=%7B60b5a8b0-b272-11e2-9e96-0800200c9a66%7D


# libraries
library(sf)
library(dplyr)
library(ggplot2)
library(maps)

# load the shapefiles

mtpublands <- st_read("datafiles/jurisdictions/MTPublicLands_SHP/Montana_PublicLands/Montana_PublicLands.shp")

# Get state outline
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states_mt <- states %>% filter(ID == "montana")

# look at all of them. heavy plot, takes a while
mtpublands %>% 
  ggplot() +
  geom_sf(data = states_mt) +
  geom_sf(aes(fill = OWNER, color = OWNER))

# To visualize independent jurisdictions
owners <- mtpublands$OWNER %>%unique()
# change the value of j to visualize a different jurisdiction
j<-18
mtpublands %>% 
  filter(OWNER == owners[j]) %>% 
  ggplot() +
  geom_sf(data = states_mt) +
  geom_sf(fill = "darkslateblue", color = "darkslateblue") +
  labs(title = owners[j])

# Missing the tribal lands
# Found a different shape file:
# https://mslservices.mt.gov/Geographic_Information/Data/DataList/datalist_Details.aspx?did={341205DA-7668-4119-9D21-0D1C8AFCF5F1}

mtreservations <- read_sf("datafiles/jurisdictions/MontanaReservations_shp/MontanaReservations.shp")

mtreservations %>% 
  ggplot() +
  geom_sf(data = states_mt) +
  geom_sf(aes(fill = NAME, color = NAME))


# Create a single dataframe with the geometries of interest
our_jurisd <- owners[c(5, 7, 10, 13, 17, 18)]

mtpublands %>% 
  select(OWNER, geometry) %>% 
  janitor::clean_names() %>% 
  filter(owner %in% our_jurisd) %>% 
  rename(jurisd = owner) -> A

mtreservations %>% 
  select(geometry) %>% 
  rename(jurisd = 'Reservations') -> B

bind_rows(A, B) -> C

C %>% 
  ggplot() +
  geom_sf(data = states_mt) +
  geom_sf(aes(fill = jurisd, color = jurisd))
