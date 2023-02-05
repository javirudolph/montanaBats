# FORMATTING BAT CALL DATA AND FITTING MULTI-STATE OCCUPANCY MODELS
# Questions: f.javierarudolph@gmail.com

# Three data products were shared by Christian Stratton and the description of those is as follows:
#
# > *SiteDataClean* are data on sites including site type, clutter distance, and water distance. Also the detector type and microphone are listed. We did not get 123 data for some sites. I manually classified them based on the site description, but didn't assign a distance to water or clutter. I also included the serial numbers for the detector and microphone. Not sure if these are useful. If you are going to use them I can clean them up as there is some variability in how they were recorded between observers.
#
# > *AllCallsCapHistory* is a summary of all calls recorded by cell, site, and year and Julian day (sample night). NA values indicate that he detector did not record on that day. 0's are no calls recorded but the detector was deployed.
#
# > *WNSCallsCapHistory* is a summary of all calls that have a mean characteristic frequency at or above 34kHz.



# 1. Libraries --------------------------------------------------
library(tidyverse)
library(jagsUI)
library(sf)
library(AHMbook) # for the standardize() function


# 2. Call Data  ---------------------------

# Read the raw files - focus on WNS susceptible species only (WNSCallsCapHistory.csv)
# Data is in wide format with columns = julian date
# CSY = cell site year, where cell corresponds to the NABat Cell ID

raw_wns_calls <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/WNSCallsCapHistory.csv", show_col_types = FALSE)

# What do we consider the cutoff for many bats?
manybats_cut <- 250

# Transform to long format:
raw_wns_calls %>%
  # create variable for julian date, formerly in the columns
  pivot_longer(., 4:last_col(), names_to = "jdate", values_to = "y") %>%
  # keep only needed variables
  select(SiteID, CSY, jdate, y) %>%
  # julian date should be a number
  mutate(jdate = as.double(jdate)) %>%
  # remove NAs - means recorders wasn't on
  drop_na(y) %>%
  # clean the names
  janitor::clean_names() %>%
  # separate CSY into three
  separate(csy, c("cell", "site", "year"), remove = FALSE) %>%
  # create binary version of detection data
  mutate(dnd = as.numeric(y > 0)) %>%
  # Create state categories, y multi state
  # Asumming a cutoff of 100 bat calls as many bats
  mutate(yms = case_when(
    y == 0 ~ 1,
    y <= manybats_cut ~ 2,
    y > manybats_cut ~ 3
  )) %>%
  # add column for number of nights of recording for each site_id
  add_count(site_id, name = "n_surveys") -> wns_calls_long

wns_calls_long %>% 
  ggplot(aes(x = y)) +
  geom_histogram()

summary(wns_calls_long$y)

wns_calls_long %>% 
  ggplot(aes(x = yms)) +
  geom_histogram()

# We have 4 sites in each cell
# Except cell 1313 which shows as having sites 1,2,3,4 in year 2021, and sites 1,3,4,5 in year 2020. 
# These sites are not always the same location, should consider changing that 5 to a 2

##2a. Quantities --------------


ncells <- length(unique(wns_calls_long$cell))
nsites <- max(wns_calls_long$site)
max_nsurveys <- as.numeric(max(wns_calls_long$n_surveys))
nyears <- length(unique(wns_calls_long$year))
minyear <- min(wns_calls_long$year)

siteidlist <- sort(unique(wns_calls_long$site_id))
sitelist <- sort(unique(wns_calls_long$site))
celllist <- sort(unique(wns_calls_long$cell))
yearlist <- sort(unique(wns_calls_long$year))

#3. Site Covariates ------------------

raw_sitecovs <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/SiteDataCleanFinal.csv")

raw_sitecovs %>% 
  janitor::clean_names() %>% 
  select(csy, site_id, jstart, jstop, site_type, other_site_type) %>% 
  mutate(site_type = str_remove(site_type, " "),
         jduration = jstop-jstart) -> sitecovs
  # separate(csy, c("cell", "site", "year"), remove = FALSE)


##3a. Standardize dates ---------

hist(sitecovs$jstart)
summary(sitecovs$jstart)
max(sitecovs$jstart) - min(sitecovs$jstart)
sd(sitecovs$jstart)

st_jstart <- (sitecovs$jstart-median(sitecovs$jstart))/sd(sitecovs$jstart)
hist(st_jstart)
sd(st_jstart)

sitecovs %>% 
  mutate(date = (jstart-median(jstart))/sd(jstart),
         date_sqrd = date^2) -> sitecovs

##3b. Factor sitetypes ----------

sitetype_list <- unique(sitecovs$site_type)

# we should try to classify these into the categories available
sitecovs %>% 
  filter(site_type == "Other")

sitecovs %>% 
  mutate(fctr_sitetype = as.numeric(factor(site_type, levels = sitetype_list))) -> sitecovs



#4. Cell Covariates -----------------

mt_covariates <- read_sf("datafiles/nabat_covariates/NABat_grid_covariates/NABat_grid_covariates.shp") %>%
  filter(., admin1 == "Montana")

mt_covariates %>% 
  filter(GRTS_ID %in% celllist) %>% 
  select(GRTS_ID, karst, p_forest, p_wetland, mean_temp,
         precip, DEM_max, physio_div, dist_mines, starts_with("eco")) %>% 
  # arrange to make sure its same order of cells as obs data
  arrange(factor(GRTS_ID, levels = celllist)) %>% 
  rename(cell = GRTS_ID) -> cell_covs



#5. Build Arrays --------------------

##5a. MultiState data ---------
# re-organize the state data for bat calls into an array for the model

yms <- array(NA, dim = c(ncells, nsites, max_nsurveys, nyears),
             dimnames = list(celllist, 1:nsites, 1:max_nsurveys, yearlist))

for(i in 1:ncells){
  for(j in 1:nsites){
    for(t in 1:nyears){
      sel_csy <- paste(celllist[i], sitelist[j], yearlist[t], sep = "_")
      tmp <- wns_calls_long[wns_calls_long$csy == sel_csy,]
      nr <- nrow(tmp)
      
      if(nr > 0){
        yms[i, j, 1:nr, t] <- tmp$yms
      }
    }
  }
}

# Example - First survey in 2020 for all cells with their sites
# order of dimensions: cell, site, survey, year
yms[ , , 1 , 1]

##5b. Site covs -------

date <- sitetype <- duration <- array(NA, dim = c(ncells, nsites, nyears),
                dimnames = list(celllist, 1:nsites, yearlist))

for(i in 1:ncells){
  for(j in 1:nsites){
    for(t in 1:nyears){
      sel_csy <- paste(celllist[i], sitelist[j], yearlist[t], sep = "_")
      tmp <- sitecovs[sitecovs$csy == sel_csy,]
      nr <- nrow(tmp)
      
      if(nr>0){
        date[i,j,t] <- tmp$date
        duration[i,j,t] <- tmp$jduration
        sitetype[i,j,t] <- tmp$fctr_sitetype
      }
    }
  }
}

# Example - some cells were not sampled that year = NA in row
date
duration
sitetype



##5c. Cell covs ---------

summary(elev <- cell_covs$DEM_max)
elev.scaled <- standardize(elev)
summary(forest <- cell_covs$p_forest)
forest.scaled <- standardize(forest)
summary(temp <- cell_covs$mean_temp) # wondering if this is already scaled or the units
temp.scaled <- standardize(temp)
summary(precip <- cell_covs$precip)
precip.scaled <- standardize(precip)
summary(wetlands <- cell_covs$p_forest)
wetlands.scaled <- standardize(wetlands)
summary(physdiv <- cell_covs$physio_div)
physdiv.scaled <- standardize(physdiv)
karst <- cell_covs$karst

# Not sure how to handle region yet
unique(cell_covs$eco3_name)
region <- cell_covs$eco3_name %>% 
  as_factor() %>% 
  as.numeric()
regionID <- cell_covs$eco3_name
cbind(region, regionID)

##5d. Nsurveys ------
# There is variation in the number of nights that microphones are recording for each cell or site


nsurveys <- array(NA, dim = c(ncells, nsites, nyears))
for(i in 1:ncells){
  for(j in 1:nsites){
    for(t in 1:nyears){
      tmp <- which(!is.na(yms[i,j,,t]))
      if(length(tmp) > 0){
        nsurveys[i,j,t] <- max(tmp)
      }
    }
  }
}
nsurveys

#6. Bundle data -----------

str(batcalldata <- list(y = yms, 
                      ncell = dim(yms)[1], nsite = dim(yms)[2], 
                      nsurveys = nsurveys, nyears = dim(yms)[4],
                      date = date,
                      duration = duration,
                      sitetype = sitetype,
                      region = region,
                      elev = elev.scaled, 
                      temp = temp.scaled,
                      physdiv = physdiv.scaled,
                      precip = precip.scaled,
                      forest = forest.scaled,
                      wetlands = wetlands.scaled,
                      karst = karst))
str(batcalldata)

#7. Specify model -----------------

cat(file = "jags_txt/msms_occu.txt", "
model {

  # Definition of state vector (Omega) and observation matrix (Theta)
  # State vector (Omega)
  for (i in 1:ncell){
    for (t in 1:nyears){
      Omega[i,t,1] <- 1 - psi[i,t]           # Prob. of no bats
      Omega[i,t,2] <- psi[i,t] * (1-r[i,t])  # Prob. of occ. by a few bats
      Omega[i,t,3] <- psi[i,t] * r[i,t]      # Prob. of occ. by many bats
    }
  }
  # Observation matrix (Theta)
  # Order of indices: true state, observed state, site, occasion, year
  for(i in 1:ncell){
    for(j in 1:nsites){
      for (t in 1:nyears){
        for (k in 1:nsurveys[i,j,t]){
          Theta[1,1,i,j,t,k] <- 1
          Theta[1,2,i,j,t,k] <- 0
          Theta[1,3,i,j,t] <- 0
          Theta[2,1,i,j,t] <- 1-p2[i,j,t]
          Theta[2,2,i,j,t] <- p2[i,j,t]
          Theta[2,3,i,j,t] <- 0
          Theta[3,1,i,j,t] <- 1-p3[2,i,j,t]-p3[3,i,j,t]
          Theta[3,2,i,j,t] <- p3[2,i,j,t]
          Theta[3,3,i,j,t] <- p3[3,i,j,t]
      }
    }
    }
  }
  # State-space likelihood
  # State equation: model of true states (z)
  for (i in 1:nsites){
    for (t in 1:nyears){
      z[i,t] ~ dcat(Omega[i,t,])
    }
  }
  # Observation equation: model for observed multistate detections
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        y[i,j,t] ~ dcat(Theta[z[i, t], ,i,j,t])
      }
    }
  }
}
")















