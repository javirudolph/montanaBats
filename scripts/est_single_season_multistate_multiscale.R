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

# 2. Load bat call data csv and modify ---------------------------

# Read the raw files - focus on WNS susceptible species only (WNSCallsCapHistory.csv)
# Data is in wide format with columns = julian date
# CSY = cell site year, where cell corresponds to the NABat Cell ID

raw_wns_calls <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/WNSCallsCapHistory.csv", show_col_types = FALSE)

# What do we consider the cutoff for many bats?
manybats_cut <- 100

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


# The spatial scale we decided to work on was the 10x10km cell grid from NABat
# The data presented here has four sites sampled within a cell grid: multi-scale
# There is a varying number of nights surveyed for each site


# create id column 
wns_calls_long %>% 
  mutate(cell_site = paste(cell, site, sep = "_")) %>% 
  filter(year == 2020) -> obs

nsites <- as.numeric(max(obs$site))
max_nsurveys <- as.numeric(max(obs$n_surveys))
ncell <- as.numeric(length(unique(obs$cell)))

siteidlist <- sort(unique(obs$site_id))
sitelist <- sort(unique(obs$site))
celllist <- sort(unique(obs$cell))


# Create empty arrays and fill them with the data
yms <- date <- y <- array(NA, dim = c(ncell, nsites, max_nsurveys),
                          dimnames = list(celllist, 1:nsites, 1:max_nsurveys))

for(i in 1:ncell){
  for(j in 1:nsites){
    sel_cell_site <- paste(celllist[i], sitelist[j], sep = "_")
    tmp <- obs[obs$cell_site == sel_cell_site,]
    nr <- nrow(tmp)
    if(nr > 0){
      yms[i, j, 1:nr] <- tmp$yms
      date[i, j, 1:nr] <- tmp$jdate
      y[i, j, 1:nr] <- tmp$dnd
    }
  }
}




## NABat covariates ---------------------------------------
# source: https://www.sciencebase.gov/catalog/item/620e6f3bd34e6c7e83baa603

mt_covariates <- read_sf("datafiles/nabat_covariates/NABat_grid_covariates/NABat_grid_covariates.shp") %>%
  filter(., admin1 == "Montana")

mt_covariates %>% 
  filter(GRTS_ID %in% celllist) %>% 
  select(GRTS_ID, karst, p_forest, p_wetland, mean_temp,
         precip, DEM_max, physio_div, dist_mines, starts_with("eco")) %>% 
  # arrange to make sure its same order of cells as obs data
  arrange(factor(GRTS_ID, levels = sitelist)) -> cell_covs


library(AHMbook) # for the standardize() function

# fit model with a couple of covariates only and in the parameters for year 1

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

region <- cell_covs$eco3_name %>% 
  as_factor() %>% 
  as.numeric()

regionID <- cell_covs$eco3_name

cbind(region, regionID)

# Bundle data -------------------------------

# There is variation in the number of nights that microphones are recording for each cell or site

nsurveys <- array(NA, dim = c(ncell, nsites))
for(i in 1:ncell){
  for(j in 1:nsites){
    tmp <- which(!is.na(yms[i,j,]))
    if(length(tmp) > 0){
      nsurveys[i,t] <- max(tmp)
    }
  }
}

str(bcalldata <- list(y = yms, ncell = dim(yms)[1], nsite = dim(yms)[2], nsurveys = dim(yms)[3], 
                      region = region,
                      elev = elev.scaled, 
                      temp = temp.scaled,
                      physdiv = physdiv.scaled,
                      precip = precip.scaled,
                      forest = forest.scaled,
                      wetlands = wetlands.scaled,
                      karst = karst))

# Specify model----------------------------
cat(file = "ss_ms_ms.txt", "
model {

  ### (1) Priors
  ## Omega
  for (i in 1:ncell){
    logit(psi[i]) <- alpha.lpsi[region[i]] + beta.lpsi[1] * elev[i] + beta.lpsi[2] * elev[i]^2 + beta.lpsi[3] * temp[i] + beta.lpsi[4] * temp[i]^2 + beta.lpsi[5] * physdiv[i] + beta.lpsi[6] * precip[i] + beta.lpsi[7] * forest[i] + beta.lpsi[8] * wetlands[i]
    logit(r[i]) <- alpha.lr[region[i]] + beta.lr[1] * karst[i] + beta.lr[2] * forest[i] + beta.lr[3] * physdiv[i]
  }
  
  # Priors for parameters in the linear models of psi and r
  # Region-specific intercepts
  for (k in 1:6){
    alpha.lpsi[k] <- logit(mean.psi[k])
    mean.psi[k] ~ dunif(0, 1)
    alpha.lr[k] <- logit(mean.r[k])
    mean.r[k] ~ dunif(0, 1)
  }
  
  # Priors for coefficients of covariates in Omega
  for (k in 1:8){
    beta.lpsi[k] ~ dnorm(0, 0.1)
  }
  
  for (k in 1:3){
    beta.lr[k] ~ dnorm(0, 0.1)
  }
  
  ## Theta 
  # Priors for parameters in local occupancy
  theta2 ~ dunif(0, 1)              # Local occupancy when cell has few bats
  for (s in 1:3) {                  # Local occupancy when cell has many bats
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    theta3[s] <- beta[s] / sum(beta[])
  }
  
  ## detP
  # Priors for parameters in observation process
  p2 ~ dunif(0, 1)              # detection with few bats locally
  for (s in 1:3) {                  # detection with many bats locally
    alpha[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- alpha[s] / sum(alpha[])
  }
  
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector
  for (i in 1:ncell){
    Omega[i,1] <- 1 - psi[i]                 # Prob. of no bats
    Omega[i,2] <- psi[i] * (1-r[i])          # Prob. of occ. by a few bats
    Omega[i,3] <- psi[i] * r[i]              # Prob. of occ. by many bats
  }
  
  # Define local occupancy probability matrix (Theta)
  # Order of indices: global state, local state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta2
  Theta[2,2] <- p2
  Theta[2,3] <- 0
  Theta[3,1] <- theta3[1]
  Theta[3,2] <- theta3[2]
  Theta[3,3] <- theta3[3]
  
  # Define observation probability matrix (detP)
  # Order of indices: true local state, observed local state
  detP[1,1] <- 1
  detP[1,2] <- 0.0000001
  detP[1,3] <- 0.0000001
  detP[2,1] <- 1-p2
  detP[2,2] <- p2
  detP[2,3] <- 0.0000001
  detP[3,1] <- p3[1]
  detP[3,2] <- p3[2]
  detP[3,3] <- p3[3]
  
  
  ### (3) Likelihood
  
  # global occupancy
  for (i in 1:ncell){
    z[i] ~ dcat(Omega[i,])
  }
  
  # local occupancy
  for (i in 1:ncell){
    for (j in 1:nsite){
      u[i,j] ~ dcat(Theta[z[i], ])
    }
  }
  
  # detection
  for (i in 1:ncell){
    for (j in 1:nsite){
      for (k in 1:nsurveys){
        y[i,j,k] ~ dcat(detP[u[i,j], ])
      }
    }
  }
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = bcalldata$ncell)
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha.lpsi", "alpha.lr", "beta.lpsi", "beta.lr", "psi", "r","theta2", "theta3", "p2", "p3", "Omega",
            "Theta", "detP", "z") # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
out_msms <- jags(bcalldata, inits, params, "ss_ms_ms.txt", n.adapt = na,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)


traceplot(out_msms)
print(out_msms, 3)


# Local probabilities
round(out_msms$mean$Theta, 2)

# Detection probabilities
round(out_msms$mean$detP, 2)





