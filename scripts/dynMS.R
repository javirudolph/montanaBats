#
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
  # create index for siteid.year, and cell.site
  mutate(site_id_year = paste(site_id, year, sep = "_"),
         cell_site = paste(cell, site, sep = "_"),
         cell_year = paste(cell, year, sep = "_")) %>%
  # create binary version of detection data
  mutate(det_non_det = as.numeric(y > 0)) %>%
  # Create state categories, y multi state
  # Asumming a cutoff of 100 bat calls as many bats
  mutate(y_multi_state = case_when(
    y == 0 ~ 1,
    y <= 100 ~ 2,
    y > 100 ~ 3
  )) %>%
  # add column for number of nights of recording for each site_id
  add_count(site_id, name = "n_nights_site_id") -> wns_calls_long


# The spatial scale we decided to work on was the 10x10km cell grid from NABat
# The data presented here has four sites sampled within a cell grid: multi-scale
# There is a varying number of nights surveyed for each site

# 3. Take medians -----------------------------

# I can't figure out how to incorporate the multi-scale component into the multi-state model
# My approach right now is to group the observations by cell (larger spatial scale)
# and then take the median for the number of calls on a given night if there is data for more than one site

# which cells have data for both years?

wns_calls_long %>%
  select(cell, year) %>%
  distinct() %>%
  count(cell) %>%
  filter(n == 2) %>%
  pull(cell) -> cell_2_years

# filter by the ones with two years of data
wns_calls_long %>%
  filter(cell %in% cell_2_years) %>%
  # group by julian date
  group_by(cell, year, jdate) %>%
  # trying the median number of calls
  summarise(y = median(y)) %>%
  # create the multi state variable using the many bats cutoff
  mutate(yms = case_when(
    y == 0 ~ 1,
    y <= 100 ~ 2,
    y > 100 ~ 3
  )) %>%
  mutate(cell.year = paste(cell, year, sep = "."),
         dnd = as.numeric(y > 0)) -> obs


# 4. Prep data to fit model ----------------------------------------

max.nsurveys <- max(as.numeric(table(obs$cell)))
sitelist <- sort(unique(obs$cell))

ncell <- length(unique(obs$cell))
nyears <- length(unique(obs$year))
yearm1 <- as.numeric(min(obs$year))

# Create empty arrays and fill them with the data
yms <- date <- y <- array(NA, dim = c(ncell, max.nsurveys, nyears),
                          dimnames = list(sitelist, 1:max.nsurveys, 1:nyears))

for(i in 1:ncell){
  for(t in 1:nyears){
    sel.cell.year <- paste(sitelist[i], t+yearm1, sep = ".")
    tmp <- obs[obs$cell.year == sel.cell.year,]
    nr <- nrow(tmp)
    if(nr > 0){
      yms[i, 1:nr,t] <- tmp$yms
      date[i, 1:nr, t] <- tmp$jdate
      y[i,1:nr,t] <- tmp$dnd
    }
  }
}


# 5. N. surveys -----------------------------------------

# There is variation in the number of nights that microphones are recording for each cell or site

nsurveys <- array(1, dim = c(ncell, nyears))
for(i in 1:ncell){
  for(t in 1:nyears){
    tmp <- which(!is.na(yms[i,,t]))
    if(length(tmp) > 0){
      nsurveys[i,t] <- max(tmp)
    }
  }
}


# 6. Fit null model ------------------------------------------
# No covariates

str(bcalldata <- list(y = yms, nsites = dim(yms)[1], nsurveys = nsurveys,
                      nyears = dim(yms)[3]))


# Dynamic multiseason model with the parameterization for the PhiMat using growth, decline, persistence, etc

# Specify model in BUGS language
cat(file = "jags_txt/nulldynMS.txt", "
model {
  ### (1) Priors for parameters
  # State process priors
  # Priors for parameters in initial state vector (Omega)
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)
  # Priors for parameters in state transition matrix (PhiMat)
  for(s in 1:2){
    phi[s] ~ dunif(0, 1)
  }
  lgamma ~ dunif(0,1)
  G ~ dunif(0,1)
  D ~ dunif(0,1)

  # Priors for parameters in observation process (Theta)
  p2 ~ dunif(0, 1)                 # Detection prob. when in state 2
  for (s in 1:3) {                 # Detection prob. when in state 3
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- beta[s] / sum(beta[])
  }
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector: Year 1
  Omega[1] <- 1 - psi              # Prob. of non-occupation
  Omega[2] <- psi * (1-r)          # Prob. of occ. by a few bats
  Omega[3] <- psi * r              # Prob. of occ. by many bats
  # Define state transition probability matrix (PhiMat): years 2:nyears
  # Define probabilities of state S(t+1) given S(t)
  # For now, constant over sites and years
  # Note conditional Bernoulli parameterization of multinomial
  # Order of indices: Departing state, arrival state
  PhiMat[1,1] <- 1 - lgamma
  PhiMat[1,2] <- lgamma
  PhiMat[1,3] <- 0
  PhiMat[2,1] <- 1 - phi[1]
  PhiMat[2,2] <- phi[1] * (1 - G)
  PhiMat[2,3] <- phi[1] * G
  PhiMat[3,1] <- 1 - phi[2]
  PhiMat[3,2] <- phi[2] * D
  PhiMat[3,3] <- phi[2] * (1 - D)
  # Define observation probability matrix (Theta)
  # Order of indices: true state, observed state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-p2
  Theta[2,2] <- p2
  Theta[2,3] <- 0
  Theta[3,1] <- p3[1]
  Theta[3,2] <- p3[2]
  Theta[3,3] <- p3[3]
  ### (3) Likelihood
  # Initial state: year 1
  for (i in 1:nsites){
    z[i,1] ~ dcat(Omega[])
  }
  # State transitions from yearly interval 1:(nyears-1)
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      z[i,t+1] ~ dcat(PhiMat[z[i,t],])
    }
  }
  # Observation equation
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        y[i,j,t] ~ dcat(Theta[z[i, t], ])
      }
    }
  }
  ### (4) Derived quantities
  # Number of sites in each state per year
  for (t in 1:nyears){
    for (i in 1:nsites){
      state1[i,t] <- equals(z[i,t], 1)   # Indicator for site in state 1
      state2[i,t] <- equals(z[i,t], 2)   # ... state 2
      state3[i,t] <- equals(z[i,t], 3)   # ... state 3
    }
    n.occ[t,1] <- sum(state1[,t])        # Number of unoccupied sites
    n.occ[t,2] <- sum(state2[,t])        # Number of sites with few bats
    n.occ[t,3] <- sum(state3[,t])        # Number of sites with many bats
    n.occ.total[t] <- n.occ[t,2] + n.occ[t, 3] # All occupied
  }
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(bcalldata$nsites, bcalldata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("psi", "r", "phi", "lgamma", "G", "D", "p2", "p3", "Omega", "PhiMat",
            "Theta", "n.occ", "n.occ.total") # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
odms_null <- jags(bcalldata, inits, params, "jags_txt/nulldynMS.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

saveRDS(odms_null, file = "modelouts/odms_null.RDS")

traceplot(odms_null)
print(odms_null, 3)



# PhiMat
round(odms_null$mean$PhiMat, 2)

# Observation probabilities
round(odms_null$mean$Theta, 2)



# 7. Fit with covariates -----------------------------------


## 7.a. Site covariates -------------------------------------
raw_sites <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/SiteDataCleanFinal.csv", show_col_types = FALSE)

# I honestly don't know what we can do with these covariates. They are also at the site level (smaller scale)
# Maybe we could use the site type, but I don't even know what these site types mean

## 7.b. NABat covariates ---------------------------------------
# source: https://www.sciencebase.gov/catalog/item/620e6f3bd34e6c7e83baa603

mt_covariates <- read_sf("datafiles/nabat_covariates/NABat_grid_covariates/NABat_grid_covariates.shp") %>%
  filter(., admin1 == "Montana")

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

mt_covariates %>% 
  filter(GRTS_ID %in% cell_2_years) %>% 
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

#bundle data
str(bcalldata <- list(y = yms, nsites = dim(yms)[1], nsurveys = nsurveys,
                      nyears = dim(yms)[3], elev = elev.scaled, forest = forest.scaled))

# not sure how you include an intercept here for the logit link

cat(file = "jags_txt/simpledynMS.txt", "
model {
  
  ### Linear models and priors 
  for (i in 1:nsites){
    logit(psi[i]) <- beta.lpsi[1] * forest[i] + beta.lpsi[2] * elev[i]
    logit(r[i]) <- beta.lr[1] * forest[i] + beta.lr[2] * elev[i]
  }
  
  # Coefficients fo the 2 covariates in Omega
  for (k in 1:2){
    beta.lpsi[k] ~ dnorm(0, 0.1)
    beta.lr[k] ~ dnorm(0, 0.1)
  }
  
  ### (1) Priors for parameters
  # State process priors
  # Priors for parameters in initial state vector (Omega)
  # psi ~ dunif(0, 1)
  # r ~ dunif(0, 1)
  
  # Priors for parameters in state transition matrix (PhiMat)
  for(s in 1:2){
    phi[s] ~ dunif(0, 1)
  }
  lgamma ~ dunif(0,1)
  G ~ dunif(0,1)
  D ~ dunif(0,1)

  # Priors for parameters in observation process (Theta)
  p2 ~ dunif(0, 1)                 # Detection prob. when in state 2
  for (s in 1:3) {                 # Detection prob. when in state 3
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- beta[s] / sum(beta[])
  }
  
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector: Year 1
  for (i in 1:nsites){
    Omega[i,1] <- 1 - psi[i]                 # Prob. of non-occupation
    Omega[i,2] <- psi[i] * (1-r[i])          # Prob. of occ. by a few bats
    Omega[i,3] <- psi[i] * r[i]              # Prob. of occ. by many bats
  }
  
  
  # Define state transition probability matrix (PhiMat): years 2:nyears
  # Define probabilities of state S(t+1) given S(t)
  # For now, constant over sites and years
  # Note conditional Bernoulli parameterization of multinomial
  # Order of indices: Departing state, arrival state
  PhiMat[1,1] <- 1 - lgamma
  PhiMat[1,2] <- lgamma
  PhiMat[1,3] <- 0
  PhiMat[2,1] <- 1 - phi[1]
  PhiMat[2,2] <- phi[1] * (1 - G)
  PhiMat[2,3] <- phi[1] * G
  PhiMat[3,1] <- 1 - phi[2]
  PhiMat[3,2] <- phi[2] * D
  PhiMat[3,3] <- phi[2] * (1 - D)
  # Define observation probability matrix (Theta)
  # Order of indices: true state, observed state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-p2
  Theta[2,2] <- p2
  Theta[2,3] <- 0
  Theta[3,1] <- p3[1]
  Theta[3,2] <- p3[2]
  Theta[3,3] <- p3[3]
  
  ### (3) Likelihood
  # Initial state: year 1
  for (i in 1:nsites){
    z[i,1] ~ dcat(Omega[i,])
  }
  # State transitions from yearly interval 1:(nyears-1)
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      z[i,t+1] ~ dcat(PhiMat[z[i,t],])
    }
  }
  # Observation equation
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        y[i,j,t] ~ dcat(Theta[z[i, t], ])
      }
    }
  }
  ### (4) Derived quantities
  # Number of sites in each state per year
  for (t in 1:nyears){
    for (i in 1:nsites){
      state1[i,t] <- equals(z[i,t], 1)   # Indicator for site in state 1
      state2[i,t] <- equals(z[i,t], 2)   # ... state 2
      state3[i,t] <- equals(z[i,t], 3)   # ... state 3
    }
    n.occ[t,1] <- sum(state1[,t])        # Number of unoccupied sites
    n.occ[t,2] <- sum(state2[,t])        # Number of sites with few bats
    n.occ[t,3] <- sum(state3[,t])        # Number of sites with many bats
    n.occ.total[t] <- n.occ[t,2] + n.occ[t, 3] # All occupied
  }
}
")


# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(bcalldata$nsites, bcalldata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("beta.lpsi", "beta.lr", "psi", "r", "phi", "lgamma", "G", "D", "p2", "p3", "Omega", "PhiMat",
            "Theta", "n.occ", "n.occ.total") # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
odms <- jags(bcalldata, inits, params, "jags_txt/simpledynMS.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

traceplot(odms)
print(odms, 3)

# PhiMat
round(odms$mean$PhiMat, 2)

# Observation probabilities
round(odms$mean$Theta, 2)


# 8. More covariates ---------------------------------------------------------

mt_covariates %>% 
  filter(GRTS_ID %in% cell_2_years) %>% 
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



#bundle data
str(bcalldata <- list(y = yms, nsites = dim(yms)[1], nsurveys = nsurveys, nyears = dim(yms)[3], 
                      region = region,
                      elev = elev.scaled, 
                      temp = temp.scaled,
                      physdiv = physdiv.scaled,
                      precip = precip.scaled,
                      forest = forest.scaled,
                      wetlands = wetlands.scaled,
                      karst = karst))

# not sure how you include an intercept here for the logit link

cat(file = "jags_txt/dynMS_wcovs.txt", "
model {
  
  ### Linear models and priors 
  for (i in 1:nsites){
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
  
  # Coefficients of covariates in Omega
  for (k in 1:8){
    beta.lpsi[k] ~ dnorm(0, 0.1)
  }
  
  for (k in 1:3){
    beta.lr[k] ~ dnorm(0, 0.1)
  }
  
  # Priors for parameters in state transition matrix (PhiMat)
  for(s in 1:2){
    phi[s] ~ dunif(0, 1)
  }
  lgamma ~ dunif(0,1)
  G ~ dunif(0,1)
  D ~ dunif(0,1)

  # Priors for parameters in observation process (Theta)
  p2 ~ dunif(0, 1)                 # Detection prob. when in state 2
  for (s in 1:3) {                 # Detection prob. when in state 3
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- beta[s] / sum(beta[])
  }
  
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector: Year 1
  for (i in 1:nsites){
    Omega[i,1] <- 1 - psi[i]                 # Prob. of non-occupation
    Omega[i,2] <- psi[i] * (1-r[i])          # Prob. of occ. by a few bats
    Omega[i,3] <- psi[i] * r[i]              # Prob. of occ. by many bats
  }
  
  
  # Define state transition probability matrix (PhiMat): years 2:nyears
  # Define probabilities of state S(t+1) given S(t)
  # For now, constant over sites and years
  # Note conditional Bernoulli parameterization of multinomial
  # Order of indices: Departing state, arrival state
  PhiMat[1,1] <- 1 - lgamma
  PhiMat[1,2] <- lgamma
  PhiMat[1,3] <- 0
  PhiMat[2,1] <- 1 - phi[1]
  PhiMat[2,2] <- phi[1] * (1 - G)
  PhiMat[2,3] <- phi[1] * G
  PhiMat[3,1] <- 1 - phi[2]
  PhiMat[3,2] <- phi[2] * D
  PhiMat[3,3] <- phi[2] * (1 - D)
  # Define observation probability matrix (Theta)
  # Order of indices: true state, observed state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-p2
  Theta[2,2] <- p2
  Theta[2,3] <- 0
  Theta[3,1] <- p3[1]
  Theta[3,2] <- p3[2]
  Theta[3,3] <- p3[3]
  
  ### (3) Likelihood
  # Initial state: year 1
  for (i in 1:nsites){
    z[i,1] ~ dcat(Omega[i,])
  }
  # State transitions from yearly interval 1:(nyears-1)
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      z[i,t+1] ~ dcat(PhiMat[z[i,t],])
    }
  }
  # Observation equation
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        y[i,j,t] ~ dcat(Theta[z[i, t], ])
      }
    }
  }
  ### (4) Derived quantities
  # Number of sites in each state per year
  for (t in 1:nyears){
    for (i in 1:nsites){
      state1[i,t] <- equals(z[i,t], 1)   # Indicator for site in state 1
      state2[i,t] <- equals(z[i,t], 2)   # ... state 2
      state3[i,t] <- equals(z[i,t], 3)   # ... state 3
    }
    n.occ[t,1] <- sum(state1[,t])        # Number of unoccupied sites
    n.occ[t,2] <- sum(state2[,t])        # Number of sites with few bats
    n.occ[t,3] <- sum(state3[,t])        # Number of sites with many bats
    n.occ.total[t] <- n.occ[t,2] + n.occ[t, 3] # All occupied
  }
}
")


# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(bcalldata$nsites, bcalldata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha.lpsi", "alpha.lr", "beta.lpsi", "beta.lr", "psi", "r", "phi", "lgamma", "G", "D", "p2", "p3", "Omega", "PhiMat",
            "Theta", "n.occ", "n.occ.total") # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
odms_wcovs <- jags(bcalldata, inits, params, "jags_txt/dynMS_wcovs.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

saveRDS(odms_wcovs, file = "modelouts/odms_wcovs.RDS")

#traceplot(odms_wcovs)
print(odms_wcovs, 3)

# PhiMat
round(odms_wcovs$mean$PhiMat, 2)

# Observation probabilities
round(odms_wcovs$mean$Theta, 2)

# 9. MultiScale model --------------------------------------------------

# We have 4 sites in each cell
# Except cell 1313 which shows as having sites 1,2,3,4 in year 2021, and sites 1,3,4,5 in year 2020. 
# These sites are not always the same location, so for now, changing that 5 to a 2

multiscale_yms <- wns_calls_long
multiscale_yms$site[multiscale_yms$site == "5"] <- "2"


# Sit on this....
# I think I need a cell and siteId variable. 
# each site is 'unique' and doesn't always get resampled...
# and the variables there that are important change from one year to the next...
#

# set up data
multiscale_yms %>% 
  mutate(cellsiteid = paste(cell, site_id, sep = "_"),
         cellsiteidyear = paste(cellsiteid, year, sep = "_")) %>% 
  # create the multi state variable using the many bats cutoff
  mutate(yms = y_multi_state,
         dnd = det_non_det) -> obs2


max.nsurveys <- max(as.numeric(table(obs2$cellsiteid)))
sitelist <- sort(unique(obs2$cellsiteid))
celllist <- sort(unique(obs2$cell))

ncell <- length(unique(obs2$cell))
nsites <- length(unique(obs2$cellsiteid))
nyears <- length(unique(obs2$year))
yearm1 <- as.numeric(min(obs2$year))-1

# Create empty arrays and fill them with the data
yms <- date <- y <- array(NA, dim = c(nsites, max.nsurveys, nyears),
                          dimnames = list(sitelist, 1:max.nsurveys, 1:nyears))

for(i in 1:nsites){
  for(t in 1:nyears){
    t<-1
    sel.cell.year <- paste(sitelist[i], t+yearm1, sep = "_")
    tmp <- obs2[obs2$cellsiteidyear == sel.cell.year,]
    nr <- nrow(tmp)
    if(nr > 0){
      yms[i, 1:nr,t] <- tmp$yms
      date[i, 1:nr, t] <- tmp$jdate
      y[i,1:nr,t] <- tmp$dnd
    }
  }
}


## N. surveys -----------------------------------------

# There is variation in the number of nights that microphones are recording for each cell or site

nsurveys <- array(1, dim = c(nsites, nyears))
for(i in 1:nsites){
  for(t in 1:nyears){
    tmp <- which(!is.na(yms[i,,t]))
    if(length(tmp) > 0){
      nsurveys[i,t] <- max(tmp)
    }
  }
}


##  Fit null model ------------------------------------------
# No covariates

str(bcalldata <- list(y = yms, nsites = dim(yms)[1], nsurveys = nsurveys,
                      nyears = dim(yms)[3]))


# Dynamic multiseason model with the parameterization for the PhiMat using growth, decline, persistence, etc

# Specify model in BUGS language
cat(file = "jags_txt/nulldynMSMS.txt", "
model {
  ### (1) Priors for parameters
  # State process priors
  # Priors for parameters in initial state vector (Omega)
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)
  # Priors for parameters in state transition matrix (PhiMat)
  for(s in 1:2){
    phi[s] ~ dunif(0, 1)
  }
  lgamma ~ dunif(0,1)
  G ~ dunif(0,1)
  D ~ dunif(0,1)
  
  # Priors for parameters in local occupancy (Theta)
  theta2 ~ dunif(0, 1)              # Detection prob. when in state 2
  for (s in 1:3) {                 # Detection prob. when in state 3
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    theta3[s] <- beta[s] / sum(beta[])
  }

  # Priors for parameters in observation process (varp)
  p2 ~ dunif(0, 1)                 # Detection prob. when in state 2
  for (s in 1:3) {                 # Detection prob. when in state 3
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- beta[s] / sum(beta[])
  }
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector: Year 1
  Omega[1] <- 1 - psi              # Prob. of non-occupation
  Omega[2] <- psi * (1-r)          # Prob. of occ. by a few bats
  Omega[3] <- psi * r              # Prob. of occ. by many bats
  # Define state transition probability matrix (PhiMat): years 2:nyears
  # Define probabilities of state S(t+1) given S(t)
  # For now, constant over sites and years
  # Note conditional Bernoulli parameterization of multinomial
  # Order of indices: Departing state, arrival state
  PhiMat[1,1] <- 1 - lgamma
  PhiMat[1,2] <- lgamma
  PhiMat[1,3] <- 0
  PhiMat[2,1] <- 1 - phi[1]
  PhiMat[2,2] <- phi[1] * (1 - G)
  PhiMat[2,3] <- phi[1] * G
  PhiMat[3,1] <- 1 - phi[2]
  PhiMat[3,2] <- phi[2] * D
  PhiMat[3,3] <- phi[2] * (1 - D)
  
  # Define local occupancy probability (Theta)
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta2
  Theta[2,2] <- theta2
  Theta[2,3] <- 0
  Theta[3,1] <- theta3[1]
  Theta[3,2] <- theta3[2]
  Theta[3,3] <- theta3[3]
  # Define observation probability matrix (p)
  # Order of indices: true state, observed state
  varp[1,1] <- 1
  varp[1,2] <- 0
  varp[1,3] <- 0
  varp[2,1] <- 1-p2
  varp[2,2] <- p2
  varp[2,3] <- 0
  varp[3,1] <- p3[1]
  varp[3,2] <- p3[2]
  varp[3,3] <- p3[3]
  
  ### (3) Likelihood
  # Initial state: year 1
  for (i in 1:nsites){
    z[i,1] ~ dcat(Omega[])
  }
  # State transitions from yearly interval 1:(nyears-1)
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      z[i,t+1] ~ dcat(PhiMat[z[i,t],])
    }
  }
  
  # local occupancy
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        u[i,j,t] ~ dcat(Theta[z[i, t], ])
      }
    }
  }
  
  
  # Observation equation
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        for(k in )
      }
    }
  }
  ### (4) Derived quantities
  # Number of sites in each state per year
  for (t in 1:nyears){
    for (i in 1:nsites){
      state1[i,t] <- equals(z[i,t], 1)   # Indicator for site in state 1
      state2[i,t] <- equals(z[i,t], 2)   # ... state 2
      state3[i,t] <- equals(z[i,t], 3)   # ... state 3
    }
    n.occ[t,1] <- sum(state1[,t])        # Number of unoccupied sites
    n.occ[t,2] <- sum(state2[,t])        # Number of sites with few bats
    n.occ[t,3] <- sum(state3[,t])        # Number of sites with many bats
    n.occ.total[t] <- n.occ[t,2] + n.occ[t, 3] # All occupied
  }
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(bcalldata$nsites, bcalldata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("psi", "r", "phi", "lgamma", "G", "D", "p2", "p3", "Omega", "PhiMat",
            "Theta", "n.occ", "n.occ.total") # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
odms_null <- jags(bcalldata, inits, params, "jags_txt/nulldynMS.txt", n.adapt = na,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

saveRDS(odms_null, file = "modelouts/odms_null.RDS")

#traceplot(odms_null)
print(odms_null, 3)



# PhiMat
round(odms_null$mean$PhiMat, 2)

# Observation probabilities
round(odms_null$mean$Theta, 2)



