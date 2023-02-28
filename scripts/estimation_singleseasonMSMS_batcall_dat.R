#   montanaBats Project
#   Management of White Nose Syndrome


#   Code written by F. Javiera Rudolph, Ph.D.

library(tidyverse)
library(cowplot)
library(sf)
library(AHMbook)

library(jagsUI)
library(mcmcOutput)

# 1. Data sources -----------------------------------

## 1a. Bat call Data -------------------------------
# Data here is not specific to any bat species, but it is a subset of calls
# the frequency used to subset corresponds to bat species known to be 
# susceptible to white nose syndrome

raw_mtcalls <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/WNSCallsCapHistory.csv")

raw_mtcalls %>% 
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
  separate(csy, c("cell", "site", "year")) %>%
  # create binary version of detection data
  mutate(dnd = as.numeric(y > 0)) -> mt_calls

# Summary of the bat call data

# bat_threshold <- median(mt_calls$y)
bat_threshold <- 100
ggplot(data = mt_calls, aes(x = y)) + 
  geom_histogram() + 
  geom_vline(xintercept = bat_threshold, color = "red") -> call_dist_plot

mt_calls %>% 
  mutate(yms = case_when(
    y == 0 ~ 1,
    y <= bat_threshold ~ 2, 
    y > bat_threshold ~ 3
  )) -> mt_calls

ggplot(data = mt_calls, aes(x = yms, fill = year)) +
  geom_bar() -> state_dist_plot

plot_grid(call_dist_plot, state_dist_plot)

# naive occupancy probability not accounting for imperfect detection or multilevel occupancy
# nested sites

sum(mt_calls$dnd)/nrow(mt_calls)
table(mt_calls$yms)
table(mt_calls$yms)/nrow(mt_calls)


# Number of cells, sites per cell, surveys per site:
cell_list <- sort(unique(mt_calls$cell))
ncells <- length(cell_list)

# distribution of number of surveys per site
mt_calls %>% 
  count(cell, site) %>% 
  ggplot(., aes(x = n)) +
  geom_histogram()

# How many sites per cell per year
mt_calls %>% 
  select(site_id, cell, year) %>% 
  distinct() %>% 
  count(cell, year) %>% 
  count(n)

max_nsites <- mt_calls %>% select(site_id, cell) %>% distinct() %>% count(cell) %>% pull(n) %>% max()

# how many surveys/nights per site:
mt_calls %>% 
  count(site_id, cell) %>% 
  pull(n) %>% 
  max() -> max_surveys

# How many cells have info for both years?
mt_calls %>% 
  select(site_id, cell, year) %>% 
  distinct() %>% 
  select(cell, year) %>% 
  distinct() %>% 
  count(cell) %>% 
  filter(n ==2) %>% 
  pull(cell) -> cell_list_2years

length(cell_list_2years)


## 1b. Site covariates -----------------

raw_sitecovs <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/SiteDataCleanFinal.csv") %>% 
                    mutate(Site_Type = str_remove(Site_Type, " "))

sitecovs <- raw_sitecovs %>% 
  janitor::clean_names() %>% 
  select(cell, site_id, year, latitude, longitude, 
         jstart, jstop,
         site_type, other_site_type) %>% 
  mutate(duration = jstop - jstart,
         duration = standardize(duration))

# interested in site type as a covariate for local availability
unique(sitecovs$site_type)

# What falls under the "other" category
sitecovs %>% 
  filter(site_type == "Other") %>% View()

sitecovs %>% 
  # filter(site_type == "Other") %>%
  mutate(other_site_type = str_to_lower(other_site_type)) %>% 
  mutate(fix_other = case_when(
    str_detect(other_site_type, "meadow|ag") ~ "Meadow",
    str_detect(other_site_type, "talus|slope") ~ "RockOut",
    str_detect(other_site_type, "lotic|river") ~ "LoticWater",
    str_detect(other_site_type, "wooded") ~ "OtherRoost"
  )) %>% 
  mutate(site_type = ifelse(site_type == "Other", fix_other, site_type),
         site_type_num = as.numeric(factor(site_type))) %>% 
  select(-fix_other) -> sitecovs

table(sitecovs$site_type)


# scale the starting survey dates
hist(sitecovs$jstart)
summary(sitecovs$jstart)
mean(sitecovs$jstart) # July 1st is 182
sd(sitecovs$jstart)

hist(standardize(sitecovs$jstart))
hist((sitecovs$jstart-182)/15) # scale so 0 is July 1st, and sd = 2 weeks

sitecovs %>% 
  mutate(st_jstart = (jstart-182)/15) -> sitecovs


## 1c. NABat covariates ------------------------
mt_covariates <- read_sf("datafiles/nabat_covariates/NABat_grid_covariates/NABat_grid_covariates.shp") %>%
  filter(., admin1 == "Montana")

mt_covariates %>% 
  filter(GRTS_ID %in% cell_list) %>% 
  select(GRTS_ID, karst, p_forest, p_wetland, mean_temp,
         precip, DEM_max, physio_div, dist_mines, starts_with("eco")) %>% 
  # arrange to make sure its same order of cells as obs data
  arrange(factor(GRTS_ID, levels = cell_list)) %>% 
  rename(cell = GRTS_ID) %>% 
  mutate(region = as.numeric(factor(eco3_name))) -> cell_covs

# standardize variables
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
region <- cell_covs$region


# Models ----------------
## Null ----------------

# Multi-state multi-scale occupancy model
# no covariates

# Cell level occupancy
# z_i ~ multinom(Omega)
# Omega = [1-psi, psi(1-r), psi*r]

# Site level availability
# u_{i,j} ~ multinom(Theta)
# Theta = [1, 0, 0
#          1-theta2, theta2, 0
#          1-theta32-theta33, theta32, theta33]

# Observation process
# y_{i,j,k} ~ multinom(p)
# p = [1,0,0
#      1-p2, p2, 0,
#      1-p32-p33, p32, p33]

# Create new site identifier that pools data from both years
# and then creates an identifier for cell_site
mt_calls %>% 
  group_split(cell) %>% 
  purrr::map_df(~.x %>% group_by(site_id) %>% mutate(id = cur_group_id())) %>% 
  mutate(cell_site = paste(cell, id, sep = "_")) -> null_msms_dat


# bundle data
yms <- array(NA, dim = c(ncells, max_nsites, max_surveys),
             dimnames = list(cell_list, 1:max_nsites, 1:max_surveys))

for(i in 1:ncells){
  for(j in 1:max_nsites){
    sel_cell_site <- paste(cell_list[i], j, sep = "_")
    tmp <- null_msms_dat[null_msms_dat$cell_site == sel_cell_site,]
    nr <- nrow(tmp)
    if(nr > 0){
      yms[i,j, 1:nr] <- tmp$yms
    }
  }
}

# There is variation in the number of nights that microphones are recording for each cell or site
# to make it faster, have the model estimate only for when there is data, and not run NA cycles
nsurveys <- array(NA, dim = c(ncells, max_nsites))
for(i in 1:ncells){
  for(j in 1:max_nsites){
    tmp <- which(!is.na(yms[i,j,]))
    if(length(tmp) > 0){
      nsurveys[i,j] <- max(tmp)
    }
  }
}

nsites <- array(NA, dim = ncells)
for(i in 1:ncells){
  tmp <- which(!is.na(yms[i,,1]))
  if(length(tmp) > 0){
    nsites[i] <- max(tmp)
  }
}


# same for number of sites. Since not all cells were sampled in both years, so some have only 4 sites instead of 8


str(batdata <- list(y = yms, ncells = dim(yms)[1], nsites = dim(yms)[2], nsurveys = dim(yms)[3]))

# Specify the model

cat(file = "jags_txt/null_msms.txt", "
model {
  # Priors
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)
  theta2 ~ dunif(0,1)
  p2 ~ dunif(0, 1)
  
  # Multinomial logit link for availability model for state 3 (= many bats)
  ltheta32 ~ dnorm(0, 0.001)
  ltheta33 ~ dnorm(0, 0.001)
  theta32 <- exp(ltheta32) / (1 + exp(ltheta32) + exp(ltheta33))
  theta33 <- exp(ltheta33) / (1 + exp(ltheta32) + exp(ltheta33))
  theta31 <- 1-theta32-theta33  
  
  # Multinomial logit link for observation model for state 3
  lp32 ~ dnorm(0, 0.001)
  lp33 ~ dnorm(0, 0.001)
  p32 <- exp(lp32) / (1 + exp(lp32) + exp(lp33))
  p33 <- exp(lp33) / (1 + exp(lp32) + exp(lp33))
  p31 <- 1-p32-p33                     
  
  # Define initial state vector (Omega)
  Omega[1] <- 1 - psi                  # Prob. of no bats
  Omega[2] <- psi * (1-r)              # Prob. of occupancy (w/ few bats)
  Omega[3] <- psi * r                  # Prob. of occupancy (with many bats)
  
  # Define availability matrix (Theta)
  # Order of indices: true state, observed state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta2
  Theta[2,2] <- theta2
  Theta[2,3] <- 0
  Theta[3,1] <- theta31                    
  Theta[3,2] <- theta32
  Theta[3,3] <- theta33
  
  # Define observation matrix (pDet)
  # Order of indices: true state, observed state
  pDet[1,1] <- 1
  pDet[1,2] <- 0.000001
  pDet[1,3] <- 0.000001
  pDet[2,1] <- 1-p2
  pDet[2,2] <- p2
  pDet[2,3] <- 0.000001
  pDet[3,1] <- p31                    # = 1-p32-p33 as per prior section
  pDet[3,2] <- p32
  pDet[3,3] <- p33
  
  # State-space likelihood
  # State equation: model of true states (z)
  for (i in 1:ncells){
    z[i] ~ dcat(Omega[])
  }
  
  # Availability equation
  for (i in 1:ncells){
    for (j in 1:nsites){
      a[i,j] ~ dcat(Theta[z[i],])
    }
  }
  
  # Observation equation
  for (i in 1:ncells){
    for (j in 1:nsites){
      for (k in 1:nsurveys){
        y[i,j,k] ~ dcat(pDet[a[i,j],])
      }
    }
  }
  
  # Derived quantities
  for (i in 1:ncells){
    occ1[i] <- equals(z[i], 1)
    occ2[i] <- equals(z[i], 2)
    occ3[i] <- equals(z[i], 3)
  }
  n.occ[1] <- sum(occ1[]) # Grid cells in state 1
  n.occ[2] <- sum(occ2[]) # Grid cells in state 2
  n.occ[3] <- sum(occ3[]) # Grid cells in state 3
}
")

# Initial values
zst <- rep(3, nrow(batdata$y)) # Initialize at highest possible state
inits <- function(){list(z = zst)}

# Parameters monitored (could add "z")
params <- c("psi", "r", "p2", "p31", "p32", "p33", 
            "theta2", "theta31", "theta32", "theta33",
            "Omega", "Theta", "pDet", "n.occ")

# MCMC settings
na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 2000 ; nc <- 3

# Call JAGS, check convergence and summarize posteriors
out_null_msms <- jags(batdata, inits, params, "jags_txt/null_msms.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
# traceplot(out_null_msms)
print(out_null_msms, 3)

diagPlot(out_null_msms)


## M1 - cellcovs --------------------
# Incorporate only covariates on psi and r

str(batdata <- list(y = yms, ncells = dim(yms)[1], nsites = dim(yms)[2], nsurveys = max_surveys, 
                      region = region,
                      elev = elev.scaled, 
                      temp = temp.scaled,
                      physdiv = physdiv.scaled,
                      precip = precip.scaled,
                      forest = forest.scaled,
                      wetlands = wetlands.scaled,
                      karst = karst))


cat(file = "jags_txt/M1_msms.txt", "
model {


  # Priors
  
  # Omega
  for (i in 1:ncells){
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
  
  # Multinomial logit link for availability model for state 3 (= many bats)
  theta2 ~ dunif(0,1)
  ltheta32 ~ dnorm(0, 0.001)
  ltheta33 ~ dnorm(0, 0.001)
  theta32 <- exp(ltheta32) / (1 + exp(ltheta32) + exp(ltheta33))
  theta33 <- exp(ltheta33) / (1 + exp(ltheta32) + exp(ltheta33))
  theta31 <- 1-theta32-theta33  
  
  # Multinomial logit link for observation model for state 3
  p2 ~ dunif(0, 1)
  lp32 ~ dnorm(0, 0.001)
  lp33 ~ dnorm(0, 0.001)
  p32 <- exp(lp32) / (1 + exp(lp32) + exp(lp33))
  p33 <- exp(lp33) / (1 + exp(lp32) + exp(lp33))
  p31 <- 1-p32-p33                     
  
  # Define initial state vector (Omega)
  for (i in 1:ncells){
    Omega[i,1] <- 1 - psi[i]                 # Prob. of no bats
    Omega[i,2] <- psi[i] * (1-r[i])          # Prob. of occ. by a few bats
    Omega[i,3] <- psi[i] * r[i]              # Prob. of occ. by many bats
  }
  
  # Define availability matrix (Theta)
  # Order of indices: true state, observed state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta2
  Theta[2,2] <- theta2
  Theta[2,3] <- 0
  Theta[3,1] <- theta31                    
  Theta[3,2] <- theta32
  Theta[3,3] <- theta33
  
  # Define observation matrix (pDet)
  # Order of indices: true state, observed state
  pDet[1,1] <- 1
  pDet[1,2] <- 0.000001
  pDet[1,3] <- 0.000001
  pDet[2,1] <- 1-p2
  pDet[2,2] <- p2
  pDet[2,3] <- 0
  pDet[3,1] <- p31                    # = 1-p32-p33 as per prior section
  pDet[3,2] <- p32
  pDet[3,3] <- p33
  
  # State-space likelihood
  # State equation: model of true states (z)
  for (i in 1:ncells){
    z[i] ~ dcat(Omega[i,])
  }
  
  # Availability equation
  for (i in 1:ncells){
    for (j in 1:nsites){
      a[i,j] ~ dcat(Theta[z[i],])
    }
  }
  
  # Observation equation
  for (i in 1:ncells){
    for (j in 1:nsites){
      for (k in 1:nsurveys){
        y[i,j,k] ~ dcat(pDet[a[i,j],])
      }
    }
  }
  
  # Derived quantities
  for (i in 1:ncells){
    occ1[i] <- equals(z[i], 1)
    occ2[i] <- equals(z[i], 2)
    occ3[i] <- equals(z[i], 3)
  }
  n.occ[1] <- sum(occ1[]) # Grid cells in state 1
  n.occ[2] <- sum(occ2[]) # Grid cells in state 2
  n.occ[3] <- sum(occ3[]) # Grid cells in state 3
}
")

# Initial values
zst <- rep(3, nrow(batdata$y)) # Initialize at highest possible state
inits <- function(){list(z = zst)}

# Parameters monitored (could add "z")
params <- c("alpha.lpsi", "alpha.lr", "beta.lpsi", "beta.lr",
            # "psi", "r", 
            "p2", "p31", "p32", "p33",
            "theta2", "theta31", "theta32", "theta33",
            # "Omega", "Theta", "pDet",
            "n.occ")

# MCMC settings
na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 2000 ; nc <- 3

# Call JAGS, check convergence and summarize posteriors
out_M1_msms <- jags(batdata, inits, params, "jags_txt/M1_msms.txt", n.adapt = na,
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
# traceplot(out_M1_msms)
print(out_M1_msms, 3)

diagPlot(out_M1_msms)

## M2 - sitecovs ---------------------------
# add julian date start, duration, and site type as availability covariates

jstart <- jstart_sqrd <- duration <- site_type <- array(NA, dim = c(ncells, max_nsites),
                                         dimnames = list(cell_list, 1:max_nsites))

# duration <- nsurveys


# match the cell_site with the one for the call data

null_msms_dat %>% 
  select(site_id, cell_site) %>% 
  left_join(., sitecovs) -> M2_sitecovs

for(i in 1:ncells){
  for(j in 1:max_nsites){
    sel_cell_site <- paste(cell_list[i], j, sep = "_")
    tmp <- M2_sitecovs[M2_sitecovs$cell_site == sel_cell_site,]
    nr <- nrow(tmp)
    if(nr > 0){
      jstart[i,j] <- tmp$st_jstart[1]
      jstart_sqrd[i,j] <- (tmp$st_jstart[1])^2 
      duration[i,j] <- tmp$duration[1]
      site_type[i,j] <- tmp$site_type_num[1]
    }
  }
}

nsites


str(batdata <- list(y = yms, ncells = dim(yms)[1], nsites = nsites, nsurveys = max_surveys, 
                    region = region,
                    elev = elev.scaled, 
                    temp = temp.scaled,
                    physdiv = physdiv.scaled,
                    precip = precip.scaled,
                    forest = forest.scaled,
                    wetlands = wetlands.scaled,
                    karst = karst,
                    date = jstart,
                    date_sqrd = jstart_sqrd,
                    site_type = site_type,
                    duration = duration))


cat(file = "jags_txt/M2_msms.txt", "
model {
  
  # Omega
  for (i in 1:ncells){
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
  
  # availability model with covariates
  for (i in 1:ncells){
    for (j in 1:nsites[i]){
      logit(theta2[i,j]) <- alpha.ltheta2 + beta.site.type.ltheta2[site_type[i,j]] +
        beta.ltheta2[1] * date[i,j] + beta.ltheta2[2] * date_sqrd[i,j] + 
        beta.ltheta2[3] * duration[i,j]
      
      mlogit.theta3[2,i,j] <- alpha.ltheta32 + beta.site.type.ltheta32[site_type[i,j]] +
        beta.ltheta32[1] * date[i,j] + beta.ltheta32[2] * date_sqrd[i,j] + 
        beta.ltheta32[3] * duration[i,j]
        
      mlogit.theta3[3,i,j] <- alpha.ltheta33 + beta.site.type.ltheta33[site_type[i,j]] +
        beta.ltheta33[1] * date[i,j] + beta.ltheta33[2] * date_sqrd[i,j] + 
        beta.ltheta33[3] * duration[i,j]
    }
  }
  
  # Priors for parameters in the linear models in Theta
  # intercepts --------
  alpha.ltheta2 <- logit(mean.alpha.theta2)
  mean.alpha.theta2 ~ dunif(0,1)
  alpha.ltheta32 ~ dnorm(0, 0.01)
  alpha.ltheta33 ~ dnorm(0, 0.01)
  
  # effects of site type -------
  beta.site.type.ltheta2[1] <- 0
  beta.site.type.ltheta32[1] <- 0
  beta.site.type.ltheta33[1] <- 0
  for (k in 2:8){
    beta.site.type.ltheta2[k] ~ dnorm(0, 0.1)
    beta.site.type.ltheta32[k] ~ dnorm(0, 0.1)
    beta.site.type.ltheta33[k] ~ dnorm(0, 0.1)
  }
  # coefficients of date and duration
  for (k in 1:3){
    beta.ltheta2[k] ~ dnorm(0, 0.1)
    beta.ltheta32[k] ~ dnorm(0, 0.1)
    beta.ltheta33[k] ~ dnorm(0, 0.1)
  }
  
  
  # Multinomial logit link for availability model for state 3 (= many bats)
  for (i in 1:ncells){
    for (j in 1:nsites[i]){
      theta3[2, i, j] <- exp(mlogit.theta3[2,i,j]) / (1 + exp(mlogit.theta3[2,i,j]) + exp(mlogit.theta3[3,i,j]))
      theta3[3, i, j] <- exp(mlogit.theta3[3,i,j]) / (1 + exp(mlogit.theta3[2,i,j]) + exp(mlogit.theta3[3,i,j]))
      theta3[1, i, j] <- 1 - theta3[2, i, j] - theta3[3, i, j]
    }
  }
  
  # Multinomial logit link for observation model for state 3
  p2 ~ dunif(0, 1)
  lp32 ~ dnorm(0, 0.001)
  lp33 ~ dnorm(0, 0.001)
  p32 <- exp(lp32) / (1 + exp(lp32) + exp(lp33))
  p33 <- exp(lp33) / (1 + exp(lp32) + exp(lp33))
  p31 <- 1-p32-p33                     
  
  # Define initial state vector (Omega)
  for (i in 1:ncells){
    Omega[i,1] <- 1 - psi[i]                 # Prob. of no bats
    Omega[i,2] <- psi[i] * (1-r[i])          # Prob. of occ. by a few bats
    Omega[i,3] <- psi[i] * r[i]              # Prob. of occ. by many bats
  }
  
  # Define availability matrix (Theta)
  # Order of indices: true state, observed state for each cell_site 
  for (i in 1:ncells){
    for (j in 1:nsites[i]){
        Theta[1,1,i,j] <- 1
        Theta[1,2,i,j] <- 0.0000001
        Theta[1,3,i,j] <- 0.0000001
        Theta[2,1,i,j] <- 1-theta2[i,j]
        Theta[2,2,i,j] <- theta2[i,j]
        Theta[2,3,i,j] <- 0.0000001
        Theta[3,1,i,j] <- theta3[1,i,j]                  
        Theta[3,2,i,j] <- theta3[2,i,j]
        Theta[3,3,i,j] <- theta3[3,i,j]
    }
  }

  
  # Define observation matrix (pDet)
  # Order of indices: true state, observed state
  pDet[1,1] <- 1
  pDet[1,2] <- 0.000001
  pDet[1,3] <- 0.000001
  pDet[2,1] <- 1-p2
  pDet[2,2] <- p2
  pDet[2,3] <- 0
  pDet[3,1] <- p31                    # = 1-p32-p33 as per prior section
  pDet[3,2] <- p32
  pDet[3,3] <- p33
  
  # State-space likelihood
  # State equation: model of true states (z)
  for (i in 1:ncells){
    z[i] ~ dcat(Omega[i,])
  }
  
  # Availability equation
  for (i in 1:ncells){
    for (j in 1:nsites[i]){
      a[i,j] ~ dcat(Theta[z[i],,i,j])
    }
  }
  
  # Observation equation
  for (i in 1:ncells){
    for (j in 1:nsites[i]){
      for (k in 1:nsurveys){
        y[i,j,k] ~ dcat(pDet[a[i,j],])
      }
    }
  }
  
  # Derived quantities
  for (i in 1:ncells){
    occ1[i] <- equals(z[i], 1)
    occ2[i] <- equals(z[i], 2)
    occ3[i] <- equals(z[i], 3)
  }
  n.occ[1] <- sum(occ1[]) # Grid cells in state 1
  n.occ[2] <- sum(occ2[]) # Grid cells in state 2
  n.occ[3] <- sum(occ3[]) # Grid cells in state 3
}
")

# Initial values
zst <- rep(3, nrow(batdata$y)) # Initialize at highest possible state
inits <- function(){list(z = zst)}

# Parameters monitored (could add "z")
params <- c("alpha.lpsi", "alpha.lr", "beta.lpsi", "beta.lr",
            "alpha.ltheta2", "alpha.ltheta32", "alpha.ltheta33", 
            "beta.site.type.ltheta2", "beta.site.type.ltheta32", "beta.site.type.ltheta33", 
            "beta.ltheta2", "beta.ltheta32", "beta.ltheta33",
            "pDet", "n.occ")

# MCMC settings
na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 2000 ; nc <- 3

# Call JAGS, check convergence and summarize posteriors
out_M2_msms <- jags(batdata, inits, params, "jags_txt/M2_msms.txt", n.adapt = na,
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
# traceplot(out_M2_msms)
print(out_M2_msms, 3)

diagPlot(out_M2_msms)


### To do on the next steps:
#1. Create a dataframe with the output of these models to compare
# 1.b I want the model algebra in a table together so that I can see exactly what I was comparing
# Which parameters are the ones we are interested in keeping track of right now?
# Because we are using the multiscale component to get a better estimation of the effects of
# environmental covariates by using the site specific covariates, but we are not truly interested in the site
# since the effects and management would probably take place at the cell level anyways

# I personally want to see this model comparison because it doesn't seem that by incorporating these
# the multiscale and such we get such differences
# Paper number 1 is the model comparison and use of these models to generate predictions on alternatives
# for the management

# Paper number 2 is the dynamic model. Again, we can do the null and one that incorporates covariates
# set up the structure for future years of data where we can estimate this correctly. 
# the importance of the dynamic component is that it would allow to incorporate real time pathogen info

# clean this model and write the output on a table. 
# Compare with the deviance
# 






