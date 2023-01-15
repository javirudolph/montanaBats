# TEST SIMULATION AND ESTIMATION FOR SINGLE SEASON MULTISTATE MULTISCALE MODEL

# F. JAVIERA RUDOLPH
# Questions: f.javierarudolph at gmail com

# ISSUE
# I can't fit the model and get errors with JAGS
# There's too much nuisance with the data
# simulate a theoretical dataset and fit model to estimate coeffs


# Libraries
library(tidyverse)
library(sf)
library(jagsUI)

# SIMULATION ----------

## COVARIATES --------
mt_covariates <- readRDS("datafiles/mt_covariates.RDS")

grid_geometry <- mt_covariates %>% 
  select(GRTS_ID)

ggplot(data = grid_geometry) + geom_sf()


# what are we trying to do here?
# Create a prediction for occupancy state given a set of covariates

## CELL LEVEL --------

# intercepts
alpha.lpsi <- c(0, 1.968, 1.719, 1.126, 0.150, 1.134, 0.332)
alpha.lr <- c(0, 0.115, -0.409, 0.321, 0.871, -0.465, 0.855)

# Coefficients
beta.lpsi <- c(0.64, 5.3, 1.2, 4.6, 0.15, 1.054, 0.051, -0.029)
beta.lr <- c(0.566, -0.43, 0.346)


# build the parameters psi and r based on covariates
# center and scale the covariates
center_scale <- function(x){
  center <- mean(x, na.rm = TRUE)
  x <- x - center
  scale <- sd(x, na.rm=TRUE)
  x <- x / scale
  x
}

mt_covariates %>% 
  st_drop_geometry() %>% 
  transmute(GRTS_ID,
            elev = center_scale(DEM_max),
            forest = center_scale(p_forest),
            temp = center_scale(mean_temp),
            precip = center_scale(precip),
            wetlands = center_scale(p_wetland),
            physdiv = center_scale(physio_div),
            karst = karst,
            region = as.numeric(as.factor(eco3_name))) -> scaled_covariates


scaled_covariates %>% 
  # st_drop_geometry() %>% 
  transmute(GRTS_ID,
            lpsi = alpha.lpsi[as.numeric(region)] + beta.lpsi[1] * elev + beta.lpsi[2] * elev^2 + beta.lpsi[3] * temp + beta.lpsi[4] * temp^2 + beta.lpsi[5] * physdiv + beta.lpsi[6] * precip + beta.lpsi[7] * forest + beta.lpsi[8] * wetlands,
            psi = exp(lpsi)/(1+exp(lpsi)),
            lr = alpha.lr[as.numeric(region)] + beta.lr[1] * karst + beta.lr[2] * forest + beta.lr[3] * physdiv,
            r = exp(lr)/(1+exp(lr))) -> psi_and_r

psi <- psi_and_r$psi
r <- psi_and_r$r
Omega <- matrix(c(1-psi, psi*(1-r), psi*r), ncol = 3)

ncells <- nrow(scaled_covariates)
nyears <- 5
z <- array(NA, dim = c(ncells, nyears))

# Calculate state for each cell in year 1
for(i in 1:ncells){
  draw1 <- rmultinom(1,1,Omega[i,])
  z[i,1] <- which(draw1 == 1)
}

# Viz year 1
cbind(grid_geometry, z1 = z[,1]) %>% 
  mutate(z1 = as.factor(z1)) %>% 
  ggplot() +
  geom_sf(aes(fill = z1, color = z1))

# Define transition probability matrix (Phi)
# for now, just using the same matrix for all sites and all years.
gamma1 <- 0.001
gamma2 <- 0.001
phi1 <- 0.5
G <- 0.01
phi2 <- 0.3
D <- 0.5

Phi <- matrix(c(1-gamma1, gamma1*(1-gamma2), gamma1*gamma2,
                1-phi1, phi1*(1-G), phi1*G,
                1-phi2, phi2*D, phi2*(1-D)), nrow = 3, byrow = TRUE)

# fill the next years
for(t in 2:nyears){
  for(i in 1:ncells){
    draw1 <- rmultinom(1,1,Phi[z[i,t-1],])
    z[i,t] <- which(draw1 == 1)
  }
}

# Viz year 1
as.data.frame(z) %>% 
  cbind(grid_geometry, .) %>% 
  pivot_longer(., cols = starts_with("V"), names_to = "year", values_to = "yms") %>% 
  mutate(yms = as.factor(yms)) -> yms_prediction

# viz all
yms_prediction %>% 
  ggplot() +
  facet_wrap(~year) +
  geom_sf(aes(fill = yms, color = yms))

## SITE LEVEL ------
# local occupancy simulation

# Take a year from the previous simulation to have the base


as.data.frame(z) %>% 
  cbind(grid_geometry, .) %>%
  select(GRTS_ID, V2) %>% 
  rename(yms = V2) %>% 
  slice_sample(n = 200) -> sample_grids

sample_grids %>% 
  ggplot() +
  geom_sf(aes(fill = as.factor(yms)))

samp_grts <- pull(sample_grids, GRTS_ID)

scaled_covariates %>% 
  filter(GRTS_ID %in% samp_grts) -> samp_covariates

# Create the observation data frame, with 4 sites per grid cell

# Probability matrix Theta
theta1 <- 0.65
theta2 <- 0.8
theta3 <- 0.7

Theta <- matrix(c(1,0,0,
                  1-theta1, theta1, 0,
                  1-theta2, theta2*(1-theta3), theta2*theta3), nrow = 3, byrow = TRUE)


# create array to fill


# build the for loop to sample for each site with the given theta matrix of probabilities
ni <- length(samp_grts)
nj <- 4

u <- array(NA, dim = c(ni, nj))

for(i in 1:ni){
  for(j in 1:nj){
    draw1 <- rmultinom(1,1,Theta[sample_grids$yms[i],])
    u[i,j] <- which(draw1 == 1)
  }
}

cbind(sample_grids, u) %>% 
  pivot_longer(cols = starts_with("X"), names_to = "site", values_to = "site_ms") %>% 
  ggplot() +
  facet_wrap(~site) +
  geom_sf(aes(fill = as.factor(site_ms)))


## DETECTION ---------

p2 <- 0.9
p3 <- c(0.1, 0.2, 0.7)

detP<- matrix(c(1,0,0,
                 1-p2, p2, 0,
                 p3[1], p3[2], p3[3]), nrow = 3, byrow = TRUE)

nk <- 5
y_obs <- array(NA, dim = c(ni, nj, nk))

for(i in 1:ni){
  for(j in 1:nj){
    for(k in 1:nk){
      draw1 <- rmultinom(1,1,detP[u[i,j],])
      y_obs[i,j,k] <- which(draw1 == 1)
    }
  }
}


# ESTIMATION -------
# now that we have a dataset of simulated observations, for a subset of grid cells
# where each grid cell has four sites surveyed on five occasions each
# we also have the associated covariates for this

# Build the jags model to fit this multiscale, multistate model for one year only

# build a null model first, don't include any covariates
# just make sure structure is well set up


## bundle data ----------------

yms <- array(NA, dim = c(ni, nj, nk),
                          dimnames = list(samp_grts, 1:nj, 1:nk))

array(y_obs, dim = c(ni, nj, nk),
      dimnames = list(samp_grts, 1:nj, 1:nk)) -> test

nregion <- length(unique(samp_covariates$region))

str(testsimData <- list(u = test[,,1], ncells = dim(yms)[1], nsites = dim(yms)[2], nsurveys = dim(yms)[3], 
                        nregion = nregion,
                      region = samp_covariates$region,
                      elev = samp_covariates$elev, 
                      temp = samp_covariates$temp,
                      physdiv = samp_covariates$physdiv,
                      precip = samp_covariates$precip,
                      forest = samp_covariates$forest,
                      wetlands = samp_covariates$wetlands,
                      karst = samp_covariates$karst))


## Model v.1----------------------------
cat(file = "testsim_msms.txt", "
model {

  ### (1) Priors
  ## Omega
  for (i in 1:ncells){
    logit(psi[i]) <- alpha.lpsi[region[i]] + beta.lpsi[1] * elev[i] + beta.lpsi[2] * elev[i]^2 + beta.lpsi[3] * temp[i] + beta.lpsi[4] * temp[i]^2 + beta.lpsi[5] * physdiv[i] + beta.lpsi[6] * precip[i] + beta.lpsi[7] * forest[i] + beta.lpsi[8] * wetlands[i]
    logit(r[i]) <- alpha.lr[region[i]] + beta.lr[1] * karst[i] + beta.lr[2] * forest[i] + beta.lr[3] * physdiv[i]
  }
  
  # Priors for parameters in the linear models of psi and r
  # Region-specific intercepts
  for (k in 1:nregion){
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
  theta1 ~ dunif(0, 1)
  theta2 ~ dunif(0, 1)
  theta3 ~ dunif(0, 1)
  
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector
  for (i in 1:ncells){
    Omega[i,1] <- 1 - psi[i]                 # Prob. of no bats
    Omega[i,2] <- psi[i] * (1-r[i])          # Prob. of occ. by a few bats
    Omega[i,3] <- psi[i] * r[i]              # Prob. of occ. by many bats
  }
  
  # Define local occupancy probability matrix (Theta)
  # Order of indices: global state, local state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta1
  Theta[2,2] <- theta1
  Theta[2,3] <- 0
  Theta[3,1] <- 1-theta2
  Theta[3,2] <- theta2 * (1-theta3)
  Theta[3,3] <- theta2 * theta3
  
  
  ### (3) Likelihood
  
  # global occupancy
  for (i in 1:ncells){
    z[i] ~ dcat(Omega[i,])
  }
  
  # local occupancy
  for (i in 1:ncells){
    for (j in 1:nsites){
      u[i,j] ~ dcat(Theta[z[i], ])
    }
  }
  
}
")


# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(testsimData$ncells))
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha.lpsi", "alpha.lr", "beta.lpsi", "beta.lr", 
            "theta1", "theta2", "theta3")

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; n.iter <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
test_fit <- jags(data = testsimData, inits = inits, parameters.to.save = params,
                 model.file = "testsim_msms.txt", 
                 n.adapt = na, n.chains = nc, n.thin = nt, n.iter = n.iter, 
                 n.burnin = nb, parallel = TRUE)

traceplot(test_fit)
print(test_fit)


## Model v.2 ----------------------------

yms <- array(NA, dim = c(ni, nj, nk),
             dimnames = list(samp_grts, 1:nj, 1:nk))

array(y_obs, dim = c(ni, nj, nk),
      dimnames = list(samp_grts, 1:nj, 1:nk)) -> test

nregion <- length(unique(samp_covariates$region))

str(testsimData <- list(y = test, ncells = dim(yms)[1], nsites = dim(yms)[2], nsurveys = dim(yms)[3], 
                        nregion = nregion,
                        region = samp_covariates$region,
                        elev = samp_covariates$elev, 
                        temp = samp_covariates$temp,
                        physdiv = samp_covariates$physdiv,
                        precip = samp_covariates$precip,
                        forest = samp_covariates$forest,
                        wetlands = samp_covariates$wetlands,
                        karst = samp_covariates$karst))

cat(file = "testsim_msms_det.txt", "
model {
  ### (1) Priors
  ## Omega
  for (i in 1:ncells){
    logit(psi[i]) <- alpha.lpsi[region[i]] + beta.lpsi[1] * elev[i] + beta.lpsi[2] * elev[i]^2 + beta.lpsi[3] * temp[i] + beta.lpsi[4] * temp[i]^2 + beta.lpsi[5] * physdiv[i] + beta.lpsi[6] * precip[i] + beta.lpsi[7] * forest[i] + beta.lpsi[8] * wetlands[i]
    logit(r[i]) <- alpha.lr[region[i]] + beta.lr[1] * karst[i] + beta.lr[2] * forest[i] + beta.lr[3] * physdiv[i]
  }
  
  # Priors for parameters in the linear models of psi and r
  # Region-specific intercepts
  for (k in 1:nregion){
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
  theta1 ~ dunif(0, 1)
  theta2 ~ dunif(0, 1)
  theta3 ~ dunif(0, 1)
  
  ## detP
  # Priors for parameters in observation process
  p2 ~ dunif(0, 1)              # detection with few bats locally
  for (s in 1:3) {                  # detection with many bats locally
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- beta[s] / sum(beta[])
  }
  
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector
  for (i in 1:ncells){
    Omega[i,1] <- 1 - psi[i]                 # Prob. of no bats
    Omega[i,2] <- psi[i] * (1-r[i])          # Prob. of occ. by a few bats
    Omega[i,3] <- psi[i] * r[i]              # Prob. of occ. by many bats
  }
  
  # Define local occupancy probability matrix (Theta)
  # Order of indices: global state, local state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta1
  Theta[2,2] <- theta1
  Theta[2,3] <- 0
  Theta[3,1] <- 1-theta2
  Theta[3,2] <- theta2 * (1-theta3)
  Theta[3,3] <- theta2 * theta3
  
  # Define observation probability matrix (detP)
  # Order of indices: true local state, observed local state
  detP[1,1] <- 1
  detP[1,2] <- 0
  detP[1,3] <- 0
  detP[2,1] <- 1-p2
  detP[2,2] <- p2
  detP[2,3] <- 0
  detP[3,1] <- p3[1]
  detP[3,2] <- p3[2]
  detP[3,3] <- p3[3]
  
  
  ### (3) Likelihood
  
  # global occupancy
  for (i in 1:ncells){
    z[i] ~ dcat(Omega[i,])
  }
  
  # local occupancy
  for (i in 1:ncells){
    for (j in 1:nsites){
      u[i,j] ~ dcat(Theta[z[i], ])
    }
  }
  
  # detection
  for (i in 1:ncells){
    for (j in 1:nsites){
      for (k in 1:nsurveys){
        y[i,j,k] ~ dcat(detP[u[i,j], ])
      }
    }
  }
  
  ## (4) Derived quantities
  # number of cells with each state
  for(i in 1:ncells){
    state1[i] <- equals(z[i],1)
    state2[i] <- equals(z[i],2)
    state3[i] <- equals(z[i],3)
  }
  
    n.occ[1] <- sum(state1)
    n.occ[2] <- sum(state2)
    n.occ[3] <- sum(state3)
  
}
")


# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(testsimData$ncells))
inits <- function(){list(z = zst)}

inits <- function(){
  list()
}

# Parameters monitored
params <- c("alpha.lpsi", "alpha.lr", "beta.lpsi", "beta.lr", 
            "theta1", "theta2", "theta3", "p2", "p3",
            "n.occ")

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; n.iter <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
test_fit_det <- jags(data = testsimData, inits = inits, parameters.to.save = params,
                 model.file = "testsim_msms_det.txt", 
                 n.adapt = na, n.chains = nc, n.thin = nt, n.iter = n.iter, 
                 n.burnin = nb, parallel = TRUE)

traceplot(test_fit_det)
print(test_fit_det)



