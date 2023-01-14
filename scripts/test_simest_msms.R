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
  slice_sample(n = 65) -> sample_grids

sample_grids %>% 
  ggplot() +
  geom_sf(aes(fill = as.factor(yms)))

samp_grts <- pull(sample_grids, GRTS_ID)

scaled_covariates %>% 
  filter(GRTS_ID %in% sampled_grts) -> samp_covariates

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

P_var<- matrix(c(1,0,0,
                 1-p2, p2, 0,
                 p3[1], p3[2], p3[3]), nrow = 3, byrow = TRUE)

nk <- 5
y_obs <- array(NA, dim = c(ni, nj, nk))

for(i in 1:ni){
  for(j in 1:nj){
    for(k in 1:nk){
      draw1 <- rmultinom(1,1,Theta[u[i,j],])
      y_obs[i,j,k] <- which(draw1 == 1)
    }
  }
}


# ESTIMATION -------
# now that we have a dataset of simulated observations, for a subset of grid cells
# where each grid cell has four sites surveyed on five occasions each
# we also have the associated covariates for this

# Build the jags model to fit this multiscale, multistate model for one year only










